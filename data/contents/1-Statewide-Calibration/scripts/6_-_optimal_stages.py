# import necessary libraries 
import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio
import rasterio.mask
import rasterio.warp
import rasterio.features
import rasterio.plot
import rasterio.merge
import rasterio.enums
from tqdm import tqdm

# import the output from script 5 - src_power_law.py
catchnc = gpd.read_file('nc_catchments_rrpct_flow_fema_iparams.shp')

# import NC shapefile and intersect with catchments
nc = gpd.read_file('nc.shp')
catchnc = catchnc.overlay(nc, how= 'intersection').drop(['STATEFP','STATENS','AFFGEOID','GEOID','NAME','LSAD','ALAND','AWATER','STUSPS'], axis = 1)

# import the HAND file for all of NC. This file is too large to include in the repository (11.35 GB), 
# refer to section 4 in the "nc full optimization.ipynb" file to see how to download and merge
# the hand layers needed to create this file in QGIS or another GIS software
hand = rasterio.open('nc_hand.tif')
handr = hand.read()

# read in the FEMA 100-year flood map for NC. This file is zipped in the repository, 
# so it will need to be unzipped first
femanc = gpd.read_file('fema_nc.gpkg')

# read in output from script 2 - src_combiner.py
src = pd.read_csv('nc_src.csv')

# calculate the initial stage value using the rating curve function for each catchment,
# if the initial stage is above the maximum stage of the original src table, set the initial stage to 
# the maximum stage to prevent stages outside the valid range of the original src data
catchnc['stg_init_m'] = catchnc['src_a1'] * catchnc['1pct_Q_m3s'] ** catchnc['src_b1']
catchnc['stg_init_m'] = np.where(catchnc['stg_init_m'] > catchnc['maxstage_m'], catchnc['maxstage_m'], catchnc['stg_init_m'])

# define the agreeance function to evaluate two flood maps covering the same extent
def agreeance(m1, m2):

    # adding the flood maps show where they agree
    add = m1 + m2
    #print(add)
    ww = np.where(add > 1, 1, 0).sum()
    dd = np.where(add == 0, 1, 0).sum()

    # subtracting the flood maps show where they disagree
    sub = m1 - m2
    #print(sub)
    wd = np.where(sub == 1, 1, 0).sum()
    dw = np.where(sub == -1, 1, 0).sum()

    # calculate agreement, underprediction and overprediction
    # for the flood map using the add and sub output from above
    agr = ww / (ww + wd + dw)
    und = dw / (ww + wd + dw)
    ovr = wd / (ww + wd + dw)
    
    #print('agreeance =', np.round(agr * 100, 1), '%')
    #print('wet in 1 dry in 2 =', np.round(und * 100, 1), '%')
    #print('dry in 1 wet in 2 =', np.round(ovr * 100, 1), '%')
    
    out = [agr,und,ovr]
    
    return(out)

# create an empty list to hold the output
out = []

# loop through every catchment and calculate the optimal stage
for catchi in tqdm(range(len(catchnc))):

    # skip catchments with a no_fema flag of 1 since we know
    # they don't overlap with the fema map already from script
    # 4 - fema_flag.py
    if catchnc['no_fema'][catchi] == 1:
        continue

    # determine the COMID of the current catchment
    catch = catchnc['comid'][catchi]

    # extract the geodataframe of the current catchment
    catchgdf = catchnc[catchnc['comid'] == catch]

    # extract HAND for the current catchment
    catchhand, catchtrans = rasterio.mask.mask(hand, catchgdf.geometry, crop=True)

    # replace the no_data value with np.nan
    catchhand = np.where(catchhand < 0, np.nan, catchhand).squeeze()

    # create a copy of HAND without np.nan for math later in the loop
    catchhand2 = np.where(np.isnan(catchhand), 0, catchhand)

    # find the initial stage value for the catchment
    stg_init = catchnc['stg_init_m'][catchi]

    # create a flood map by taking every cell with HAND values less
    # than the stage value
    fldmap = np.where(catchhand <= stg_init, 1, 0).squeeze()

    # clip the fema flood map to the catchment boundary
    femafldmapv = gpd.clip(femanc, catchgdf)

    # check if the clipped fema flood map is empty, skip this catchment if so
    if len(femafldmapv) == 0:
        catchnc.loc[catchi,'no_fema'] = 1
        continue

    # rasterize the fema flood map for the catchment and ensure it is the same shape 
    # as the HAND layer for the catchment
    geom_value = ((geom,value) for geom, value in zip(femafldmapv.geometry, [1]*len(femafldmapv.geometry)))
    femafldmap = rasterio.features.rasterize(geom_value,
                                             out_shape= catchhand.shape,
                                             transform= catchtrans,
                                             all_touched= False,
                                             fill = 0,
                                             merge_alg= rasterio.enums.MergeAlg.replace,
                                             dtype= rasterio.uint8)

    # check if the rasterized fema flood map is all zeros, if so it is empty
    # and the catchment is skipped
    if femafldmap.sum() == 0:
        catchnc.loc[catchi,'no_fema'] = 1
        continue

    # calculate the initial agreement between the flood maps
    agree_init = agreeance(fldmap, femafldmap) 

    # use catchhand2 which does not have np.nan to check the max HAND value for the catchment.
    # some catchments have extremely high maximum HAND values, so for these catchments the maximum
    # HAND value is reduced to 150 m, otherwise the next loop will run extremely long. This is only
    # used for a small number of catchments 
    if catchhand2.max() > 150:
        handmax = 150
    else:
        handmax = catchhand2.max()

    # create a sequence of stages from 0 to the maximum stage value for the catchment
    stgs = np.arange(0, handmax + 0.025, 0.025)

    # create an empty list to store the agreement values 
    ags = []

    # loop through the stage values, calculate a flood map for each stage, calculate
    # agreement for each stage, save the agreement values and stages to the list created 
    for stg in stgs:
        stg_fldmap = np.where(catchhand <= stg, 1, 0).squeeze() 
        ag = agreeance(stg_fldmap, femafldmap)
        ag.append(stg)
        ags.append(ag)

    # find the optimal stage value associated with the maximum agreement value
    ags = pd.DataFrame(ags, columns= ['agreeance','wet1dry2','dry1wet2', 'stage (m)'])
    ags_opt = ags[ags['agreeance'] == np.max(ags['agreeance'])].reset_index()
    ags_opt = ags_opt['stage (m)']
    if len(ags_opt) != 1:
        stg_opt = float(ags_opt[0])
    else:
        stg_opt = float(ags_opt)
    agmax = ags['agreeance'].max()

    # add output to the out list
    out.append([catch, stg_init, agree_init[0], stg_opt, agmax])

# convert output to a dataframe
outdf = pd.DataFrame(out, columns=['catchmask', 'stg_init_m', 'agree_init', 'stg_opt_m', 'agree_opt'])

# This creates a figure that shows the histogram of the agreement values and the optimal stage values across all catchments
fig, (ax1, ax2) = plt.subplots(1,2,figsize = (10,5))

sns.histplot(data = outdf['agree_opt'], ax = ax1, color = 'maroon', kde = True)
ax1.set_title('Optimal agreement across catchments')
ax1.set_xlabel('HAND vs. FEMA agreement')

sns.histplot(data = outdf['stg_opt_m'], ax = ax2, color = 'olive', kde = True)
ax2.set_title('Optimal stage across catchments')
ax2.set_xlabel('Optimal stage (m)')

fig.tight_layout()

# save the output to csv
outdf.to_csv('statewide_opt_output.csv')
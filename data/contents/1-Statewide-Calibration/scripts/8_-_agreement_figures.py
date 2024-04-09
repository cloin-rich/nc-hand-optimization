# import necessary libraries
import geopandas as gpd
import rasterio
import rasterio.mask
import rasterio.warp
import rasterio.features
import rasterio.plot
import rasterio.merge
import rasterio.enums
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import seaborn as sns

# import the output from script 5 - src_power_law.py
catchnc = gpd.read_file('nc_catchments_rrpct_flow_fema_iparams.shp')

# import NC shapefile and intersect with catchments
nc = gpd.read_file('nc.shp')
catchnc = catchnc.overlay(nc, how= 'intersection').drop(['STATEFP','STATENS','AFFGEOID','GEOID','NAME','LSAD','ALAND','AWATER','STUSPS'], axis = 1)

# read in output from script 7 - scale_factors.py
out = pd.read_csv('statewide_opt_output_final.csv')

# import the HAND file for all of NC. This file is too large to include in the repository (11.35 GB), 
# refer to section 4 in the "nc full optimization.ipynb" file to see how to download and merge
# the hand layers needed to create this file in QGIS or another GIS software
hand = rasterio.open('nc_hand.tif')
handr = hand.read()

# read in the FEMA 100-year flood map for NC. This file is zipped in the repository, 
# so it will need to be unzipped first
femanc = gpd.read_file('fema_nc.gpkg')
femanc = gpd.clip(femanc, nc)


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

# define a raster legend function for the final comparison map
def raster_legend(raster,labels,colmap):
    values = np.unique(raster)
    fig = plt.figure()
    im = plt.imshow(comp, cmap= colmap)
    colors = [im.cmap(im.norm(value)) for value in values]
    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(values))]
    plt.close(fig)
    return(patches)

# determine which catchments were calibrated and which were not. This is 
# used to filter out catchments that were not calibrated from the final agreement 
# calculation
calcatch = catchnc.merge(out, how='left', on= 'comid', indicator = True)
nocalcatch = calcatch[calcatch['_merge'] == 'left_only']
calcatch = calcatch[calcatch['_merge'] == 'both']

# replace the no_data value in the HAND layer with np.nan
handr = np.where(handr <0, np.nan, handr).squeeze()

# rasterize the calibrated catchments with optimal stage values as the cell values
geom_value = ((geom,value) for geom, value in zip(calcatch.geometry, calcatch['stg_opt_m']))

calstg = rasterio.features.rasterize(geom_value,
                                     out_shape = hand.shape,
                                     transform = hand.transform,
                                     all_touched = False,
                                     fill = -9999,   # background value
                                     merge_alg = rasterio.enums.MergeAlg.replace,
                                     dtype = rasterio.float32)

# rasterize the uncalibrated catchments with a value of 1 as the cell values
geom_value2 = ((geom,value) for geom, value in zip(nocalcatch.geometry, [1]*len(nocalcatch.geometry)))

nocal = rasterio.features.rasterize(geom_value2,
                                    out_shape= hand.shape,
                                    transform= hand.transform,
                                    all_touched= False,
                                    fill = 0,
                                    merge_alg= rasterio.enums.MergeAlg.replace,
                                    dtype= rasterio.uint8)

# calculate the calibrated HAND flood map by finding anywhere HAND is less than the 
# optimal stage values from each catchment
handfldmap = calstg-handr
handfldmap = np.where(np.isnan(handfldmap), 0, handfldmap)
handfldmap = np.where(handfldmap <= 0, 0, 1)

# rasterize the FEMA flood map 
geom_value3 = ((geom,value) for geom, value in zip(femanc.geometry, [1]*len(femanc.geometry)))
femafldmap = rasterio.features.rasterize(geom_value3,
                                         out_shape= hand.shape,
                                         transform= hand.transform,
                                         all_touched= False,
                                         fill = 0,
                                         merge_alg= rasterio.enums.MergeAlg.replace,
                                         dtype= rasterio.uint8)

# turn off parts of the FEMA floodmap that are present in catchments that 
# were not calibrated. This ensures only catchments that were calibrated are
# considered in the final agreement map
femafldmapfilt = femafldmap - nocal
femafldmapfilt = np.where(femafldmapfilt != 1, 0, femafldmapfilt)

# print the agreement statistics for the calibrated HAND map
print(agreeance(handfldmap, femafldmapfilt))

# plot the statewide HAND vs. FEMA agreement map showing all cells where the maps
# agree and where they dont
fig,ax = plt.subplots(figsize = (10,4))

add = handfldmap + femafldmapfilt
add = np.where(add == 1, 0, add)
sub = handfldmap - femafldmapfilt
comp = add + sub

cmap = colors.ListedColormap(['#872325','#ededed','#C7BD0E','#84bccb'])

rasterio.plot.show(comp, ax = ax, transform=hand.transform, cmap = cmap)

ax.set_ylabel('Latitude ($^\circ$)')
ax.set_xlabel('Longitude ($^\circ$)')
ax.legend(handles = raster_legend(comp, ['HAND Not Flooded / FEMA Flooded', 'Both Not Flooded', 'HAND Flooded / FEMA Not Flooded', 'Both Flooded'], cmap),
          loc = 'lower left')

# plt.rcParams['savefig.dpi']=450
# plt.savefig('statewide_AOU_map.tif')

# create histogram figures to show agreement and scale factors across all catchments
fig, (ax1, ax2) = plt.subplots(1,2,figsize = (10,5))

sns.histplot(data = out['agree_opt'], ax = ax1, color = 'maroon', kde = False)
# ax1.set_title('Optimal agreement across catchments')
ax1.set_xlabel('HAND vs. FEMA Agreement')

# only plot the bottom 95% of scale factors as extreme outliers make this plot 
# illegible otherwise
scalequant = out['scale_factor'].quantile(0.95)

sns.histplot(data = out[out['scale_factor'] < scalequant]['scale_factor'], ax = ax2, color = 'olive', kde = False)
# ax2.set_title('Optimal stage across catchments')
ax2.set_xlabel('Optimal Scale Factor')

fig.tight_layout()

# plt.rcParams['savefig.dpi']=450
# plt.savefig('statewide_agree_scalefac_hists.tif')
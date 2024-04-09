# import libraries 
import geopandas as gpd
import pandas as pd
from tqdm import tqdm
import simpledbf

# read in regression region & catchment shapefiles 
regreg = gpd.read_file('nc_regression_regions.shp')
catchnc = gpd.read_file('nc_catchments.shp')

# create new columns for regression region areas
catchnc['A_r1_sqkm'] = 0.0
catchnc['A_r2_sqkm'] = 0.0
catchnc['A_r3_sqkm'] = 0.0
catchnc['A_r4_sqkm'] = 0.0
catchnc['A_r5_sqkm'] = 0.0
catchnc['Atot_sqkm'] = 0.0

# determine how much area of each regression region falls in each catchment
for catchi in tqdm(range(len(catchnc['comid']))):
    # find the COMID for the current catchment
    catch = catchnc['comid'][catchi]
    
    # get the geodataframe for the current catchment
    catchgdf = catchnc[catchnc['comid']==catch]
    
    # calculate the total area of the catchment and save
    tot_area = float(catchgdf.to_crs(5070).area * 1e-6)
    catchnc['Atot_sqkm'][catchi] = tot_area

    # clip the regression regions to the catchment, skip catchment if the clip returns empty
    catchclip = gpd.clip(regreg, catchgdf).reset_index()
    if len(catchclip) == 0:
        continue

    # project catchment geodataframe and calculate areas of each regression region within the catchment
    catchclip = catchclip.to_crs(5070)
    areas = catchclip.area * 1e-6

    # if only one regression region in the catchment, save the area to the output dataframe
    if len(catchclip) == 1:
        if catchclip['FF_REGION'][0] == 'Region 1':
            catchnc['A_r1_sqkm'][catchi] = tot_area
        if catchclip['FF_REGION'][0] == 'Region 2':
            catchnc['A_r2_sqkm'][catchi] = tot_area
        if catchclip['FF_REGION'][0] == 'Region 3':
            catchnc['A_r3_sqkm'][catchi] = tot_area
        if catchclip['FF_REGION'][0] == 'Region 4':
            catchnc['A_r4_sqkm'][catchi] = tot_area
        if catchclip['FF_REGION'][0] == 'Region 5':
            catchnc['A_r5_sqkm'][catchi] = tot_area

    # if more than one regression region in the catchment, save the area of each region to the output dataframe
    else:
        catchclip['area_sqkm'] = areas
        for regi in range(len(catchclip['FF_REGION'])):
            reg = catchclip['FF_REGION'][regi]
            if reg == 'Region 1':
                catchnc['A_r1_sqkm'][catchi] = float(catchclip[catchclip['FF_REGION'] == reg]['area_sqkm'])
            if reg == 'Region 2':
                catchnc['A_r2_sqkm'][catchi] = float(catchclip[catchclip['FF_REGION'] == reg]['area_sqkm'])
            if reg == 'Region 3':
                catchnc['A_r3_sqkm'][catchi] = float(catchclip[catchclip['FF_REGION'] == reg]['area_sqkm'])
            if reg == 'Region 4':
                catchnc['A_r4_sqkm'][catchi] = float(catchclip[catchclip['FF_REGION'] == reg]['area_sqkm'])
            if reg == 'Region 5':
                catchnc['A_r5_sqkm'][catchi] = float(catchclip[catchclip['FF_REGION'] == reg]['area_sqkm'])
                
# read in the flow line data for NC, extract the hydrosequence data which lays out the flow network
flows = simpledbf.Dbf5('nc_flows.dbf').to_dataframe()
flows = flows[['COMID','Hydroseq', 'UpHydroseq', 'DnHydroseq', 'StreamOrde','TotDASqKM']]

# set up columns for each catchment that will be updated in the for loop below
# if UpHydroSeq is 0, it is a headwater catchment, so the total upstream regression 
# region area will just be the local regression region area. Set up the "used" flag
# to indicate if a catchment has been cylced through already or not
joined = catchnc.merge(flows, left_on = 'comid', right_on = 'COMID', validate = '1:1')
joined['used'] = np.where(joined['UpHydroseq'] == 0, True, False)
joined['A_r1_tot'] = np.where(joined['UpHydroseq'] == 0, joined['A_r1_sqkm'], 0.0)
joined['A_r2_tot'] = np.where(joined['UpHydroseq'] == 0, joined['A_r2_sqkm'], 0.0)
joined['A_r3_tot'] = np.where(joined['UpHydroseq'] == 0, joined['A_r3_sqkm'], 0.0)
joined['A_r4_tot'] = np.where(joined['UpHydroseq'] == 0, joined['A_r4_sqkm'], 0.0)
joined['A_r5_tot'] = np.where(joined['UpHydroseq'] == 0, joined['A_r5_sqkm'], 0.0)

# sort the geodataframe by descending UpHydroSeq, allowing us to cycle through the catchments 
# from upstream to downstream
joined = joined.sort_values(by = 'UpHydroseq', ascending=False).reset_index()

# for each catchment, update the regression region area totals for that catchment by adding the
# value from the upstream catchment to the current catchment. Update the "used" flag to make sure that
# all catchments have been used to calculate the regression region area totals by the end of the loop
for catchi in tqdm(range(len(joined))):
    preind = joined[joined['Hydroseq'] == joined['UpHydroseq'][i]].index
    joined.loc[catchi,'A_r1_tot'] = float(joined.loc[preind,'A_r1_tot']) + float(joined.loc[catchi,'A_r1_sqkm'])
    joined.loc[catchi,'A_r2_tot'] = float(joined.loc[preind,'A_r2_tot']) + float(joined.loc[catchi,'A_r2_sqkm'])
    joined.loc[catchi,'A_r3_tot'] = float(joined.loc[preind,'A_r3_tot']) + float(joined.loc[catchi,'A_r3_sqkm'])
    joined.loc[catchi,'A_r4_tot'] = float(joined.loc[preind,'A_r4_tot']) + float(joined.loc[catchi,'A_r4_sqkm'])
    joined.loc[catchi,'A_r5_tot'] = float(joined.loc[preind,'A_r5_tot']) + float(joined.loc[catchi,'A_r5_sqkm'])
    joined.loc[catchi,'used'] = True

# Calculate the total upstream regression region percentages using the totals for each catchment
joined['A_r_sum'] = joined['A_r1_tot'] + joined['A_r2_tot'] + joined['A_r3_tot'] + joined['A_r4_tot'] + joined['A_r5_tot'] 
joined['r1_%'] = np.where(joined['A_r_sum'] != 0, joined['A_r1_tot']/joined['A_r_sum']*100, 0.0)
joined['r2_%'] = np.where(joined['A_r_sum'] != 0, joined['A_r2_tot']/joined['A_r_sum']*100, 0.0)
joined['r3_%'] = np.where(joined['A_r_sum'] != 0, joined['A_r3_tot']/joined['A_r_sum']*100, 0.0)
joined['r4_%'] = np.where(joined['A_r_sum'] != 0, joined['A_r4_tot']/joined['A_r_sum']*100, 0.0)
joined['r5_%'] = np.where(joined['A_r_sum'] != 0, joined['A_r5_tot']/joined['A_r_sum']*100, 0.0)

# save the output to shapefile, this shapefile serves as an input to script 3 - cath_DA_flow.py
out = joined[['comid', 'r1_%', 'r2_%', 'r3_%', 'r4_%', 'r5_%', 'TotDASqKM', 'geometry']]
out.to_file('nc_catchments_rrpct.shp')
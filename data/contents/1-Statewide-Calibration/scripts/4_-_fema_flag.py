# import libraries
import geopandas as gpd
import numpy as np

# read in the output from script 3 - catch_DA_flow.py
catchnc = gpd.read_file('nc_catchments_rrpct_flow.shp')

# read in the FEMA 100-year flood map for NC. This file is zipped in the repository, 
# so it will need to be unzipped first
fema = gpd.read_file('fema_nc.gpkg')

# determine which catchments overlap with the FEMA flood maps 
catchncfema = catchnc.overlay(fema, how='intersection')
femacatch = catchncfema['comid'].unique()
femaind = np.where(catchnc.comid.isin(femacatch) == False)

# label catchments overlapping FEMA with 1, others as 0
catchnc['no_fema'] = 0
catchnc['no_fema'][femaind] = 1

# save the output to shapefile. This file serves as the input to script 5 - src_power_lay.py
catchnc.to_file('nc_catchments_rrpct_flow_fema.shp')
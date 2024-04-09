# import libraries
import geopandas as gpd

# read in output from script 1 - regression_region_percentages.py
catch = gpd.read_file('nc_catchments_rrpct.shp')

# calculate the 100-year flow according to the regional regression equations
catch['TotDASqMi'] = catch['TotDASqKM'] * 0.386102
catch['1pct_Q_m3s'] = 0.028316847*(10**(2.64 + 0.00218*catch['r1_%'] - 0.00311*catch['r3_%'] + 0.00309*catch['r5_%']))*(catch['TotDASqMi']**(0.605 + 0.00161*catch['r2_%']))

# save the output to shapefile, this serves as an input to script 4 - fema_flag.py
catch.to_file('/Users/colin/Desktop/nc_catchments_rrpct_flow.shp')
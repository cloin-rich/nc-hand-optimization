# import necessary libraries 
import geopandas as gpd
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from tqdm import tqdm

# import the output from script 5 - src_power_law.py
catchnc = gpd.read_file('nc_catchments_rrpct_flow_fema_iparams.shp')

# import NC shapefile and intersect with catchments
nc = gpd.read_file('nc.shp')
catchnc = catchnc.overlay(nc, how= 'intersection').drop(['STATEFP','STATENS','AFFGEOID','GEOID','NAME','LSAD','ALAND','AWATER','STUSPS'], axis = 1)

# read in output from script 2 - src_combiner.py
src = pd.read_csv('nc_src.csv')

# read in output from script 6 - optimal_stages.py
out = pd.read_csv('statewide_opt_output.csv')

# join the out csv with the catchments geodataframe
stgoptgdf = catchnc.merge(out, how = 'right', left_on ='comid', right_on = 'catchmask').drop(['Unnamed: 0', 'catchmask'], axis = 1)

# create columns for the output dataframe
stgoptgdf['src_a2'] = 0.0
stgoptgdf['src_b2'] = 0.0
stgoptgdf['src_R2_2'] = 0.0

# define the general form of the power law function used to relate stage, wetted area (A), 
# and hydraulic radius (R)
def f(x, a, b):
    return a*x**b

# for each catchment, fit the relationship between stage and A*R^(2/3) using curve_fit(),
# calculate the R squared value of the relationship, save the a and b coefficient values
# and the R squared value for each catchment
for catchi in tqdm(range(len(stgoptgdf['comid']))):
    
    # find the COMID of the current catchment
    catch = stgoptgdf['comid'][catchi]

    # extract the SRC for the current catchment
    csrc = src[src['CatchId'] == catch]

    # subset the SRC to include rows only just higher than the optimal stage value 
    # determined in script 6 - optimal_stages.py
    stgsneed = csrc[csrc['Stage'] < (float(stgoptgdf[stgoptgdf['comid']==catch]['stg_opt_m'])+2)].reset_index()
    stgsneed = stgsneed.drop(np.where(np.isnan(stgsneed['WetArea (m2)']*stgsneed['HydraulicRadius (m)']**(2/3)))[0].tolist())

    # some src tables have invalid A and R data equal to zero for all stages, in these cases skip these catchments
    if sum(stgsneed['WetArea (m2)']*stgsneed['HydraulicRadius (m)']**(2/3) == 0) > 10:
        continue

    # fit the relationship between stage and A*R^(2/3) for each catchment using curve_fit()
    popt, pcov = curve_fit(f, stgsneed['Stage'],(stgsneed['WetArea (m2)']*stgsneed['HydraulicRadius (m)']**(2/3)))

    # save the coefficient values and the R squared values for each catchment
    stgoptgdf.loc[catchi,'src_a2'] = popt[0]
    stgoptgdf.loc[catchi,'src_b2'] = popt[1]
    r2_2 = r2_score((stgsneed['WetArea (m2)']*stgsneed['HydraulicRadius (m)']**(2/3)), f(stgsneed['Stage'], *popt))
    stgoptgdf.loc[catchi,'src_R2_2'] = r2_2

# filter out the handful of catchments that had invalid SRC data
stgoptgdf = stgoptgdf.drop(np.where(stgoptgdf['src_R2_2'] == 0)[0]).reset_index(drop=True)

# calculate the calibrated manning's n value for each catchment
stgoptgdf['mannings_n_opt'] = (stgoptgdf['src_a2']*stgoptgdf['stg_opt_m']**stgoptgdf['src_b2'])/(stgoptgdf['src_a2']*stgoptgdf['stg_init_m']**stgoptgdf['src_b2'])*0.05

# calculate the rating curve scale factor for each catchment 
stgoptgdf['scale_factor'] = (stgoptgdf['mannings_n_opt']/0.05)**stgoptgdf['src_b1']

# calculate the calibrated a coefficient for each catchment using the scale_factor
stgoptgdf['src_a1_opt'] = stgoptgdf['src_a1']*stgoptgdf['scale_factor']

# drop the catchment geometry and save the output to csv. Geometry can interfere with
# the csv export so it is dropped. This serves as the input to script 8 - agreement_figures.py
stgoptgdf.drop('geometry', axis = 1).to_csv('statewide_opt_output_final.csv')
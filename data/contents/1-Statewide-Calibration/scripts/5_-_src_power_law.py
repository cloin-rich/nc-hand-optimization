# import libraries
import geopandas as gpd
from scipy.optimize import curve_fit
from tqdm import tqdm

# read in output from script 4 - fema_flag.py
catchnc = gpd.read_file('nc_catchments_rrpct_flow_fema.shp')

# read in output from script 2 - src_combiner.py
src = pd.read_csv('nc_src.csv')

# create columns for output from for loop
catchnc['src_a1'] = 0.0
catchnc['src_b1'] = 0.0
catchnc['src_R2_1'] = 0.0
catchnc['maxstage_m'] = 0.0

# define structure of rating curve power law functions
def f(x,a,b):
    return a*x**b

# for each catchment, extract the corresponding SRC table, determine the a and b
# coefficient values using curve_fit(), calculate the R squared value for the rating 
# curve function for each catchment, save the coefficient values and the R squared values 
# to the output dataframe for each catchment. Also save the max stage value from each src
# to cap the stages produced by each src in case they generate stages above the max
for catchi in tqdm(range(len(catchnc['comid']))):
    catch = catchnc['comid'][catchi]
    csrc = src[src['CatchId'] == catch]
    csrc['Discharge (m3s-1)'] = np.where(np.isnan(csrc['Discharge (m3s-1)']), 0, csrc['Discharge (m3s-1)'])
    popt, pcov = curve_fit(f, csrc['Discharge (m3s-1)'], csrc['Stage'])
    catchnc.loc[catchi,'src_a1'] = popt[0]
    catchnc.loc[catchi,'src_b1'] = popt[1]
    catchnc.loc[catchi,'maxstage_m'] = csrc['Stage'].max()
    src_R2_1 = r2_score(csrc['Stage'], f(csrc['Discharge (m3s-1)'], *popt))
    catchnc.loc[catchi,'src_R2_1'] = src_R2_1

# save the output to shapefile, this is the input to script 6 - optimal_stages.py
catchnc.to_file('nc_catchments_rrpct_flow_fema_iparams.shp')
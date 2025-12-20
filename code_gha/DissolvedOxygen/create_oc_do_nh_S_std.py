import re
import xarray as xr
import numpy as np
import pandas as pd
from fun_pd_df2csvR_time import fun_pd_df2csvR_time



# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# ie year
iea_yr = 2025

# --station wanted, two stations newport hydrographic at 5 and 25 km
sttn_wnt = ['NH05', 'NH25']
num_sttn = len(sttn_wnt)

# --input directory
dir_in = './data_x13/DissolvedOxygen/'

# --input file type extension
file_type_in = '.nc'
file_pre_in = 'do_ts_'

# --variable name
var_name = 'Dissolved Oxygen'

# --depth want
depth_wnt = np.float64([50, 150])

# --variable xr.ds name
var_pre_wnt = 'Oxygen'

# --name of useful attributes
att_wnt = ['station', 'lat', 'lon']
num_att = len(att_wnt)

# --IEA file names
file_pre = 'oc_do_Newport'

# --IEA labels
y_lbl = 'DO (ml/L)'
ttl_lbl = var_name

# --names of time, depth, lat dimensions
var_time = 'time'

# --output CSV file
dir_out = './data_x13/csv_files/'
yr_csv_bgn = 1998
yr_csv_end = iea_yr
# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

for i in range(0, num_sttn):
    # --depth
    depth = depth_wnt[i]
    var_wnt = var_pre_wnt + depth.astype('int').astype('str') + 'm'

    # --create input filename
    fn_in = dir_in + file_pre_in + sttn_wnt[i] + file_type_in

    # --open the netcdf files as an xarray
    dall = xr.open_dataset(fn_in)

    # --get the attributes
    station = dall.attrs[att_wnt[0]]
    lat = dall.attrs[att_wnt[1]]
    lon = dall.attrs[att_wnt[2]]

    # --station label
    line_sttn_lbl = station

    # --create location label for the IEA title
    loc_lbl = ' ({:4.2f}N {:4.2f}W)'.format(lat, lon)

    # --variable as xr.da and pd.df
    dai = dall[var_wnt]
    dfi = dall[var_wnt].to_dataframe()

    # --resample by taking 'M' means, i.e. monthly means
    dfMmn = dfi.resample('ME').mean()

    # --resample by taking 'Q' means, i.e. monthly means
    dfQmn = dfMmn.resample('QE-MAR').mean()
    dfQsd = dfMmn.resample('QE-MAR').std(ddof=0)

    # --combine mean and sd into one pd.df
    dfQ = dfQmn.rename(index=str, columns={var_wnt: "data"})
    dfQ['Datetime'] = pd.to_datetime(dfQ.index)
    dfQ = dfQ.set_index('Datetime')

    # dfQ['sd'] = dfQsd.values
    dfQ['sd'] = dfQsd.values

    in1 = dfQ.index.month == 3
    in2 = dfQ.index.month == 6
    in3 = dfQ.index.month == 9
    in4 = dfQ.index.month == 12
    dfQ1 = dfQ[in1]
    dfQ2 = dfQ[in2]
    dfQ3 = dfQ[in3]
    dfQ4 = dfQ[in4]
    ts_lbl1 = 'Winter ' + ttl_lbl + \
        ' at {}'.format(depth.astype(int)) + ' m: ' + line_sttn_lbl + loc_lbl
    ts_lbl2 = 'Spring ' + ttl_lbl + \
        ' at {}'.format(depth.astype(int)) + ' m: ' + line_sttn_lbl + loc_lbl
    ts_lbl3 = 'Summer ' + ttl_lbl + \
        ' at {}'.format(depth.astype(int)) + ' m: ' + line_sttn_lbl + loc_lbl
    ts_lbl4 = 'Fall ' + ttl_lbl + \
        ' at {}'.format(depth.astype(int)) + ' m: ' + line_sttn_lbl + loc_lbl

    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfQ, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list = [dfQ1, dfQ2, dfQ3, dfQ4]
    num_order = len(df_list)

    fn_out = file_pre + '_' + re.sub(r'[^\w]', '', sttn_wnt[i]) + '_S.csv'

    metric_lbl = y_lbl
    ts_lbl = [ts_lbl1, ts_lbl2, ts_lbl3, ts_lbl4]
    clmns_iea = ['year', 'time', 'index', 'error', 'SElo', 'SEup',
                 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

    fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_list, lat, lon, depth,
                                metric_lbl, ts_lbl, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

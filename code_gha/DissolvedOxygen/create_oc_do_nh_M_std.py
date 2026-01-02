import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib import interactive
from fun_pd_df2csvR_time import fun_pd_df2csvR_time
interactive(True)


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
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

# --dpth want
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
dir_out = './csv_for_erddap/'
yr_csv_bgn = 1998

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

df_list = []
ts_lbl = []
depth_list = []
lon_list = []
lat_list = []
for i in range(0, num_sttn):
    # --depth
    depth = depth_wnt[i]
    depth_list.append(depth)

    # var want at depth
    var_wnt = var_pre_wnt + depth.astype('int').astype('str') + 'm'

    # --create input filename
    fn_in = dir_in + file_pre_in + sttn_wnt[i] + file_type_in

    # --open the netcdf files as an xarray
    dall = xr.open_dataset(fn_in)

    # --get the last year of the dataset
    yr_csv_end = dall.time.dt.year.data[-1]

    # --get the attributes
    station = dall.attrs[att_wnt[0]]
    lat = dall.attrs[att_wnt[1]]
    lon = dall.attrs[att_wnt[2]]
    lon_list.append(lon)
    lat_list.append(lat)

    # --station label
    line_sttn_lbl = station

    # --create location label for the IEA title
    loc_lbl = ' ({:4.2f}N {:4.2f}W)'.format(lat, lon)

    # --variable as a xr.da and pd.df
    dai = dall[var_wnt]
    dfi = dall[var_wnt].to_dataframe()

    # --resample by taking 'M' means, i.e. monthly means
    dfMmn = dfi.resample('ME').mean()

    # --get std from xr.da not pd.df, think I'm doing something wrong with std and nan and pandas
    daMsd = dai.resample(time='ME').std(skipna=1)

    # --combine mean and sd into one pd.df
    dfM = dfMmn.rename(index=str, columns={var_wnt: "data"})
    dfM['Datetime'] = pd.to_datetime(dfM.index)
    dfM = dfM.set_index('Datetime')

    # dfM['sd'] = dfMsd.values
    dfM['sd'] = daMsd.values

    # title
    ts_lbl1 = 'Monthly ' + ttl_lbl + \
        ' at {}'.format(depth.astype(int)) + ' m: ' + line_sttn_lbl + loc_lbl
    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfM, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list.append(dfM)
    ts_lbl.append(ts_lbl1)

# num_order = len(df_list)
# fn_out = file_pre + '_' + re.sub(r'[^\w]', '', sttn_wnt[i]) + '.csv'

# filename is based on Lynn's Uploader tool name
fn_out = '{}_M.csv'.format(file_pre)

metric_lbl = y_lbl

clmns_iea = ['year', 'month', 'time', 'index', 'error', 'SElo', 'SEup',
             'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_list, lat_list, lon_list, depth_list,
                                 metric_lbl, ts_lbl, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

import itertools
import os
import xarray as xr
import pandas as pd
import numpy as np
from fun_pd_df2csvR_time import fun_pd_df2csvR_time


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# --UI lat wanted
lat_wnt = [45, 39, 33]

# input filename
dir_out = './data_gha/bakunUI/'
dir_in = dir_out
fn_in = '{}UI_monthly.nc'.format(dir_in)

# --input file type extension
file_type_in = '.nc'

# --variable xr.ds name
var_wnt = 'ui_mon'

# --attrributes xr.ds name
att_wnt = ['title']

# --IEA file names
file_pre = 'oc_ui'

# --IEA file columns
y_lbl = 'UI (m^3/s/100 m coastline)'
lat = lat_wnt
depth = [np.nan, np.nan, np.nan]

# --names of time, depth, lat dimensions
var_time = 'time'

# --output CSV file
dir_out = './csv_for_erddap/'
yr_csv_bgn = 1967


# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# size of input variables
num_lat = len(lat_wnt)

# open monthly UI netcdf
ds1 = xr.open_dataset(fn_in)

df_all_list = []
ts_lbl_list = []
metric_list = []
lon = []
for i in range(0, num_lat):
    # open the netcdf files as an xarray
    dall = xr.open_dataset(fn_in).sel(lat=lat_wnt[i])

    # get longitude
    lon.append(int(dall['lon'].data))
    
    # --get the attributes
    title = dall.attrs[att_wnt[0]]

    # xr.dataset as pd.df
    dfi = dall[var_wnt].to_dataframe()

    # remove the lat column
    dfi = dfi.drop('lat', axis=1)

    # the last year of the csv file
    yr_csv_end = dfi.index.year[-1]

    # resample by taking 'M' means, i.e. monthly means
    dfMmn = dfi.resample('ME').mean()
    dfMsd = dfi.resample('ME').std()

    # combine mean and sd into one pd.df
    dfM = dfMmn.rename(index=str, columns={var_wnt: "data"})
    dfM['Datetime'] = pd.to_datetime(dfM.index)
    dfM = dfM.set_index('Datetime')
    dfM['sd'] = dfMsd.values

    ts_lbl1 = 'Monthly ' + title + ' ({} °N)'.format(int(lat_wnt[i]))
    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfM, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list = [dfM]
    num_order = len(df_list)

    fn_out = '{}_M.csv'.format(file_pre)

    metric_lbl = y_lbl
    ts_lbl = [ts_lbl1]
    clmns_iea = ['year', 'month', 'time', 'index', 'error', 'SElo', 'SEup',
                 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

    df_all_list.append(df_list)
    ts_lbl_list.append(ts_lbl)
    metric_list.append([str(metric_lbl)])

df_flat = list(itertools.chain(*df_all_list))
ts_lbl_flat = list(itertools.chain(*ts_lbl_list))
metric_flat = list(itertools.chain(*metric_list))

fn_out = '{}_M.csv'.format(file_pre)
fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_flat, lat, lon, depth, metric_flat, ts_lbl_flat, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

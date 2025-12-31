import os
import itertools
import xarray as xr
import pandas as pd
import numpy as np
from fun_pd_df2csvR_time import fun_pd_df2csvR_time
# -*- coding: utf-8 -*-
# pylint: disable=C0103


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------

# --UI lat wanted
lat_wnt = [33, 39, 45]

# --input directory
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
depth = np.nan

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
lon_list = []
lat_list = []
depth_list = []
for i in range(0, num_lat):
    # open the netcdf files as an xarray
    dall = xr.open_dataset(fn_in).sel(lat=lat_wnt[i])

    # get lat and lon
    lon1 = int(dall['lon'].data)
    lat1 = int(dall['lat'].data)

    # --get the attributes
    title = dall.attrs[att_wnt[0]]

    # --xr.dataset as pd.df
    dfi = dall[var_wnt].to_dataframe()

    # remove the lat column
    dfi = dfi.drop('lat', axis=1)

    # the last year of the csv file
    yr_csv_end = dfi.index.year[-1]

    # --resample by taking 'M' means, i.e. monthly means
    dfQmn = dfi.resample('QE-MAR').mean()
    dfQsd = dfi.resample('QE-MAR').std()

    # --combine mean and sd into one pd.df
    dfQ = dfQmn.rename(index=str, columns={var_wnt: "data"})
    dfQ['Datetime'] = pd.to_datetime(dfQ.index)
    dfQ = dfQ.set_index('Datetime')

    dfQ['sd'] = dfQsd.values

    in1 = dfQ.index.month == 3
    in2 = dfQ.index.month == 6
    in3 = dfQ.index.month == 9
    in4 = dfQ.index.month == 12
    dfQ1 = dfQ[in1]
    dfQ2 = dfQ[in2]
    dfQ3 = dfQ[in3]
    dfQ4 = dfQ[in4]
    ts_lbl1 = 'Winter ' + title + ' ({} °N)'.format(int(lat_wnt[i]))
    ts_lbl2 = 'Spring ' + title + ' ({} °N)'.format(int(lat_wnt[i]))
    ts_lbl3 = 'Summer ' + title + ' ({} °N)'.format(int(lat_wnt[i]))
    ts_lbl4 = 'Fall ' + title + ' ({} °N)'.format(int(lat_wnt[i]))

    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfQ, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list = [dfQ1, dfQ2, dfQ3, dfQ4]

    metric_lbl = y_lbl.format(int(lat_wnt[i]))
    ts_lbl = [ts_lbl1, ts_lbl2, ts_lbl3, ts_lbl4]
    clmns_iea = ['year', 'time', 'index', 'error', 'SElo', 'SEup',
                 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

    df_all_list.append(df_list)
    ts_lbl_list.append(ts_lbl)
    metric_list.append([metric_lbl, metric_lbl, metric_lbl, metric_lbl])
    lon_list.append([lon1, lon1, lon1, lon1])
    lat_list.append([lat1, lat1, lat1, lat1])
    depth_list.append([depth, depth, depth, depth])

fn_out = '{}_S.csv'.format(file_pre)
df_flat = list(itertools.chain(*df_all_list))
ts_lbl_flat = list(itertools.chain(*ts_lbl_list))
metric_flat = list(itertools.chain(*metric_list))
lon_flat = list(itertools.chain(*lon_list))
lat_flat = list(itertools.chain(*lat_list))
depth_flat = list(itertools.chain(*depth_list))



fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_flat, lat_flat, lon_flat, depth_flat, metric_flat, ts_lbl_flat, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

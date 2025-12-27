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

# --MHW index wanted
mhw_wnt = ['heatwave_cover', 'sum_area', 'intensity']

# --input directory
dir_in = './data_gha/MHW/'

# --input file type extension
file_type_in = '.nc'

# --variable xr.ds name
var_wnt = mhw_wnt

# --attrributes xr.ds name
att_wnt = ['title']

# --IEA file names
file_pre = 'oc_mhw'

# --IEA file columns
y_lbl = mhw_wnt
lat = np.nan
lon = np.nan
depth = np.nan

# --names of time, depth, lat dimensions
var_time = 'time'

# --output CSV file
dir_out = './csv_for_erddap/'
yr_csv_bgn = 1983

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# size of input variables
num_mhw = len(mhw_wnt)

# list basin files that have been downloaded from erddap
files = os.listdir(dir_in)

indx = []
for i in range(num_mhw):
    for j in range(len(files)):
        if mhw_wnt[i] in files[j]:
            indx.append(j)

df_all_list = []
ts_lbl_list = []
metric_list = []
for i in range(0, num_mhw):
    # --create input filename
    file_mhw = files[indx[i]]

    fn_in = '{}{}'.format(dir_in, file_mhw)
    print(fn_in)

    # --open the netcdf files as an xarray
    dall = xr.open_dataset(fn_in)

    # --get the attributes
    title = dall.attrs[att_wnt[0]]

    # --xr.dataset as pd.df
    dfi = dall.to_dataframe()
    dfi = dfi.set_index('time')

    # the last year of the csv file
    yr_csv_end = dfi.index.year[-1]

    # --resample by taking 'M' means, i.e. monthly means
    dfQmn = dfi.resample('QE-MAR').mean()
    dfQsd = dfi.resample('QE-MAR').std()

    # --combine mean and sd into one pd.df
    dfQ = dfQmn.rename(index=str, columns={var_wnt[i]: "data"})
    # dfQ['Datetime'] = pd.to_datetime(dfQ.index, infer_datetime_format=True)
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
    ts_lbl1 = 'Winter ' + title
    ts_lbl2 = 'Spring ' + title
    ts_lbl3 = 'Summer ' + title
    ts_lbl4 = 'Fall ' + title

    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfQ, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list = [dfQ1, dfQ2, dfQ3, dfQ4]



    metric_lbl = y_lbl[i]
    ts_lbl = [ts_lbl1, ts_lbl2, ts_lbl3, ts_lbl4]
    clmns_iea = ['year', 'time', 'index', 'error', 'SElo', 'SEup',
                 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

    df_all_list.append(df_list)
    ts_lbl_list.append(ts_lbl)
    metric_list.append([metric_lbl, metric_lbl, metric_lbl, metric_lbl])


fn_out = '{}_S.csv'.format(file_pre)
df_flat = list(itertools.chain(*df_all_list))
ts_lbl_flat = list(itertools.chain(*ts_lbl_list))
metric_flat = list(itertools.chain(*metric_list))

fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_flat, lat, lon, depth, metric_flat, ts_lbl_flat, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

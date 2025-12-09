import itertools
import os
import xarray as xr
import pandas as pd
import numpy as np
from fun_pd_df2csvR_time import fun_pd_df2csvR_time
# -*- coding: utf-8 -*-
# pylint: disable=C0103


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# --basin index wanted
basin_wnt = ['ONI', 'PDO', 'NPGO']

# --input directory
dir_in = './data_gha/BasinScaleIndicators/'

# --input file type extension
file_type_in = '.nc'

# --variable xr.ds name
var_wnt = basin_wnt

# --attrributes xr.ds name
att_wnt = ['title']

# --IEA file names
file_pre = 'oc_'

# --IEA file columns
y_lbl = basin_wnt
lat = np.nan
lon = np.nan
depth = np.nan

# --names of time, depth, lat dimensions
var_time = 'time'

# --output CSV file
dir_out = './csv_for_erddap/'
yr_csv_bgn = 1900


# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# size of input variables
num_basin = len(basin_wnt)

# list basin files that have been downloaded from erddap
files = os.listdir(dir_in)

indx = []
for i in range(num_basin):
    for j in range(len(files)):
        if basin_wnt[i] in files[j]:
            indx.append(j)

df_all_list = []
ts_lbl_list = []
metric_list = []
for i in range(0, num_basin):
    # --create input filename
    file_basin = files[indx[i]]

    fn_in = '{}{}'.format(dir_in, file_basin)

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
    dfMmn = dfi.resample('ME').mean()
    dfMsd = dfi.resample('ME').std()

    # --combine mean and sd into one pd.df
    dfM = dfMmn.rename(index=str, columns={var_wnt[i]: "data"})
    # dfM['Datetime'] = pd.to_datetime(dfM.index, infer_datetime_format=True)
    dfM['Datetime'] = pd.to_datetime(dfM.index)
    dfM = dfM.set_index('Datetime')
    dfM['sd'] = dfMsd.values

    ts_lbl1 = 'Monthly ' + title
    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfM, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list = [dfM]
    num_order = len(df_list)

    # fn_out = file_pre + re.sub(r'[^\w]', '', basin_wnt[i]) + '.csv'
    fn_out = '{}{}_M.csv'.format(file_pre, basin_wnt[i])

    metric_lbl = y_lbl[i]
    ts_lbl = [ts_lbl1]
    clmns_iea = ['year', 'month', 'time', 'index', 'error', 'SElo', 'SEup',
                 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

    fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_list, lat, lon, depth,
                                     metric_lbl, ts_lbl, dir_out, fn_out,
                                     yr_csv_bgn, yr_csv_end)

metric_list.append([metric_lbl])
df_all_list.append(df_list)
ts_lbl_list.append(ts_lbl)

df_flat = list(itertools.chain(*df_all_list))
ts_lbl_flat = list(itertools.chain(*ts_lbl_list))
metric_flat = list(itertools.chain(*metric_list))

fn_out = '{}_M_all.csv'.format(file_pre)
fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_flat, lat, lon, depth, metric_flat, ts_lbl_flat, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

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
# lat_wnt = [33, 39, 45]
lat_wnt = [45, 39, 33]

# --input directory
dir_out = './data_gha/newUI/'
# dir_out = './data_x13/newUI/'
dir_in = dir_out
fn_in = '{}{}_monthly.nc'

# --input file type extension
file_type_in = '.nc'

# --variable xr.ds name
var_wnt = 'ui_mon'

# --dataset
dataset_wnt = ['CUTI', 'BEUTI']

# --attrributes xr.ds name
att_wnt = ['lat']

# --IEA file names
file_pre = ['oc_cuti', 'oc_beuti']

# --IEA file columns
y_lbl = ['CUTI (m^2 s^{-1})', 'BEUTI (mmol s^-1 m^-1)']
lat = lat_wnt

# --names of time, depth, lat dimensions
var_time = 'time'

# --output CSV file
dir_out = './csv_for_erddap/'
yr_csv_bgn = 1988

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# size of input variables
num_dataset_wnt = len(dataset_wnt)
num_lat = len(lat_wnt)

for iii in range(num_dataset_wnt):
    # open monthly netcdf
    ds1 = xr.open_dataset(fn_in.format(dir_in, dataset_wnt[iii]))

    df_all_list = []
    ts_lbl_list = []
    metric_list = []
    lon_list = []
    lat_list = []
    depth_list = []
    for i in range(0, num_lat):
        # open the netcdf files as an xarray
        dall = ds1[var_wnt].sel(lat=lat_wnt[i])

        # --get the attributes
        lat = dall[att_wnt[0]].data

        # set lon to 0
        lon = np.nan

        # set depth to 0
        depth = np.nan

        # --xr.dataset as pd.df
        dfi = dall.to_dataframe()

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
        title = dataset_wnt[iii]
        ts_lbl1 = 'Winter ' + title + ' ({} °N)'.format(int(lat_wnt[i]))
        ts_lbl2 = 'Spring ' + title + ' ({} °N)'.format(int(lat_wnt[i]))
        ts_lbl3 = 'Summer ' + title + ' ({} °N)'.format(int(lat_wnt[i]))
        ts_lbl4 = 'Fall ' + title + ' ({} °N)'.format(int(lat_wnt[i]))

        # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
        # --the case of monthly means this list only consists of dfQ, but a
        # --CSV for seasonal means will need a list of 4 pd.df for each season.
        df_list = [dfQ1, dfQ2, dfQ3, dfQ4]

        metric_lbl = y_lbl[iii]
        ts_lbl = [ts_lbl1, ts_lbl2, ts_lbl3, ts_lbl4]
        clmns_iea = ['year', 'time', 'index', 'error', 'SElo', 'SEup',
                     'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

        df_all_list.append(df_list)
        ts_lbl_list.append(ts_lbl)
        metric_list.append([metric_lbl, metric_lbl, metric_lbl, metric_lbl])
        lon_list.append([lon, lon, lon, lon])
        lat_list.append([lat, lat, lat, lat])
        depth_list.append([depth, depth, depth, depth])

    fn_out = '{}_S.csv'.format(file_pre[iii])

    df_flat = list(itertools.chain(*df_all_list))
    ts_lbl_flat = list(itertools.chain(*ts_lbl_list))
    metric_flat = list(itertools.chain(*metric_list))
    lon_flat = list(itertools.chain(*lon_list))
    lat_flat = list(itertools.chain(*lat_list))
    depth_flat = list(itertools.chain(*depth_list))



    fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_flat, lat_flat, lon_flat, depth_flat, metric_flat, ts_lbl_flat, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

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
depth = [np.nan, np.nan, np.nan]

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
    lat_list = []
    lon_list = []
    depth_list = []
    for i in range(0, num_lat):
        dall = ds1[var_wnt].sel(lat=lat_wnt[i])

        # --get the attributes
        lat = dall[att_wnt[0]].data
        lat_list.append(lat)

        # set lon to 0
        lon = np.nan
        lon_list.append(lon)

        # set depth to 0
        depth = np.nan
        depth_list.append(depth)

        # --get the attributes
        title = dataset_wnt[iii]

        # xr.dataset as pd.df
        dfi = dall.to_dataframe()

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

        metric_lbl = y_lbl[iii]
        ts_lbl = [ts_lbl1]

        df_all_list.append(df_list)
        ts_lbl_list.append(ts_lbl)
        metric_list.append([str(metric_lbl)])

    df_flat = list(itertools.chain(*df_all_list))
    ts_lbl_flat = list(itertools.chain(*ts_lbl_list))
    metric_flat = list(itertools.chain(*metric_list))

    clmns_iea = ['year', 'month', 'time', 'index', 'error', 'SElo', 'SEup',
                 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']
    fn_out = '{}_M.csv'.format(file_pre[iii])
    fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_flat, lat_list, lon_list, depth_list, metric_flat, ts_lbl_flat, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

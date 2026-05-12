import os
import calendar as clndr
import xarray as xr
import pandas as pd
import numpy as np


# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------

# .1) sst time series
fn_sst = './data_gha/SSTspatial/TS_monthly.nc'
fn1_hci = './data_gha/HCI/sst_oi_lat_30_48_xdis_75km.nc'
fn2_hci = './data_gha/HCI/sst_oi_lat_30_48_xdis_150km.nc'
var_list = ['sst', 'sst_oi', 'sst_oi']
lon_list = ['lon_vec', 'longitude', 'longitude']
lat_list = ['lat_vec', 'latitude', 'latitude']
lbl_list = ['NEP', '75km', '150km']

# clim time period
yr_bgn = 1982
yr_clim_bgn = 1991
yr_clim_end = 2020


roll_vec = np.arange(1, 31)

# .2) basin indices
dir_basin = './data_gha/BasinScaleIndicators/'

basin_wnt = ['PDO', 'ONI', 'NPGO']

# .3) output directory
dir_out = './data_gha/CalCofiSoccer/'

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# input variable size
num_roll = len(roll_vec)
num_basin = len(basin_wnt)


ds1 = xr.open_dataset(fn_sst)
ds2 = xr.open_dataset(fn1_hci)
ds3 = xr.open_dataset(fn2_hci)

yr1_end = ds1.time.dt.year.data[-1]
yr2_end = ds2.time.dt.year.data[-1]
yr3_end = ds3.time.dt.year.data[-1]

yr_end = np.max([yr1_end, yr2_end, yr3_end])
tt_wnt = np.arange(np.datetime64('{}-01'.format(yr_bgn), 'M'), np.datetime64('{}-01'.format(yr_end+1), 'M'))
nt = len(tt_wnt)

ds_list = [ds1, ds2, ds3]

for iii in range(len(ds_list)):
    da1 = ds_list[iii][var_list[iii]].sel(time=slice('{}-01-01'.format(yr_bgn), '{}-12-31'.format(yr_end)))

    da1_clim = da1.sel(time=slice('{}-01-01'.format(yr_clim_bgn), '{}-12-31'.format(yr_clim_end))).groupby('time.month').mean('time')
    da1_anom = da1.groupby('time.month') - da1_clim

    da1_mn = da1_anom.mean(dim=[lat_list[iii], lon_list[iii]])

    data1_mtrx = np.zeros([nt, num_roll])*np.nan
    for i in range(num_roll):
        da1_roll = da1_mn.rolling(time=roll_vec[i], center=True).mean().dropna("time")
        ia1 = np.isin(tt_wnt, da1_roll.time.data)
        ib1 = np.isin(da1_roll.time.data, tt_wnt)
        data1_mtrx[ia1, i] = da1_roll.data[ib1]

    da1_out = xr.DataArray(data1_mtrx, coords=[tt_wnt.astype('datetime64[ns]'), roll_vec], dims=['time', 'day_roll'])

    if iii == 0:
        ds1_out = da1_out.to_dataset(name=lbl_list[iii])
    else:
        ds1_out[lbl_list[iii]] = da1_out


# list basin files that have been downloaded from erddap
files = os.listdir(dir_basin)

indx = []
for i in range(num_basin):
    for j in range(len(files)):
        if basin_wnt[i] in files[j]:
            indx.append(j)

for i in range(num_basin):
    ds1b = xr.open_dataset('{}{}'.format(dir_basin, files[indx[i]]))
    da1b = ds1b[basin_wnt[i]]

    tt1b = ds1b['time'].data.astype('datetime64[M]')
    ia1 = np.isin(tt1b, tt_wnt)
    ib1 = np.isin(tt_wnt, tt1b)

    data_vec = np.zeros(nt)*np.nan
    data_vec[ib1] = da1b.data[ia1]
    
    da2_out = xr.DataArray(data_vec, coords=[tt_wnt.astype('datetime64[ns]')], dims=['time'])

    ds1_out[basin_wnt[i]] = da2_out

# make the output directory to save monthly and daily datasets
try:
    os.makedirs(dir_out)
except OSError:
    if not os.path.isdir(dir_out):
        raise

fn1_out = '{}basin_5rows.nc'.format(dir_out)
ds1_out.to_netcdf(fn1_out)

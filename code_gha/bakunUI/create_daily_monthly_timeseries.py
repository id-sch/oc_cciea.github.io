import os
import numpy as np
import xarray as xr
import pandas as pd
# import matplotlib as mpl


# Note: M-x pyvenv-workon py_cart
#       This creates daily, monthly means and cui_mtrx from the 6hr UI data
#       Save it to netcdf file


# -------------------------------------------------------
# -- Input variables, change these
# -------------------------------------------------------
# set lat range
lat_bgn = 24
lat_end = 54
dlat = 3

# title for the output dataset attributes
title_lbl = 'Upwelling Index'

# UI file type
file_type = 'nc'

# UI files have starting month and years
mon_bgn = 1
yr_bgn = 1967

# variable name in the xr.ds
var_lbl = ['upwelling_index']

# directory to save the netcdf files
file_type = 'nc'
dir_out = './data_gha/bakunUI/'
dir_in = dir_out

# -------------------------------------------------------
# -- END: Input variables, change these
# -------------------------------------------------------

# lat range
lat_wnt = np.arange(lat_bgn, lat_end, dlat)
num_lat = len(lat_wnt)

# get the final set of files that have the 6 hr data
files = os.listdir(dir_in)

indx = []
for i in range(num_lat):
    for j in range(len(files)):
        if str(lat_wnt[i]) in files[j]:
            indx.append(j)

files_wnt = np.array(files)[indx]

# get the end year and month
yr_end = int(files_wnt[0][6:10])
mon_end = int(files_wnt[0][11:13])

# create daily and monthly values for the input time range
date_bgn = '{}-{:02d}-01'.format(yr_bgn, mon_bgn)
if mon_end == 12:
    date_end = '{}-{:02d}-01'.format(yr_end+1, 1)
else:
    date_end = '{}-{:02d}-01'.format(yr_end, mon_end+1)

dateD = np.arange(date_bgn, date_end, dtype='datetime64[D]')
dateM = np.arange(date_bgn, date_end, dtype='datetime64[M]')

ntD = len(dateD)
ntM = len(dateM)

# place dateD into pandas datetime to extract the year
dateD_yy = pd.to_datetime(dateD).year
yrs = np.unique(dateD_yy)
num_yrs = len(yrs)

# loop over all lats
dataD_mtrx = np.zeros([num_lat, ntD])*np.nan
dataM_mtrx = np.zeros([num_lat, ntM])*np.nan
cui_mtrx = np.zeros([num_lat, num_yrs, 365])
lon_vec = np.zeros([num_lat])
for i in range(0, num_lat):
    dir_fn = '{}{}'.format(dir_in, files_wnt[i])
    # --open the netcdf files as an xarray
    ds = xr.open_dataset(dir_fn)

    # longitude
    lon_vec[i] = int(ds['longitude'].data[0])

    # --since time is not a dimension (it is a variable) in ERDDAP's netcdf,
    # --we need to create a new xarray that has time as a dimension
    ds1 = xr.Dataset({'upwelling_index': (('time'), ds.upwelling_index.data)},
                     {'time': ds.time.data, 'location': lat_wnt[i]})

    # --create daily means
    ds1D = ds1.resample(time='D').mean('time')

    # get index for output matrixD
    iaD = np.isin(dateD, ds1D.time)
    ibD = np.isin(ds1D.time, dateD)
    dataD_mtrx[i, iaD] = ds1D[var_lbl[0]].data[ibD]

    # --create monthly means
    ds1M = ds1.resample(time='ME').mean('time')

    # get index for output matrixM
    iaM = np.isin(dateM, ds1M.time.data.astype('datetime64[M]'))
    ibM = np.isin(ds1M.time.data.astype('datetime64[M]'), dateM)
    dataM_mtrx[i, iaM] = ds1M[var_lbl[0]].data[ibM]

    # create CUI matrix
    for j in range(0, num_yrs):
        in_yr = np.where(dateD_yy == yrs[j])[0]
        ui_yr = dataD_mtrx[i, in_yr].T
        # --check size of in_yr, can be 365, 366 or
        # --less (depending on mon_wnt1,mon_wnt2)
        num_in = np.size(in_yr)

        in_end = 365
        if num_in < in_end:
            in_end = num_in-1

        # calculate cui on the first 365 days
        ui_365 = np.zeros(365)*np.nan
        ui_365[0:in_end] = ui_yr[0:in_end]
        cui = np.nancumsum(ui_365)

        # nancumsum treats NaN as 0, but want NaN in output
        in_nan = np.isnan(ui_365)
        cui[in_nan] = np.nan

        # place in final matrix
        cui_mtrx[i, j, :] = cui

# put into xr.da
da1 = xr.DataArray(dataD_mtrx, coords=[lat_wnt, dateD.astype(
    'datetime64[ns]')], dims=['lat', 'time'])
da2 = xr.DataArray(dataM_mtrx, coords=[lat_wnt, dateM.astype(
    'datetime64[ns]')], dims=['lat', 'time'])
days = np.arange(1, 366)
da3 = xr.DataArray(cui_mtrx, coords=[lat_wnt, yrs, days],
                   dims=['lat', 'year', 'days'])
da4 = xr.DataArray(lon_vec, coords=[lat_wnt], dims=['lat'])

# put into xr.ds
ds1 = da1.to_dataset(name='ui_day')
ds2 = da2.to_dataset(name='ui_mon')
ds3 = da3.to_dataset(name='cui_mtrx')

ds1['lon'] = da4
ds2['lon'] = da4
ds3['lon'] = da4

# add the title
ds1.attrs['title'] = title_lbl
ds2.attrs['title'] = title_lbl
ds3.attrs['title'] = title_lbl


# --check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_out)
except OSError:
    if not os.path.isdir(dir_out):
        raise

# --Save Dataset to a netcdf file
print(fn1_nc)
print(fn2_nc)
print(fn3_nc)

fn1_nc = '{}/UI_daily.nc'.format(dir_out)
ds1.to_netcdf(fn1_nc)

fn2_nc = '{}/UI_monthly.nc'.format(dir_out)
ds2.to_netcdf(fn2_nc)

fn3_nc = '{}/UI_cui_mtrx.nc'.format(dir_out)
ds3.to_netcdf(fn3_nc)

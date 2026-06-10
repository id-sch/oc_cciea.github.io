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
lat_bgn = 31
lat_end = 47
dlat = 1

# variable name in the xr.d
ds1_var = ['year', 'month', 'day']
ds1_coord = ['latitude', 'time']

# input netcdf file are daily data, set the minumum number of days
# in the last month to create a monthly mean, otherwise report NaN
ndm_cutoff = 15

# directory to save the netcdf files
file_type = 'nc'
dir_out = './data_gha/newUI/'
# dir_out ='./data_x13/newUI/'
dir_in = dir_out

# -------------------------------------------------------
# -- END: Input variables, change these
# -------------------------------------------------------

# lat range
lat_wnt = np.arange(lat_bgn, lat_end, dlat)
num_lat = len(lat_wnt)

# get the final set of files that have the 6 hr data
files = os.listdir(dir_in)
num_files = len(files)

#
for iii in range(num_files):
    var_wnt = files[iii].split('_')[0]
    fn1 = '{}{}'.format(dir_in, files[iii])
    ds1 = xr.open_dataset(fn1)

    # create proper time coord
    year = ds1[ds1_var[0]].data
    month = ds1[ds1_var[1]].data
    day = ds1[ds1_var[2]].data

    yrs = np.unique(year)
    num_yrs = len(yrs)

    # create dictionary of dates
    time_dic = {}
    time_dic[ds1_var[0]] = year
    time_dic[ds1_var[1]] = month
    time_dic[ds1_var[2]] = day

    # time pd.dataframe
    df_vec = pd.DataFrame(time_dic)
    pd_dt = pd.DatetimeIndex(pd.to_datetime(df_vec[ds1_var[0:3]]))
    dt_vec = pd.DatetimeIndex(pd.to_datetime(df_vec[ds1_var[0:3]])).values

    # lat range
    lat_rng = np.arange(lat_bgn, lat_end+1, dlat)
    num_lat_rng = len(lat_rng)


    # get index of lat range
    cc, ia_lat, ib_lat = np.intersect1d(
        ds1[ds1_coord[0]].data, lat_rng, return_indices=True)

    # daily matrix
    dataD_mtrx = ds1[var_wnt].data[:, ia_lat].T

    # monthly matrix
    da1 = xr.DataArray(ds1[var_wnt].data[:, ia_lat], coords=[
                       dt_vec, lat_rng], dims=['time', 'latitude'])
    da1M = da1.resample(time='ME').mean('time')
    dataM_mtrx = da1M.data.T
    dateM = da1M.time.data.astype('datetime64[M]')

    # check to see if number of days in last month is less than ndm_cutoff
    ndm_last_mon = pd_dt.day.values[-1]
    if ndm_last_mon < ndm_cutoff:
        dataM_mtrx[:, -1] = np.nan

    # create CUI matrix
    dateD_yy = da1.time.dt.year.data
    cui_mtrx = np.zeros([num_lat_rng, num_yrs, 365])
    for i in range(0, num_lat_rng):
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
    da1 = xr.DataArray(dataD_mtrx, coords=[lat_rng, dt_vec], dims=['lat', 'time'])
    da2 = xr.DataArray(dataM_mtrx, coords=[lat_rng, dateM.astype(
        'datetime64[ns]')], dims=['lat', 'time'])
    days = np.arange(1, 366)
    da3 = xr.DataArray(cui_mtrx, coords=[lat_rng, yrs, days],
                       dims=['lat', 'year', 'days'])

    # put into xr.ds
    ds1_out = da1.to_dataset(name='ui_day')
    ds2_out = da2.to_dataset(name='ui_mon')
    ds3_out = da3.to_dataset(name='cui_mtrx')

    # add the title
    ds1_out.attrs['title'] = var_wnt
    ds2_out.attrs['title'] = var_wnt
    ds3_out.attrs['title'] = var_wnt

    # --check if directory exist, if it doesn't then create
    try:
        os.makedirs(dir_out)
    except OSError:
        if not os.path.isdir(dir_out):
            raise

    # --Save Dataset to a netcdf file
    fn1_nc = '{}/{}_daily.nc'.format(dir_out, var_wnt)
    ds1_out.to_netcdf(fn1_nc)

    fn2_nc = '{}/{}_monthly.nc'.format(dir_out, var_wnt)
    ds2_out.to_netcdf(fn2_nc)

    fn3_nc = '{}/{}_cui_mtrx.nc'.format(dir_out, var_wnt)
    ds3_out.to_netcdf(fn3_nc)

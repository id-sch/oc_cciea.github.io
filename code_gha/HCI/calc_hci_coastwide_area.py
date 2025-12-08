import calendar as clndr
import xarray as xr
import pandas as pd
import numpy as np


# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------
# x distance and lat range list
xdis = 150

test_rgn = 0

if test_rgn == 0:
    lat_rgn = [[43.5, 48], [40, 43.5], [35.5, 40], [30, 35.5]]
else:
    lat_rgn = [[43, 48], [39, 43], [34.5, 39], [30, 34.5]]

# variable wanted
var_wnt = 'sst_oi'

# x distance threshold
xdis_thrsh = 75

# clim threshold
clim_bgn = 1982
clim_end = 2010


# Months that define each season
mon_ssn = [np.array([12, 1, 2]), np.array([3, 4, 5]),
           np.array([6, 7, 8]), np.array([9, 10, 11])]

# ouput directory
dir_out = './data_gha/HCI/'

# directory of temperature data
dir_in = dir_out

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# input variable size
num_lat_rgn = len(lat_rgn)

for i in range(num_lat_rgn):
    # 1) open threshold data
    lat_rgn1 = lat_rgn[i]
    fn1 = '{}/{}_lat_{}_{}_xdis_{}km.nc'.format(
        dir_in, var_wnt, lat_rgn1[0], lat_rgn1[1], xdis)
    ds1 = xr.open_dataset(fn1)

    # var and coords
    var_ds1 = list(ds1.keys())
    coord_ds1 = list(ds1.coords.keys())

    # get temp_mtrx
    da1 = ds1[var_wnt]

    # get area matrix
    da1_area = ds1['area']

    # temp_mtrx is over a large domain find x limits
    da1_mn = da1.mean('time').mean('latitude')
    ind = np.isfinite(da1_mn)
    da1_lon_rng = da1_mn[ind]
    lon_bgn = da1_lon_rng.longitude[0].data
    lon_end = da1_lon_rng.longitude[-1].data

    da1f = da1.sel(longitude=slice(lon_bgn, lon_end))
    da1f_area = da1_area.sel(longitude=slice(lon_bgn, lon_end))

    # size of coords
    nt, ny, nx = da1f.shape

    # lat and lon
    lat1f = da1f.latitude.data
    lon1f = da1f.longitude.data

    # yr, mon vectors
    time1 = da1f.time.data
    yr_vec = da1f.time.dt.year.data
    mon_vec = da1f.time.dt.month.data

    # 2) open the threshold data
    fn2 = '{}/threshold_{}_lat_{}_{}_xdis_{}km_clim_{}_{}.nc'.format(
        dir_in, var_wnt, lat_rgn1[0], lat_rgn1[1], xdis_thrsh, clim_bgn, clim_end)

    print(fn2)
    ds2 = xr.open_dataset(fn2)

    # there are a monthly and seasonal threshold to calc the HCI for
    var_ds2 = list(ds2.keys())
    coord_ds2 = list(ds2.coords.keys())

    # save of ds2 variable and coords
    num_thrsh = len(var_ds2)
    num_ssn = len(ds2[coord_ds2[0]])

    # loop over all dates
    thrshld_vec = np.zeros([nt, num_thrsh])
    lvlM_old_mtrx = np.zeros([nt, ny, nx, num_thrsh])
    lvlM_new_mtrx = np.zeros([nt, ny, nx, num_thrsh])
    ts_num_below_old = np.zeros([nt, num_thrsh])
    ts_num_below_new = np.zeros([nt, num_thrsh])
    for j in range(nt):
        dM = da1f.data[j, :, :]

        # dM has nan for land and data over xdis from shore, get the total
        # number of possible data points in the domain
        ind = np.isfinite(dM).nonzero()
        num_normalize = len(ind[0])
        area_normalize = np.nansum(da1f_area.data)

        # for this month find the correct season, ie see if it is in winter, spring, summer or fall
        in_ssn = np.where(mon_ssn - mon_vec[j] == 0)[0]

        # loop over all thresholds, need to check to see if seasonal or monthly
        for k in range(num_thrsh):
            # check to see if seasonal and get threshold
            if coord_ds2[k] == 'season':
                in_ssn = np.where(mon_ssn - mon_vec[j] == 0)[0]
                thrsh_wnt = ds2[var_ds2[k]].data[in_ssn]

            # check to see if monthlly and get threshold
            if coord_ds2[k] == 'month':
                mon_lbl = clndr.month_name[mon_vec[j]][0:3]
                in_mon = np.where(ds2['month'].data == mon_lbl)[0]
                thrsh_wnt = ds2[var_ds2[k]].data[in_mon]

            # place in a vector, so can have time series of thresholds
            thrshld_vec[j, k] = thrsh_wnt[0]

            # 1) old method, summing up grid cells
            # create a matrix, 1=below threshold, 0=below threshold
            ones_mtrx = np.zeros([ny, nx])
            chck_lvl = dM - thrsh_wnt
            in_neg = np.where(chck_lvl < 0)
            ones_mtrx[in_neg] = 1
            lvlM_old_mtrx[j, :, :, k] = ones_mtrx

            # find where ones_mtrx is 1, this will be the HCI
            in_num = np.where(ones_mtrx == 1)
            # ts_num_below[j, k] = len(in_num[0])
            ts_num_below_old[j, k] = len(in_num[0])/num_normalize

            # 2) new method, summing up areas of grid cells
            ones_area_mtrx = np.zeros([ny, nx])
            ones_area_mtrx[in_neg] = da1f_area.data[in_neg]
            lvlM_new_mtrx[j, :, :, k] = ones_area_mtrx

            # sum grid cells that have area and normalize by total area
            ts_num_below_new[j, k] = np.sum(ones_area_mtrx)/area_normalize

    # dataarray
    da1_out = xr.DataArray(
        thrshld_vec,
        coords=[time1.astype('datetime64'), var_ds2],
        dims=['time', 'threshold'])

    da2_old_out = xr.DataArray(
        ts_num_below_old,
        coords=[time1.astype('datetime64'), var_ds2],
        dims=['time', 'threshold'])

    da4_old_out = xr.DataArray(
        lvlM_old_mtrx,
        coords=[time1.astype('datetime64'), lat1f, lon1f, var_ds2],
        dims=['time', 'latitude', 'longitude', 'threshold'])

    da2_new_out = xr.DataArray(
        ts_num_below_new,
        coords=[time1.astype('datetime64'), var_ds2],
        dims=['time', 'threshold'])

    da4_new_out = xr.DataArray(
        lvlM_new_mtrx,
        coords=[time1.astype('datetime64'), lat1f, lon1f, var_ds2],
        dims=['time', 'latitude', 'longitude', 'threshold'])

    # dataset, save both old and new method
    ds1_out = da1f.to_dataset(name=var_wnt)
    ds1_out['thresholdM_vec'] = da1_out
    ds1_out['hci_old'] = da2_old_out
    ds1_out['lvlM_old_mtrx'] = da4_old_out
    ds1_out['hci'] = da2_new_out
    ds1_out['lvlM_mtrx'] = da4_new_out

    # place the thresholds into the dataset
    for k in range(num_thrsh):
        ds1_out[var_ds2[k]] = ds2[var_ds2[k]]

    # save to netcdf
    fn_out = '{}hci_{}_lat_{}_{}_xdis_{}km_thresh_xdis_{}km.nc'.format(
        dir_out, var_wnt, lat_rgn1[0], lat_rgn1[1], xdis, xdis_thrsh)
    ds1_out.to_netcdf(fn_out)

    # 2) CSV files
    # NUM BELOW
    for iii in range(num_thrsh):
        pd1 = pd.DataFrame({'fraction below threshold "{}"'.format(var_ds2[iii]): ts_num_below_new[:, iii]})
        pd1['time'] = pd.to_datetime(time1)
        pd1 = pd1.set_index('time')
        fn_csv = '{}fraction_below_Threshold_{}_Lat{}_{}_X{}_Clim{}_{}.csv'.format(
            dir_out, var_ds2[iii], lat_rgn1[0], lat_rgn1[1], xdis, clim_bgn, clim_end)
        pd1.to_csv(fn_csv)

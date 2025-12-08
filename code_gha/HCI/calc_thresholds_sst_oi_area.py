import xarray as xr
import numpy as np


# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------
# second x_pnts is used to define the clim threshold
lat_rgn = [[43.5, 48], [40, 43.5], [35.5, 40], [30, 35.5]]


# define the clim based on the 75 km region
xdis = 75

# variable at z_wnt
var_wnt = 'sst_oi'

# clim time period
date_clim = [np.datetime64('1982-01-01'), np.datetime64('2011-01-01')]

# season order wanted
ssn_order = ['DJF', 'MAM', 'JJA', 'SON']
tt_ssn = [1, 4, 7, 10]
indx_ssn = [np.array([11, 0, 1]), np.array([2, 3, 4]),
            np.array([5, 6, 7]), np.array([8, 9, 10])]

# months
mon_order = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
tt_mon = np.arange(1, 13)
indx_ssn = np.arange(0, 12)

# ouput directory
dir_out = './data_gha/HCI/'

# directory of temperature data
dir_in = dir_out

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# len input variables
num_lat_rgn = len(lat_rgn)
num_ssn = len(ssn_order)
num_mon = len(mon_order)

# Months covering the clim time period
dates_wnt = np.arange(date_clim[0], date_clim[1], dtype='datetime64[M]')

# Threshold scenario: 
# a) seasonal averages of the regional mean (denoted by _mn_lon)
# Here calculate these values, will have final dataarray with size (time X lat)

for i in range(num_lat_rgn):
    lat_rgn1 = lat_rgn[i]
    fn1 = '{}/{}_lat_{}_{}_xdis_{}km.nc'.format(
        dir_in, var_wnt, lat_rgn1[0], lat_rgn1[1], xdis)

    ds1 = xr.open_dataset(fn1)

    # dataarray with the data we want
    da1 = ds1[var_wnt]

    # get da1 over clim time period
    da2 = da1.sel(time = dates_wnt)

    # spatial mean (in the future might want to explore weighted means, but don't think it will make much of a difference due to small latitude extend)
    # da2_mn = da2.mean(('lon', 'lat'))
    da2_mn = da2.mean(('longitude', 'latitude'))

    # 1) Seasonal
    # seasonal mean of the spatial clim time series
    da2_mn_ssn = da2_mn.groupby('time.season').mean('time')

    # get the correct order and round to the nearest half-degree
    threshold1_vec = np.zeros(num_ssn)
    for iii in range(num_ssn):
        in_ssn = np.where(da2_mn_ssn.season.data == ssn_order[iii])[0]
        threshold1_vec[iii] = np.round(2*da2_mn_ssn.data[in_ssn])/2

    # write the seasonal thresholds to dataarrays
    da1S = xr.DataArray(threshold1_vec, coords=[ssn_order], dims=['season'])

    # write to dataset
    ds1_out = da1S.to_dataset(name='mn_{}kmS'.format(xdis))

    # 2) Monthly
    # monthly mean of the spatial clim time series
    da2_mn_mon = da2_mn.groupby('time.month').mean('time')

    # get the monthly data to the nearest half-degree
    threshold1_mon_vec = np.zeros(num_mon)
    for iii in range(num_mon):
        in_mon = np.where(da2_mn_mon.month.data == tt_mon[iii])[0]
        threshold1_mon_vec[iii] = np.round(2*da2_mn_mon[in_mon])/2

    # write the monthly thresholds ot dataarrays
    da1M = xr.DataArray(threshold1_mon_vec, coords=[mon_order], dims=['month'])

    # write to dataset
    ds1_out['mn_{}kmM'.format(xdis)] = da1M

    # 3) Save the netcdf file
    # clim years
    clim_bgn = da2_mn.time.dt.year[0].data
    clim_end = da2_mn.time.dt.year[-1].data

    # print
    print(ds1_out.mn_75kmS.data)
    print(ds1_out.mn_75kmM.data)
    print(da2_mn_mon.data)

    # save to netcdf
    fn_out = '{}threshold_{}_lat_{}_{}_xdis_{}km_clim_{}_{}.nc'.format(
        dir_out, var_wnt, lat_rgn1[0], lat_rgn1[1], xdis, clim_bgn, clim_end)
    ds1_out.to_netcdf(fn_out)

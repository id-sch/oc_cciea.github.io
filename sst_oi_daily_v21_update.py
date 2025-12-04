import xarray as xr
import numpy as np
import os
from datetime import date
import urllib.request  # py36
import requests  # Aug 9, 2024 -- changed to request as urlib gave 404 error
from shutil import copyfile
import calendar as clndr



# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------
name_wnt = 'TS_monthly'

today = date.today()


mon_end = today.month
yr_end = today.year


file_type = 'nc'

dir_in ='./'


url_nrt_time = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21NrtAgg.nc?time'
fn1_time = 'time_nrt.nc'

# lat, lon box
lat1 = 20
lat2 = 61
lon1 = 180
lon2 = 250

# erddap url
url = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/'
var_name_erddap = 'ncdcOisst21NrtAgg'
file_type_in = 'nc'
var_wnt = 'sst'

# only download if there are these number of days available in the month
day_check = 24

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# len input variables


dir_list = os.listdir()
print("START -------------------------------")
print("Files and directories in  :")
print(dir_list)


# Get the time available in the NRT SST OI data product
ddd = requests.get(url_nrt_time, timeout=20)

with open(fn1_time, 'wb') as file_time:
    file_time.write(ddd.content)

ds1nrt = xr.open_dataset(fn1_time)
time_nrt = ds1nrt.time.data.astype('datetime64[D]')
yy_nrt = ds1nrt.time.dt.year.data
mm_nrt = ds1nrt.time.dt.month.data

yr_end_nrt = ds1nrt.time.dt.year.data[-1]
mon_end_nrt = ds1nrt.time.dt.month.data[-1]
day_end_nrt = ds1nrt.time.dt.day.data[-1]

date2 = np.datetime64('{}-{:02d}'.format(yr_end_nrt, mon_end_nrt), 'M')

# Open the old monthly data SST OI data
file_in = dir_in + '{}.{}'.format(name_wnt, file_type)


# set string for altitude, lon, lat ranges
alt_str = '[(0.0):1:(0.0)]'
lon_str = '[({:.1f}):1:({:.1f})]'.format(lon1, lon2)
lat_str = '[({:.1f}):1:({:.1f})]'.format(lat1, lat2)


# Check if the file exists and is a file
if os.path.isfile(file_in):
    print("File exists")

    ds1 = xr.open_dataset(file_in)

    yr_end = ds1.time.dt.year.data[-1]
    mon_end = ds1.time.dt.month.data[-1]
    date1 = np.datetime64('{}-{:02d}'.format(yr_end, mon_end), 'M')

    # check to see if number of days downloaded for the last month is not complete
    # redownload the last month if it isn't, otherwise skip to the next month
    num_day_miss_end = ds1.day_missing.data[-1]
    if num_day_miss_end > 0:
        dates_wnt = np.arange(date1, date2+1)
    else:
        dates_wnt = np.arange(date1+1, date2+1)

    # loop over the dates and save the daily data for each month
    ntM = len(dates_wnt)
    date_final = dates_wnt[0]
    for i in range(ntM):
        yri = dates_wnt[i].astype(object).year
        moni = dates_wnt[i].astype(object).month
        in_yr = np.where(yy_nrt == yri)[0]
        in_mon = np.where(mm_nrt[in_yr] == moni)[0]
        time_nrt1 = time_nrt[in_yr][in_mon]

        yr1 = time_nrt1[0].astype(object).year
        mon1 = time_nrt1[0].astype(object).month
        day1 = time_nrt1[0].astype(object).day
        yr2 = time_nrt1[-1].astype(object).year
        mon2 = time_nrt1[-1].astype(object).month
        day2 = time_nrt1[-1].astype(object).day

        num_days_available = day2 - day1
        if num_days_available > day_check:
            date_final = dates_wnt[i]

            time_str = '[({}-{:02d}-{:02d}T00:00:00Z):1:({}-{:02d}-{:02d}T00:00:00z)]'.format(yr1, mon1, day1, yr2, mon2, day2)

            # ERDDAP restful html
            urlf = '{}{}.{}?{}{}{}{}{}'.format(
                url, var_name_erddap, file_type_in, var_wnt, time_str, alt_str,
                lat_str, lon_str)

            # --create output directory
            dir_out = '{}/{:02d}/'.format(yri, moni)

            # --check if directory exist, if it doesn't then create
            try:
                os.makedirs(dir_out)
            except OSError:
                if not os.path.isdir(dir_out):
                    raise

            # download the data
            data = urllib.request.urlretrieve(urlf)  # Python

            # create a xarray dataset
            ds1e = xr.open_dataset(data[0])
            in_mon = np.where(ds1e.time.dt.month == moni)[0]
            data1 = np.squeeze(ds1e['sst'].data[in_mon, :, :, :])
            time1 = ds1e.time.data[in_mon]
            lat1 = ds1e.latitude.data
            lon1 = ds1e.longitude.data
            daf = xr.DataArray(data1, coords=[time1.astype('datetime64'), lat1, lon1],
                               dims=['time', 'latitude', 'longitude'])

            dsf = daf.to_dataset(name='sst')

            # save to final netcdf file
            fn_out_final = '{}{}_daily_final.{}'.format(dir_out, var_wnt, file_type)
            dsf.to_netcdf(fn_out_final)

    # create monthly means and append to existing dataset
    nt, ny, nx = ds1[var_wnt].shape
    time1M = ds1.time.data.astype('datetime64[M]')
    timefM = np.arange(time1M[0], date_final + 1)
    ntf = len(timefM)
    data_final = np.zeros([ntf, ny, nx])*np.nan
    day_missing_final = np.zeros(ntf)*np.nan

    ia = np.isin(timefM, time1M)
    ib = np.isin(time1M, timefM)

    data_final[ia, :, :] = ds1[var_wnt].data[ib, :, :]
    day_missing_final[ia] = ds1['day_missing'].data

    dates_calcM = np.setdiff1d(timefM,time1M)
    for i in range(len(dates_calcM)):
        yr1 = dates_calcM[i].astype(object).year
        mon1 = dates_calcM[i].astype(object).month
        dir1 = '{}/{:02d}/'.format(yr1, mon1)
        fn1 = '{}{}_daily_final.{}'.format(dir1, var_wnt, file_type)
        ds1D = xr.open_dataset(fn1)

        in1 = np.where(timefM == dates_calcM[i])[0]
        data_final[in1, :, :] = ds1D[var_wnt].mean('time')

        num_day_data = len(ds1D.time.data)
        num_day_clndr = clndr.monthrange(yr1, mon1)[1]
        num_day_diff = num_day_clndr - num_day_data

        day_missing_final[in1] = num_day_diff

    da1_out = xr.DataArray(
        data_final, coords=[timefM.astype('datetime64[ns]'), lat1, lon1],
        dims=['time', 'lat_vec', 'lon_vec'])
    da2_out = xr.DataArray(
        day_missing_final, [timefM.astype('datetime64[ns]')], dims=['time']) 

    ds1_out = da1_out.to_dataset(name=var_wnt)
    ds1_out['day_missing'] = da2_out

    # overwrite the existing file_in with this new updated dataset
    ds1.close()
    file_out = 'TS_monthly_update.nc'
    ds1_out.to_netcdf(file_out)
else:
    print('file not found: {}'.format(file_in))

dir_list_end = os.listdir()
print("END -------------------------------")
print("Files and directories in  :")
print(dir_list_end)

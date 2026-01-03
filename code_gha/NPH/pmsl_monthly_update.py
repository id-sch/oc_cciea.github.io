import calendar as clndr
import os
import shutil
import urllib.request
import requests  
import xarray as xr
import numpy as np


# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------
# ouput dataset name, file type, and directory
name_wnt = 'pmsl_monthly_update'
file_type = 'nc'
dir_in ='./data_gha/NPH/'

# erddap url for all dates of the SST OI
url_nph_time = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdlasFnWPr.nc?time'
fn1_time = 'time_nph.nc'

# lat, lon box
lat1 = 0.5
lat2 = 65.5
lon1 = 150.5
lon2 = 250.5

# erddap variables
url = 'https://upwell.pfeg.noaa.gov/erddap/griddap/'
var_wnt = 'pmsl'
var_name_erddap = 'erdlasFnWPr'
file_type_in = 'nc'

# dir out, will use artifacts to download this data
dir_out = './data_gha/NPH/'
dir_download_out = './data_x13/NPH/'

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# len input variables


dir_list = os.listdir()
print("START -------------------------------")
print("Files and directories in  :")
print(dir_list)


# Get the time available in the NRT SST OI data product
ddd = requests.get(url_nph_time, timeout=800)

with open(fn1_time, 'wb') as file_time:
    file_time.write(ddd.content)

ds1nph = xr.open_dataset(fn1_time)
time_nph = ds1nph.time.data.astype('datetime64[D]')
yy_nph = ds1nph.time.dt.year.data
mm_nph = ds1nph.time.dt.month.data

yr_end_nph = ds1nph.time.dt.year.data[-1]
mon_end_nph = ds1nph.time.dt.month.data[-1]
day_end_nph = ds1nph.time.dt.day.data[-1]

date2 = np.datetime64('{}-{:02d}'.format(yr_end_nph, mon_end_nph), 'M')

# Open the old monthly data SST OI data
file_in = dir_in + '{}.{}'.format(name_wnt, file_type)

# set string for lon, lat ranges
lon_str = '[({:.1f}):1:({:.1f})]'.format(lon1, lon2)
lat_str = '[({:.1f}):1:({:.1f})]'.format(lat1, lat2)


# Check if the file exists and is a file
if os.path.isfile(file_in):
    print("File exists")

    ds1 = xr.open_dataset(file_in)
    lat_vec = ds1['lat_vec'].data
    lon_vec = ds1['lon_vec'].data

    yr_end = ds1.time.dt.year.data[-1]
    mon_end = ds1.time.dt.month.data[-1]
    date1 = np.datetime64('{}-{:02d}'.format(yr_end, mon_end), 'M')

    # make the output directory to save monthly and daily datasets
    try:
        os.makedirs(dir_out)
    except OSError:
        if not os.path.isdir(dir_out):
            raise

    # download these months
    dates_wnt = np.arange(date1+1, date2+1)

    # loop over the dates and save the daily data for each month
    ntM = len(dates_wnt)

    if ntM > 0:
        date_final = dates_wnt[0]
        for i in range(ntM):
            yri = dates_wnt[i].astype(object).year
            moni = dates_wnt[i].astype(object).month
            in_yr = np.where(yy_nph == yri)[0]
            in_mon = np.where(mm_nph[in_yr] == moni)[0]
            time_nrt1 = time_nph[in_yr][in_mon]
            
            yr1 = time_nrt1[0].astype(object).year
            mon1 = time_nrt1[0].astype(object).month
            day1 = time_nrt1[0].astype(object).day
            yr2 = time_nrt1[-1].astype(object).year
            mon2 = time_nrt1[-1].astype(object).month
            day2 = time_nrt1[-1].astype(object).day

            date_final = dates_wnt[i]

            time_str = '[({}-{:02d}-{:02d}T00:00:00Z):1:({}-{:02d}-{:02d}T00:00:00z)]'.format(yr1, mon1, day1, yr2, mon2, day2)

            # print the year and month that will be downloaded
            print('pmsl download: {} {}, {}: {}'.format(yr1, mon1, i, ntM))

            # ERDDAP restful html
            urlf = '{}{}.{}?{}{}{}{}'.format(
                url, var_name_erddap, file_type_in, var_wnt, time_str,
                lat_str, lon_str)

            # --create output directory for daily data
            dir1 = '{}/{}/{:02d}/'.format(dir_download_out, yri, moni)

            # --check if directory exist, if it doesn't then create
            try:
                os.makedirs(dir1)
            except OSError:
                if not os.path.isdir(dir1):
                    raise

            # download the data
            fn1_pmsl='pmsl.nc'
            ddd = requests.get(urlf, timeout=800)
            with open(fn1_pmsl, 'wb') as file_pmsl:
                file_pmsl.write(ddd.content)

            # create a xarray dataset
            dsf = xr.open_dataset(fn1_pmsl)

            # save to final netcdf file
            fn_out_final = '{}{}_monthly_new.{}'.format(dir1, var_wnt, file_type)
            dsf.to_netcdf(fn_out_final)

        # create monthly means and append to existing dataset
        nt, ny, nx = ds1[var_wnt].shape
        time1M = ds1.time.data.astype('datetime64[M]')
        timefM = np.arange(time1M[0], date_final + 1)
        ntf = len(timefM)
        data_final = np.zeros([ntf, ny, nx])*np.nan

        ia = np.isin(timefM, time1M)
        ib = np.isin(time1M, timefM)

        data_final[ia, :, :] = ds1[var_wnt].data[ib, :, :]

        dates_calcM = np.setdiff1d(timefM,time1M)
        for i in range(len(dates_calcM)):
            yr1 = dates_calcM[i].astype(object).year
            mon1 = dates_calcM[i].astype(object).month
            dir1 = '{}/{}/{:02d}/'.format(dir_download_out, yr1, mon1)
            fn1 = '{}{}_monthly_new.{}'.format(dir1, var_wnt, file_type)
            ds1D = xr.open_dataset(fn1)

            in1 = np.where(timefM == dates_calcM[i])[0]
            data_final[in1, :, :] = ds1D[var_wnt].mean('time')
        da1_out = xr.DataArray(
            data_final, coords=[timefM.astype('datetime64[ns]'), lat_vec, lon_vec],
            dims=['time', 'lat_vec', 'lon_vec'])

        ds1_out = da1_out.to_dataset(name=var_wnt)

        # overwrite the existing file_in with this new updated dataset
        ds1.close()
        file_out = file_in
        ds1_out.to_netcdf(file_out)
    else:
        print('Dataset complete, no data downnloaded')

else:
    print('file not found: {}'.format(file_in))

# remove the directory and files that has the downloaded
shutil.rmtree('./data_x13/')

# list files
dir_list_end = os.listdir()
print("END -------------------------------")
print("Files and directories in  :")
print(dir_list_end)

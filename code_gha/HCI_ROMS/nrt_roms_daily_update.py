import calendar as clndr
import os
import urllib.request  # py36
import requests  # Aug 9, 2024 -- changed to request as urlib gave 404 error
import xarray as xr
import numpy as np
from scipy import interpolate


def fun_interp_NRT(ds_in, var_ds, z_int, intrp_mtrx):
    # Funtion to extrapolate data to 2 m using interpolate.interp1d.
    # Only grid points marked with a 1 in intrp_mtrx are interpolated.

    # num depths
    z1 = ds_in.z.data
    nd = len(z1)
    
    # size of input variables
    ny, nx = np.shape(intrp_mtrx)
    num_var_ds = len(var_ds)

    # interpolate ds_var at each grid point
    for j in range(0, num_var_ds):
        var1 = ds_in[var_ds[j]].data

        var_int = np.zeros([nd, ny, nx])
        for ii in range(0, ny):
            for jj in range(0, nx):
                if intrp_mtrx[ii, jj] == 1:
                    f = interpolate.interp1d(z1[0:-1], np.squeeze(var1[0:-1, ii, jj]), fill_value="extrapolate")
                    yf = f(z_int)
                else:
                    # yf = np.squeeze(ds1[var_ds[j]].data[:, ii, jj])
                    yf = var1[:, ii, jj]
                var_int[:, ii, jj] = yf

        # save as dataarray
        da1 = xr.DataArray(var_int, coords=[z1, lat_rho, lon_rho], dims=[
                           'depth', 'lat_rho', 'lon_rho'])

        # save as dataset
        if j == 0:
            ds_intrp = da1.to_dataset(name=var_ds[j])
        else:
            ds_intrp[var_ds[j]] = da1

    return ds_intrp


# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------
# ouput dataset name, file type, and directory
name_wnt = 'sal_temp_nrt'
file_type = 'nc'

# erddap url for all dates of the SST OI
# url = 'https://oceanmodeling.ucsc.edu/thredds/dodsC/ccsra_2016a_phys_agg_zlevs/fmrc/'
# fn1 = 'CCSRA_2016a_Phys_ROMS_z-level_(depth)_Aggregation_best.ncd'
url_fn1 = 'https://thredds.cencoos.org/thredds/dodsC/UCSC.nc?'

var_wnt = ['temp', 'salt']

# label of ROMS data
out_lbl = 'TS_nrt'

# file name of hist and nrt roms
fn1_in_hist = './data_gha/HCI_ROMS/sal_temp_hist.nc'
fn1_in_nrt = './data_gha/HCI_ROMS/sal_temp_nrt.nc'
# fn1_in_hist = './oc_cciea.github.io/data_gha/HCI_ROMS/sal_temp_hist.nc'
# fn1_in_nrt = './oc_cciea.github.io/data_gha/HCI_ROMS/sal_temp_nrt.nc'

# only download if there are these number of days available in the month
day_check = 24

# depth to interpolate
z_int = np.array([-250., -200., -150., -100.,  -80.,  -60.,  -40.,  -20., -10., -5., -2.])

# subset by distance from shore
xdis = 55

# dir out, will use artifacts to download this data
dir_out = './data_gha/HCI_ROMS/'
# dir_out = './data_x13/HCI_ROMS/'
dir_in = dir_out

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# len input variables
num_var_wnt = len(var_wnt)


dir_list = os.listdir()
print("START -------------------------------")
print("Files and directories in  :")
print(dir_list)
print("END -------------------------------")

# open link with NO SSL (security risk)
# url_fn1 = '{}{}'.format(url, fn1)
#store = xr.backends.PydapDataStore.open(url_fn1, verify=False)
#ds1 = xr.open_dataset(store)

ds1 = xr.open_dataset(url_fn1)


lat_rho = ds1['lat_rho'].data[:, 0]
lon_rho = ds1['lon_rho'].data[0, :]
depth = ds1['z'].data

ny = len(lat_rho)
nx = len(lon_rho)
nz = len(z_int)

# -- 2 m interpolate grid
# now create grid of positions that need to be interpolated, only
# needs to be done once
intrp_mtrx = np.zeros([ny, nx])
da_sample = ds1[var_wnt[0]][0, :, :, :]
data_sample = da_sample.data

for ii in range(ny):
    for jj in range(nx):
        # y1 = np.squeeze(ds1[var_wnt[0]][0, :, ii, jj]).data
        y1 = data_sample[:, ii, jj]
        ind = np.isfinite(y1)

        # if 2 m depth is nan (ie last z) then interpolate
        if not ind[-1]:
            intrp_mtrx[ii, jj] = 1
        if len(ind.nonzero()[0]) == 0:
            intrp_mtrx[ii, jj] = 0

# -- xdis from shore grid
data_sample = da_sample.sum('z').data

# get lat and lon
ny_rng = ny
nx_rng = xdis
lon_mtrx = np.zeros([ny_rng, nx_rng])
lat_mtrx = np.zeros([ny_rng, nx_rng])
y_wnt = np.zeros([ny_rng, nx_rng])
x_wnt = np.zeros([ny_rng, nx_rng])

# find shore and get the lon points
for i in range(ny_rng):
    in_lat = i
    data1 = np.squeeze(data_sample[in_lat, :])
    # can have islands, calculate cumsum and find first nan
    data1_cs = np.nancumsum(data1)
    ind_nan = np.where(data1_cs == data1_cs[-1])[0]
    in2 = ind_nan[0] + 1
    in1 = in2 - xdis
    lon1 = lon_rho[in1:in2]
    lat1 = np.ones(xdis)*lat_rho[i]
    lon_mtrx[i, :] = lon1
    lat_mtrx[i, :] = lat1
    y_wnt[i, :] = in_lat
    x_wnt[i, :] = np.arange(in1, in2)


# Get the time available in the NRT ROMS data product
time_nrt = ds1.time.data.astype('datetime64[D]')
yy_nrt = ds1.time.dt.year.data
mm_nrt = ds1.time.dt.month.data

yr_end_nrt = ds1.time.dt.year.data[-1]
mon_end_nrt = ds1.time.dt.month.data[-1]
day_end_nrt = ds1.time.dt.day.data[-1]

date2 = np.datetime64('{}-{:02d}'.format(yr_end_nrt, mon_end_nrt), 'M')

# Open the hist roms data
if os.path.isfile(fn1_in_nrt):
    print("1, File exists: {}".format(fn1_in_nrt))
    ds2 = xr.open_dataset(fn1_in_hist)


# Check if the file exists and is a file
if os.path.isfile(fn1_in_nrt):
    print("1, File exists: {}".format(fn1_in_nrt))
    ds3 = xr.open_dataset(fn1_in_nrt)

    yr_end = ds3.time.dt.year.data[-1]
    mon_end = ds3.time.dt.month.data[-1]
    date1 = np.datetime64('{}-{:02d}'.format(yr_end, mon_end), 'M')

    # make the output directory to save monthly and daily datasets
    try:
        os.makedirs(dir_out)
    except OSError:
        if not os.path.isdir(dir_out):
            raise

    # check to see if number of days downloaded for the last month is not complete
    # redownload the last month if it isn't, otherwise skip to the next month
    num_day_miss_end = ds3.day_missing.data[-1]
    if num_day_miss_end > 0:
        dates_wnt = np.arange(date1, date2+1)
    else:
        dates_wnt = np.arange(date1+1, date2+1)

    # loop over the dates and save the daily data for each month
    ntM = len(dates_wnt)
    date_final = dates_wnt[0]
    ds_append = []
    date_append = []
    num_day_diff_append = []
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

        num_day_clndr = clndr.monthrange(yr1, mon1)[1]
        
        print('yr1={}, mon1={}, day1={}, yr2={}, mon2={}, day2={}'.format(yr1, mon1, day1, yr2, mon2, day2))
        print('num_days_available={}, day_check={}'.format(num_days_available, day_check))
        if num_days_available > day_check:
            date_final = dates_wnt[i]
            date_append.append(date_final)
            num_day_diff = num_day_clndr - day2
            num_day_diff_append.append(num_day_diff)

            # --create output directory for daily data
            dir1 = '{}{}/{:02d}/'.format(dir_out, yri, moni)

            # --check if directory exist, if it doesn't then create
            try:
                os.makedirs(dir1)
            except OSError:
                if not os.path.isdir(dir1):
                    raise

            # download the data
            da1yr = ds1.sel(time=ds1.time.dt.year.isin(yri)).time
            da1mon = da1yr.sel(time=da1yr.time.dt.month.isin(moni)).time
            ds1i = ds1.sel(time=ds1.time.isin(da1mon.time.data))
   
            # extract the dataarray
            for j in range(num_var_wnt):
                da1j = ds1i[var_wnt[j]].squeeze()

                # create dataset
                if j == 0:
                    ds1_out = da1j.to_dataset(name=var_wnt[j])
                else:
                    ds1_out[var_wnt[j]] = da1j

            # save to final netcdf file
            fn_out_daily = '{}{}_daily.{}'.format(dir1, out_lbl, file_type)
            print('Daily Download and Save, i={}, {}'.format(i, fn_out_daily))
            ds1_out.to_netcdf(fn_out_daily)

            # 2 m interp
            timeD = ds1_out.time.data
            ntD = len(timeD)

            data_intrp_mtrx = np.zeros([num_var_wnt, ntD, nz, ny, nx])
            for iii in range(ntD):
                ds1_in = ds1_out.sel(time=timeD[iii])
                da_intrp = fun_interp_NRT(ds1_in, var_wnt, z_int, intrp_mtrx)

                for jjj in range(num_var_wnt):
                    data_intrp_mtrx[jjj, iii, :, :, :] = da_intrp[var_wnt[jjj]].data

            eta_rho = ds1_out.eta_rho.data
            xi_rho = ds1_out.xi_rho.data
            for jjj in range(num_var_wnt):
                da1_out = xr.DataArray(data_intrp_mtrx[jjj, :, :, :, :], coords=[timeD, z_int, eta_rho, xi_rho], dims=['time', 'z', 'eta_rho', 'xi_rho'])
                if jjj == 0:
                    ds1_out_intrp = da1_out.to_dataset(name=var_wnt[jjj])
                else:
                    ds1_out_intrp[var_wnt[jjj]] = da1_out

            # monthly means of the 2 m interp
            yr_intrp = ds1_out_intrp.time.dt.year.data[0]
            ds1M = ds1_out_intrp.groupby('time.month').mean('time')
            ds1M.attrs['year'] = yr_intrp

            ds_append.append(ds1M)

    varM_mtrx = np.zeros([num_var_wnt, len(ds_append), nz, ny_rng, nx_rng])
    for i in range(len(ds_append)):
        for j in range(ny_rng):
            inx = x_wnt[j, :].astype('int')
            iny = y_wnt[j, :].astype('int')

            # save to matrix in order to calculate monthly mean
            for k in range(num_var_wnt):
                varM_mtrx[k, i, :, j, :] = ds_append[i][var_wnt[k]][0, :, iny[0], inx]

    # create nd.dataarray
    lat_indx = np.arange(0, ny_rng)
    lon_indx = np.arange(0, nx_rng)
    time_out = np.array(date_append).astype('datetime64[ns]')
    num_day_missing = np.array(num_day_diff_append)

    da1 = xr.DataArray(lon_mtrx, coords=[lat_indx, lon_indx],
                       dims=['lat_indx', 'lon_indx'])
    da2 = xr.DataArray(lat_mtrx, coords=[lat_indx, lon_indx],
                       dims=['lat_indx', 'lon_indx'])
    da3 = xr.DataArray(num_day_missing, coords=[time_out], dims=['time'])

    ds_out = da1.to_dataset(name='lon_mtrx')
    ds_out['lat_mtrx'] = da2
    ds_out['day_missing'] = da3

    for i in range(num_var_wnt):
        da1i = xr.DataArray(varM_mtrx[i, :, :, :].squeeze(),
                            coords=[time_out, depth, lat_indx, lon_indx],
                            dims=['time', 'depth', 'lat_indx', 'lon_indx'])

        ds_out['{}_mtrx'.format(var_wnt[i])] = da1i

    # filename out
    fn_out = '{}sal_temp_nrt.nc'.format(dir_out)
    ds_out.to_netcdf(fn_out)

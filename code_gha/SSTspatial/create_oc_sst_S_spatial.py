import os
import numpy as np
import xarray as xr
import pandas as pd
from C_iea import C_iea


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
iea_yr = 2025

fn_in = 'TS_monthly.nc'

ds1_dim = ['lat_vec', 'lon_vec', 'time']
ds1_var = ['sst']

lon_bgn = 195.125
lon_end = 244.625
lat_bgn = 30.125
lat_end = 60.625

dx = 1.0

freq_wnt = 'YS'
season = [3, 6, 9, 12]
season_lbl = ['winter', 'spring', 'summer', 'fall']

wndw = 5
yy_end = iea_yr
yr_clim_bgn = 1982
yr_clim_end = yy_end

att_coord = ['lon', 'lat']
data_wnt = ['anom_end', 'mn5_sd', 'trnd5_sd']
mrkr_wnt = ['map_marker_anom', 'map_marker_mn5', 'map_marker_trnd5']

dir_out = './data_gha/SSTspatial/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

# length of input variables
num_season = len(season)

# open the data
ds1 = xr.open_dataset(fn_in)
lat = ds1[ds1_dim[0]].data
lon = ds1[ds1_dim[1]].data
time = ds1[ds1_dim[2]].data

# print the time dates
print(ds1.time.data)

# years
yrs = np.unique(ds1['time'].dt.year.data)
dt_yrs = pd.to_datetime(yrs)
num_yrs = len(yrs)

# setup a yearly pd datetime index
dt_yrs = pd.to_datetime(
    {'year': yrs, 'month': np.ones(len(yrs)), 'day': np.ones(len(yrs))})

# create lon and lat grid for spatial averages
lon_wnt = np.arange(lon_bgn, lon_end+dx, dx)
lat_wnt = np.arange(lat_bgn, lat_end+dx, dx)

long, latg = np.meshgrid(lon_wnt, lat_wnt)
ny, nx = long.shape

# create spatial means
mrkr = 0
iea_list = list()
dataS_mtrx = np.zeros([ny, nx, num_yrs, num_season])*np.nan
for i in range(ny):
    lat_bttm = lat_wnt[i]-dx/2.0
    lat_top = lat_wnt[i]+dx/2.0
    in_lat = np.logical_and(lat > lat_bttm, lat < lat_top)

    lat_cntr = np.nanmean(lat[in_lat])
    for j in range(nx):
        lon_bttm = lon_wnt[j]-dx/2.0
        lon_top = lon_wnt[j]+dx/2.0
        in_lon = np.logical_and(lon > lon_bttm, lon < lon_top)
        lon_cntr = np.nanmean(lon[in_lon])

        data_box = ds1[ds1_var[0]][:, in_lat, in_lon].mean(
            axis=(1, 2), skipna=True).data

        # monthly values of spatially averaged data in box
        # ps1M = pd.Series(data=data_box, index=time).asfreq(freq_wnt)
        ps1M = pd.Series(data=data_box, index=time)

        # seasonal values of spatially averaged data in box
        # 03=Jan-Mar, 06=Apr-Jun, 09=Jul-Sep, 12=Oct-Dec
        # ps1S = ps1M.resample('Q-Mar').mean()
        ps1S = ps1M.resample('QE-MAR').mean()

        for k in range(num_season):
            in_ssn = np.where(ps1S.index.month == season[k])
            yrs_ssn = ps1S.index.year.values[in_ssn]
            data_ssn = ps1S.values[in_ssn]
            ia = np.isin(yrs, yrs_ssn)
            ib = np.isin(yrs_ssn, yrs)
            dataS_mtrx[i, j, ia, k] = data_ssn[ib]

# make the output directory to save monthly and daily datasets
try:
    os.makedirs(dir_out)
except OSError:
    if not os.path.isdir(dir_out):
        raise

# netcdf file of seasonal IEA stats for all grid locations
for i in range(num_season):
    print('season: {}'.format(i+1))
    mrkr = 0
    # iea_list = list()
    ts_mtrx = np.zeros([ny*nx, num_yrs])
    coord_mtrx = np.zeros([ny*nx, len(att_coord)])
    data_mtrx = np.zeros([ny*nx, len(data_wnt)])
    mrkr_mtrx = np.zeros([ny*nx, len(data_wnt)], dtype='str')
    for j in range(ny):
        for k in range(nx):
            data1 = np.squeeze(dataS_mtrx[j, k, :, i])
            ps1 = pd.Series(data=data1, index=dt_yrs).asfreq(freq_wnt)
            z1 = C_iea(ps1, yr_clim_bgn=yr_clim_bgn, yr_clim_end=yr_clim_end,
                       wndw=wndw, yy_end=yy_end,
                       lon=long[j, k], lat=latg[j, k])

            # place in matrices
            ts_mtrx[mrkr, :] = data1

            for iii in range(len(att_coord)):
                coord_mtrx[mrkr, iii] = getattr(z1, att_coord[iii])

            for iii in range(len(data_wnt)):
                data_mtrx[mrkr, iii] = getattr(z1, data_wnt[iii])

            for iii in range(len(mrkr_wnt)):
                mrkr_mtrx[mrkr, iii] = getattr(z1, mrkr_wnt[iii])

            mrkr = mrkr + 1

    # put into xr.da
    loc_indx = np.arange(0, mrkr)
    da1 = xr.DataArray(coord_mtrx, coords=[loc_indx, att_coord],
                       dims=['index', 'coord'])
    da2 = xr.DataArray(data_mtrx, coords=[loc_indx, data_wnt],
                       dims=['index', 'data'])
    da3 = xr.DataArray(mrkr_mtrx, coords=[loc_indx, mrkr_wnt],
                       dims=['index', 'marker'])
    da4 = xr.DataArray(ts_mtrx, coords=[loc_indx, dt_yrs],
                       dims=['index', 'time'])

    # put into xr.ds
    ds1 = da1.to_dataset(name='coord_mtrx')
    ds1['data_mtrx'] = da2
    ds1['mrkr_mtrx'] = da3
    ds1['ts_mtrx'] = da4

    fn_out_clim = '{}/anom_mn5_trnd5_{}_clim_{}_{}.nc'.format(
        dir_out, season_lbl[i], yr_clim_bgn, yr_clim_end)

    # --Save Dataset to a netcdf file
    ds1.to_netcdf(fn_out_clim)

dir_list_end = os.listdir()
print("END -------------------------------")
print("Files and directories in  :")
print(dir_list_end)

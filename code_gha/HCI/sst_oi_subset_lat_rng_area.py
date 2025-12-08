import os
import xarray as xr
import numpy as np


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# 1) Distance from shore netcdf data done the correct way
# distance
dis_wnt = [75, 150]

# 2) OI SST monthly means
fn_sst = '~TS_monthly.nc'

# 3)
# x distance and lat range list
xdis_km = [75, 150]

test_rgn = 0
if test_rgn == 0:
    lat_rgn = [[43.5, 48], [40, 43.5], [35.5, 40], [30, 35.5]]
else:
     lat_rgn = [[30, 48]]

var_wnt = 'sst_oi'

# 4) ouput directory


# 5) dir input
dir_in = dir_out
# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# input variable size
num_xdis_km = len(xdis_km)
num_lat_rgn = len(lat_rgn)

# A) Get the distance from shore data
dis_str = '_'.join(list(map(str,  dis_wnt)))
fn_in = '{}{}_distance_to_shore_dis_{}km.nc'.format(dir_in, var_wnt, dis_str)
ds1 = xr.open_dataset(fn_in)

# 2) Open SST OI
ds1_sst = xr.open_dataset(fn_sst)
da1 = ds1_sst['sst']

time1 = da1.time.data
ntf = time1.shape[0]

# Subset by latitude region
lat_vec = da1['lat_vec'].data
lon_vec = da1['lon_vec'].data - 360

lon_mtrx, lat_mtrx = np.meshgrid(lon_vec, lat_vec)


# Use the mask_mtrx to get the locations of distance from shore wanted.
da1_mask = ds1['mask_mtrx']
da1_area = ds1['area_mtrx']
distance = da1_mask['distance'].data


# subset by lat and distance to shore
for i in range(num_lat_rgn):
    # lat range
    lat_rgn1 = lat_rgn[i]
    in_lat_rgn1 = np.logical_and(lat_vec >= lat_rgn1[0], lat_vec <= lat_rgn1[1])
    da1_lat = da1[:, in_lat_rgn1, :]
    lat_mtrx1 = lat_mtrx[in_lat_rgn1, :]
    lon_mtrx1 = lon_mtrx[in_lat_rgn1, :]

    for j in range(num_xdis_km):
        xdis1 = xdis_km[j]

        in_dis = np.where(distance == xdis1)[0]
        da1_mask_lat = np.squeeze(da1_mask[:, :, in_dis].sel(
            latitude=slice(lat_rgn1[0], lat_rgn1[1])))
        da1_area_lat = np.squeeze(da1_area[:, :, in_dis].sel(
            latitude=slice(lat_rgn1[0], lat_rgn1[1])))

        ny1, nx1 = da1_mask_lat.shape

        data_mtrx = np.zeros([ntf, ny1, nx1])*np.nan
        area_mtrx = np.zeros([ny1, nx1])*np.nan
        for iii in range(ny1):
            for jjj in range(nx1):
                mask1 = da1_mask_lat.data[iii, jjj].astype('int')

                if mask1 == 1:
                    in_lon1 = np.where(lon_mtrx1[iii, :] == da1_mask_lat['longitude'].data[jjj])[0]
                    data1 = np.squeeze(da1_lat.data[:, iii, in_lon1])
                    data_mtrx[:, iii, jjj] = data1

                    area1 = da1_area_lat[iii, jjj]
                    area_mtrx[iii, jjj] = area1

        # xarray dataarray and dataset
        da1_data_out = xr.DataArray(data_mtrx, coords=[time1.astype('datetime64'), da1_mask_lat['latitude'].data, da1_mask_lat['longitude'].data], dims=['time', 'latitude', 'longitude'])
        da1_area_out = xr.DataArray(area_mtrx, coords=[da1_mask_lat['latitude'].data, da1_mask_lat['longitude'].data], dims=['latitude', 'longitude'])

        ds1_out = da1_data_out.to_dataset(name=var_wnt)
        ds1_out['area'] = da1_area_out

        # save to netcdf
        fn_out = '{}{}_lat_{}_{}_xdis_{}km.nc'.format(dir_out, var_wnt, lat_rgn1[0], lat_rgn1[1], xdis1)
        ds1_out.to_netcdf(fn_out)

# remove some large files that can not be commit to github
os.remove("TS_monthly.nc")

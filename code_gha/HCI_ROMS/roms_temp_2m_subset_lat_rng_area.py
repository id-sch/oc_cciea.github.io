import os
import xarray as xr
import numpy as np

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# 1) Distance from shore netcdf data done the correct way
# distance
dis_wnt = [75, 150]

# 2) ROMS 2m temp monthly means
# fn_nrt_in = './oc_cciea.github.io/data_gha/HCI_ROMS/sal_temp_nrt.nc'
# fn_hist_in = './oc_cciea.github.io/data_gha/HCI_ROMS/sal_temp_hist.nc'
fn_nrt_in = './data_gha/HCI_ROMS/sal_temp_nrt.nc'
fn_hist_in = './data_gha/HCI_ROMS/sal_temp_hist.nc'

# 3)
# x distance and lat range list
xdis_km = [75, 150]
# fn_dis_shore_in = './data_x13/HCI_ROMS/roms_distance_to_shore_dis_75_100_150_175_200_225_250_275_300_325_350km.nc'
fn_dis_shore_in = './data_gha/HCI_ROMS/HCI_ROMS/roms_distance_to_shore_dis_75_100_150_175_200_225_250_275_300_325_350km.nc'

lat_rgn = [[43.5, 48], [40, 43.5], [35.5, 40], [30, 35.5], [30, 48]]

var_wnt = 'temp_mtrx'

z_wnt = -2

# 4) ouput directory
dir_out = './data_gha/HCI_ROMS/'
# dir_out = './data_x13/HCI_ROMS/'

# 5) dir input
dir_in = dir_out
# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

dir_list = os.listdir()
print("START -------------------------------")
print("Files and directories in  :")
print(dir_list)


# input variable size
num_xdis_km = len(xdis_km)
num_lat_rgn = len(lat_rgn)

# A) Get the distance from shore data
ds1 = xr.open_dataset(fn_dis_shore_in)

# 2) Open temp 2 m ROMS, hist and nrt
ds1_hist = xr.open_dataset(fn_hist_in)
ds1_nrt = xr.open_dataset(fn_nrt_in)

# 2) get var want and 
# get lat and lon matrix (same for both hist or nrt)
lon_mtrx = ds1_nrt['lon_mtrx'].data
lat_mtrx = ds1_nrt['lat_mtrx'].data

# index for depth wanted
in_z = np.where(ds1_nrt['depth'].data == z_wnt)[0]

# create temp dataarray by concat hist and nrt, for var_wnt and z_wnt
da1 = np.squeeze(xr.concat(
    (ds1_hist[var_wnt][:, in_z, :, :], ds1_nrt[var_wnt][:, in_z, :, :]),
    'time'))
time1 = da1.time.data
ntf = time1.shape[0]

# Subset by latitude region, however lat in not a variable in the dataarray.
# Instead find index of lats using the lat_mtrx
lat_vec = lat_mtrx[:, 0]

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


dir_list = os.listdir()
print("END -------------------------------")
print("Files and directories in  :")
print(dir_list)


# # remove some large files that can not be commit to github
# os.remove("TS_monthly.nc")


dir_list = os.listdir()
print("END: removed large file  -------------------------------")
print("Files and directories in  :")
print(dir_list)

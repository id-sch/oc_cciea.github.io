import os
import xarray as xr
import pandas as pd
import numpy as np
from fun_hampel_outlier import hampel


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# --input directory
dir_in = './data_gha/CalCofiCSV/'
# dir_in = './oc_cciea.github.io/data_gha/CalCofiCSV/'

# station file of 66 CalCOFI stations
fn_sttn = '{}calcofi_stations66.csv'.format(dir_in)
# fn_sttn = './data_x13/CalCofiCSV/calcofi_stations66.csv'
var_fn_sttn = ['Line', 'Station', 'Lat', 'Lon', 'Distance', 'Depth']


# --input file type extension
file_type_in = 'nc'

# --names of salinity and temperature in the netcdf file
var_wnt = ['R_O2']
num_var = len(var_wnt)

# --names of time, depth, lat, lon dimensions
var_dpth = 'depth'
var_time = 'time'
var_lat = 'Lat_Dec'
var_lon = 'Lon_Dec'

# --For "Step 1" outliers removed by Hampel filter that requires 3 inputs
# --input1 = lenght of window around sample point, eg 7 will have 3
# --points before and 3 points after
# --input2 = number of standard deviation threshold
# --input3 = if set then a cast with only one value will be set to have
#            no values, this is needed if you are interpolating the output
hf_in1 = 11
hf_in2 = 4
hf_intrp_flag = 1

# missing DO value (summer 2020 86.7-45 cast has values of 0 for much of the cast, remove these)
missing_value = 0

# --Since the Hample filter uses "pandas series" we can easily interpolate
# --with different types of interpolation methods built into pandas
intrp_mthd = 'slinear'

# --names of output vectors, will have dimension of "nt"
dpth_wnt = np.array([50, 150, 500])
num_dpth = len(dpth_wnt)
out_vec_wnt = [var_wnt[0]+'{}m'.format(x) for x in dpth_wnt]
num_out_vec = len(out_vec_wnt)

# --names of output matrices, will have dimension of "nd x nt"
out_mtrx_wnt = var_wnt
num_out_mtrx = len(out_mtrx_wnt)

# --output filename and extention
file_type_out = 'nc'
file_name_out = 'do'

# depth range to find min value, data stored in min_mtrx and z_min_mtrx
z1 = 0
z2 = 500

# dir out
dir_out = './data_x13/DissolvedOxygen/'

# --END: Change These
# ----------------------------------------------------------------------

# stations
df_sttn = pd.read_csv(fn_sttn)
num_sttn = len(df_sttn)
sttn_wnt = []
for i in range(num_sttn):
    line1 = df_sttn[var_fn_sttn[0]].values[i]*10
    sttn1 = df_sttn[var_fn_sttn[1]].values[i]*10
    sttn_wnt.append('{:3d}_{:3d}'.format(
        line1.astype('int'), sttn1.astype('int')))

# bottom depth
bottom_depth = df_sttn[var_fn_sttn[5]].values

# loop over all stations
for i in range(num_sttn):
    # bottom_depth at station
    bottom_depth_station = bottom_depth[i]

    # cast only go to 500m, if bottom depth is greater then set to 500m
    if bottom_depth_station > 500:
        bottom_depth_station = 500

    # set range that cast has to be inorder for it to be considered bottom
    bottom_check = bottom_depth_station - bottom_depth_station*0.3

    # --create input filename
    fn_in = dir_in + '/' + sttn_wnt[i] + '.' + file_type_in

    # --open the netcdf files as an xarray
    dall = xr.open_dataset(fn_in)

    # --extract var_wnt as DataArray, temp=da0, sal=da1
    da0 = dall[var_wnt[0]]

    # --get the time and depth dimension
    z = dall[var_dpth].data
    tt = dall[var_time].data
    lat = dall[var_lat].data
    lon = dall[var_lon].data
    station_lat = np.nanmean(lat)
    station_lon = np.nanmean(lon)

    # get index of depth values between z1 and z2
    inz = np.logical_and(z >= z1, z <= z2)

    # --get dimensions of time and depth
    nt = len(dall[var_time])
    nd = len(dall[var_dpth])

    # initialize arrays for output vectors and matrices
    out_vec = np.zeros([num_out_vec, nt])*np.nan
    out_mtrx = np.zeros([num_out_mtrx, nd, nt])*np.nan
    min_mtrx = np.zeros([num_out_mtrx, nt])*np.nan
    z_min_mtrx = np.zeros([num_out_mtrx, nt])*np.nan
    bottom_mtrx = np.zeros([num_out_mtrx, nt])*np.nan

    for j in range(0, nt):
        # ---------------------------------------------------------------------
        # --Filter each profile by: 1) removing outliers
        # --------------------------------------------------------------------
        # (1) Apply Hampel filter to remove outliers, this needs "pandas
        # series" that can be created by xr.DataArray.to_pandas
        s0 = xr.DataArray.to_pandas(da0[:, j])
        s_hf0 = hampel(s0, hf_in1, hf_in2, hf_intrp_flag).interpolate(
            method=intrp_mthd)
        da0[:, j] = s_hf0.values

        # chang missing values to np.nan
        data_hampel = np.copy(s_hf0.values)
        ind_miss = np.where(data_hampel == missing_value)[0]
        data_hampel[ind_miss] = np.nan

        # --place the filtered data in a matrix
        out_mtrx[:, :, j] = data_hampel

        # --get data at dpth_wnt, can be NaN if cast is not as
        # --deep as depth_wnt
        ia = np.isin(da0.depth, dpth_wnt)
        ib = np.isin(dpth_wnt, da0.depth)
        # data_dpth = s_hf0.values[ia]
        data_dpth = data_hampel[ia]
        out_vec[ib, j] = data_dpth

        # check to see if cast is all nan
        data_z1_z2 = data_hampel[inz]
        z1_z2 = z[inz]
        num_nan = len(np.isnan(data_z1_z2).nonzero()[0])
        num_all = len(z1_z2)

        # find missing values in the cast
        ind_z1_z2 = np.isfinite(data_z1_z2).nonzero()[0]
        numd = len(ind_z1_z2)

        # get min and bottom data if there is data
        if numd > 0:
            # find min value between z1 and z2
            in_min = np.nanargmin(data_z1_z2)
            data_min = data_z1_z2[in_min]
            z_data_min = z1_z2[in_min]

            # bottom
            z1_z2d = z1_z2[ind_z1_z2]
            data_z1_z2d = data_z1_z2[ind_z1_z2]
            z_bottom = z1_z2d[-1]

            if z_bottom > bottom_check:
                data_bottom = data_z1_z2d[-1]
            else:
                data_bottom = np.nan
        else:
            data_min = np.nan
            z_data_min = np.nan
            data_bottom = np.nan

        min_mtrx[:, j] = data_min
        z_min_mtrx[:, j] = z_data_min
        bottom_mtrx[:, j] = data_bottom

    # --create DataArray from the out_mtrx_wnt, and add to Dataset, but first we need to create it,
    for j in range(0, num_out_mtrx):
        mtrxj = np.squeeze(out_mtrx[j, :, :])
        ma1 = xr.DataArray(mtrxj, coords=[z, tt], dims=['depth', 'time'])
        # ds1[out_mtrx_wnt[j]] = da1
        # --otherwise just add additional varialbes to the Dataset
        if j == 0:
            ds1 = ma1.to_dataset(name=out_mtrx_wnt[j])
        else:
            ds1[out_mtrx_wnt[j]] = ma1

    # --create DataArray from the min_mtrx_wnt, and add to Dataset
    for j in range(0, num_out_mtrx):
        min_mtrx_wnt = '{}_min'.format(out_mtrx_wnt[j])
        mtrxj = np.squeeze(min_mtrx[j, :])
        ma1 = xr.DataArray(mtrxj, coords=[tt], dims=['time'])
        # ds1[out_mtrx_wnt[j]] = da1
        # --otherwise just add additional varialbes to the Dataset
        ds1[min_mtrx_wnt] = ma1

    # --create DataArray from the z_min_mtrx_wnt, and add to Dataset
    for j in range(0, num_out_mtrx):
        z_min_mtrx_wnt = '{}_z_min'.format(out_mtrx_wnt[j])
        mtrxj = np.squeeze(z_min_mtrx[j, :])
        ma1 = xr.DataArray(mtrxj, coords=[tt], dims=['time'])
        # ds1[out_mtrx_wnt[j]] = da1
        # --otherwise just add additional varialbes to the Dataset
        ds1[z_min_mtrx_wnt] = ma1

    # --create DataArray from the bottom_mtrx_wnt, and add to Dataset
    for j in range(0, num_out_mtrx):
        bottom_mtrx_wnt = '{}_bottom'.format(out_mtrx_wnt[j])
        mtrxj = np.squeeze(bottom_mtrx[j, :])
        ma1 = xr.DataArray(mtrxj, coords=[tt], dims=['time'])
        # ds1[out_mtrx_wnt[j]] = da1
        # --otherwise just add additional varialbes to the Dataset
        ds1[bottom_mtrx_wnt] = ma1

    # --create DataArray from the out_vec_wnt
    for j in range(0, num_out_vec):
        vecj = np.squeeze(out_vec[j, :])
        va1 = xr.DataArray(vecj, coords=[tt], dims=['time'])
        # --add this new DataArray to the Dataset
        ds1[out_vec_wnt[j]] = va1

    # --add quarter DataArray to the Dataset
    ds1['Quarter'] = dall['Quarter']

    # --add the station name, lat and lon of the station as an attribute
    ds1.attrs['station'] = sttn_wnt[i]
    ds1.attrs['lat'] = station_lat
    ds1.attrs['lon'] = station_lon
    ds1.attrs['sttn_depth'] = bottom_depth[i]

    # -------------------------------------------------------

    # --check if directory exist, if it doesn't then create
    try:
        os.makedirs(dir_out)
    except OSError:
        if not os.path.isdir(dir_out):
            raise

    # --create output directory
    dir_fn = '{}/{}_{}_min_between_{}_{}m.nc'.format(
        dir_out, file_name_out, sttn_wnt[i], z1, z2)

    # --Save Dataset to a netcdf file
    ds1.to_netcdf(dir_fn)

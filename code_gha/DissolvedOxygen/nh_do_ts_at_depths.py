import re
import os
import xarray as xr
import numpy as np
from fun_hampel_outlier import hampel


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# --station wanted, two stations newport hydrographic at 5 and 25 km
sttn_wnt = ['NH05', 'NH25']
num_sttn = len(sttn_wnt)

# --input directory
dir_in = './'

# --input files were created from excel files, two of them
#   so specifiy which one
file_excel_in = '_CTD'

# --file input extension
file_type_in = '.nc'

# --names of salinity and temperature in the netcdf file
var_wnt = ['Oxygen']
num_var = len(var_wnt)

# --names of time, depth, lat, lon dimensions
var_dpth = 'depth'
var_time = 'time'
var_lat = 'Lat'
var_lon = 'Long'

# --For "Step 1" outliers removed by Hampel filter that requires 3 inputs
# --input1 = lenght of window around sample point, eg 7 will have 3
# --points before and 3 points after
# --input2 = number of standard deviation threshold
# --input3 = if set then a cast with only one value will be set to have
#            no values, this is needed if you are interpolating the output
hf_in1 = 11
hf_in2 = 4
hf_intrp_flag = 1

# --Since the Hample filter uses "pandas series" we can easily interpolate
# --with different types of interpolation methods built into pandas
intrp_mthd = 'slinear'

# --names of output vectors, will have dimension of "nt"
dpth_wnt = np.array([50, 150])
num_dpth = len(dpth_wnt)
out_vec_wnt = [var_wnt[0]+'{}m'.format(x) for x in dpth_wnt]
num_out_vec = len(out_vec_wnt)

# --names of output matrices, will have dimension of "nd x nt"
out_mtrx_wnt = var_wnt
num_out_mtrx = len(out_mtrx_wnt)

# --output filename and extention
file_type_out = 'nc'
file_name_out = 'do_ts'

# dir out
dir_out = './data_x13/DissolvedOxygen/'

# --END: Change These
# ----------------------------------------------------------------------

for i in range(0, num_sttn):
    # --create input filename
    fn_in = dir_in + '/' + sttn_wnt[i] + file_excel_in + file_type_in

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

    # --get dimensions of time and depth
    nt = len(dall[var_time])
    nd = len(dall[var_dpth])

    # initialize arrays for output vectors and matrices
    out_vec = np.zeros([num_out_vec, nt])*np.nan
    out_mtrx = np.zeros([num_out_mtrx, nd, nt])*np.nan

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

        # --place the filtered data in a matrix
        out_mtrx[:, :, j] = s_hf0.values

        # --get data at dpth_wnt, can be NaN if cast doesn't in not as
        # --deep as depth_wnt
        ia = np.isin(da0.depth, dpth_wnt)
        ib = np.isin(dpth_wnt, da0.depth)
        data_dpth = s_hf0.values[ia]
        out_vec[ib, j] = data_dpth

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

    # --create DataArray from the out_vec_wnt
    for j in range(0, num_out_vec):
        vecj = np.squeeze(out_vec[j, :])
        va1 = xr.DataArray(vecj, coords=[tt], dims=['time'])
        # --add this new DataArray to the Dataset
        ds1[out_vec_wnt[j]] = va1

    # --add the station name, lat and lon of the station as an attribute
    ds1.attrs['station'] = sttn_wnt[i]
    ds1.attrs['lat'] = station_lat
    ds1.attrs['lon'] = station_lon

    # -------------------------------------------------------
    # --create output directory
    fn_out = re.sub(r'[^\w]', '', file_name_out + '_' +
                    sttn_wnt[i]) + '.' + file_type_out
    dir_fn = dir_out + '/' + fn_out

    # --check if directory exist, if it doesn't then create
    try:
        os.makedirs(dir_out)
    except OSError:
        if not os.path.isdir(dir_out):
            raise

    # --Save Dataset to a netcdf file
    ds1.to_netcdf(dir_fn)

# remove NH25 netcdf, keep this private
os.remove("./NH05_CTD.nc")
os.remove("./NH25_CTD.nc")

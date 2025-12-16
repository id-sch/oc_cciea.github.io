import os
import xarray as xr
import numpy as np
import pandas as pd
from C_iea_calcofi import C_iea


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------

# IEA year
iea_yr = 2025

# --input directory
dir_in = './data_x13/DissolvedOxygen/'

# station file of 66 CalCOFI stations
dir_sttn_in = './oc_cciea.github.io/data_gha/CalCofiCSV/'
fn_sttn = '{}calcofi_stations66.csv'.format(dir_sttn_in)
var_fn_sttn = ['Line', 'Station', 'Lat', 'Lon', 'Distance']

# variable wanted
var_fn_wnt = 'do'

# calcofi variable name
var_calcofi_wnt = 'R_O2'

# depth of CalCOFI variables wanted
dpth_wnt = [50, 150, 500]

# min value between these depths
z1 = 0
z2 = 500

# --name of useful attributes
att_wnt = ['lat', 'lon', 'station']
num_att = len(att_wnt)

# --input file type extension
file_type_in = 'nc'
file_name_in = 'do'

# --IEA file names
file_pre = 'oc_do'

# --IEA labels
y_lbl = 'DO (ml/L)'
ttl_lbl = 'Dissolved Oxygen'

# --names of time, depth, lat dimensions
var_dpth = 'depth'
var_time = 'time'
var_qrtr = 'Quarter'

# --quarter names
ssn_lbl = ['Winter', 'Spring', 'Summer', 'Fall']
# ssn_lbl = ['Winter']

# pandas frequency for time series
# freq_wnt = 'AS'
freq_wnt = 'YS'

# iea window
wndw = 5
yy_end = iea_yr

yr_clim_bgn = 1984
yr_clim_end = yy_end

# ds labels
att_coord = ['lon', 'lat']
att_sttn = ['line', 'sttn']
data_wnt = ['data_end', 'mn5_sd', 'trnd5_sd']

# names of variables in the IEA matrix, these mark if present value is
# higher/lower than 1sd, mean 5 is higher/lower than 1sd, and if trend 5 is
# higher/lower than 1sd
mrkr_wnt = ['map_marker_anom', 'map_marker_mn5', 'map_marker_trnd5']

# output directory
dir_out = './data_x13/DissolvedOxygen/'
# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# input variable sizes
num_dpth = len(dpth_wnt)
num_quarter = len(ssn_lbl)

# variable want and label are kind of complicated, combination of dpth_wnt, min, z_min, bottom, z1-z2
var_wnt = list()
var_lbl_wnt = list()
for i in range(num_dpth):
    var_wnt.append('{}{}m'.format(var_calcofi_wnt, dpth_wnt[i]))
#    var_lbl_wnt.append('at {} m'.format(dpth_wnt[i]))
var_wnt.append('{}_min'.format(var_calcofi_wnt))
var_wnt.append('{}_z_min'.format(var_calcofi_wnt))
var_wnt.append('{}_bottom'.format(var_calcofi_wnt))
# var_lbl_wnt.append('min value between {}-{} m'.format(z1, z2))
# var_lbl_wnt.append('depth of min value between {}-{} m'.format(z1, z2))

num_var = len(var_wnt)

# stations
df_sttn = pd.read_csv(fn_sttn)
num_sttn = len(df_sttn)
sttn_wnt = list()
for i in range(num_sttn):
    line1 = df_sttn[var_fn_sttn[0]].values[i]*10
    sttn1 = df_sttn[var_fn_sttn[1]].values[i]*10
    sttn_wnt.append('{:3d}_{:3d}'.format(
        line1.astype('int'), sttn1.astype('int')))

# clim years
yrs_clim = np.arange(yr_clim_bgn, yr_clim_end+1)
num_yrs_clim = len(yrs_clim)

# create time series matrix for IEA data and symbols
dataS_mtrx = np.zeros([num_sttn, num_var, num_quarter, num_yrs_clim])*np.nan
sttn_depth_vec = np.zeros(num_sttn)

for i in range(0, num_sttn):
    # --create input filename
    fn_in = dir_in + file_name_in + '_' + sttn_wnt[i] + '.' + file_type_in
    fn_in = '{}{}_{}_min_between_{}_{}m.nc'.format(
        dir_in, var_fn_wnt, sttn_wnt[i], z1, z2)

    # --open the netcdf files as an xarray
    ds_in = xr.open_dataset(fn_in)

    # station depth
    sttn_depth_vec[i] = ds_in.sttn_depth

    # year and month
    yrs_in = ds_in.time.dt.year.data
    mon_in = ds_in.time.dt.month.data

    # loop over var
    for j in range(num_var):
        da1 = ds_in[var_wnt[j]]
        for k in range(num_quarter):
            in_qrtr = np.where(ds_in[var_qrtr].data == k+1)[0]
            da1_qrtr = da1[in_qrtr]
            yrs_qrtr = yrs_in[in_qrtr]
            mon_qrtr = mon_in[in_qrtr]

            # find the years that in the clim
            in_clim = np.logical_and(
                yrs_qrtr >= yrs_clim[0], yrs_qrtr <= yrs_clim[-1]).nonzero()[0]

            da1_clim = da1_qrtr[in_clim]

            # station can be sampled multiple times, for example during 1997 and 1998
            da1_clim_year_mn = da1_clim.groupby('time.year').mean('time')
            yrs_sampled = da1_clim_year_mn.year.data

            # place in matrix
            ia = np.isin(yrs_clim, yrs_sampled)
            ib = np.isin(yrs_sampled, yrs_clim)
            dataS_mtrx[i, j, k, ia] = da1_clim_year_mn.data[ib]

# IEA data and symbols
dt_yrs = pd.to_datetime({'year': yrs_clim, 'month': np.ones(
    num_yrs_clim), 'day': np.ones(num_yrs_clim)})

for i in range(num_var):
    for j in range(num_quarter):
        # iea_list = list()
        ts_mtrx = np.zeros([num_sttn, num_yrs_clim])
        coord_mtrx = np.zeros([num_sttn, len(att_coord)])
        sttn_mtrx = np.zeros([num_sttn, len(att_sttn)])
        data_mtrx = np.zeros([num_sttn, len(data_wnt)])
        mrkr_mtrx = np.zeros([num_sttn, len(data_wnt)], dtype='str')

        for k in range(num_sttn):
            print('{}: {}, {}: {}, {}: {}'.format(i, num_var, j, num_quarter, k, num_sttn))
            data1 = np.squeeze(dataS_mtrx[k, i, j, :])
            ps1 = pd.Series(data=data1, index=dt_yrs).asfreq(freq_wnt)
            lon1 = df_sttn[var_fn_sttn[3]]+360
            lat1 = df_sttn[var_fn_sttn[2]]
            line1 = df_sttn[var_fn_sttn[0]]
            sttn1 = df_sttn[var_fn_sttn[1]]
            zd1 = C_iea(ps1, yr_clim_bgn=yr_clim_bgn,
                        yr_clim_end=yr_clim_end,
                        wndw=wndw, yy_end=yy_end,
                        lon=lon1[k], lat=lat1[k],
                        line=line1[k], sttn=sttn1[k])
            # place in matrices
            ts_mtrx[k, :] = data1

            for iii in range(len(att_coord)):
                coord_mtrx[k, iii] = getattr(zd1, att_coord[iii])

            for iii in range(len(att_sttn)):
                sttn_mtrx[k, iii] = getattr(zd1, att_sttn[iii])

            for iii in range(len(data_wnt)):
                data_mtrx[k, iii] = getattr(zd1, data_wnt[iii])

            for iii in range(len(mrkr_wnt)):
                mrkr_mtrx[k, iii] = getattr(zd1, mrkr_wnt[iii])

        # put into xr.da
        loc_indx = np.arange(0, num_sttn)
        da1 = xr.DataArray(coord_mtrx, coords=[loc_indx, att_coord],
                           dims=['index', 'coord'])
        da2 = xr.DataArray(data_mtrx, coords=[loc_indx, data_wnt],
                           dims=['index', 'data'])
        da3 = xr.DataArray(mrkr_mtrx, coords=[loc_indx, mrkr_wnt],
                           dims=['index', 'marker'])
        da4 = xr.DataArray(ts_mtrx, coords=[loc_indx, dt_yrs],
                           dims=['index', 'time'])
        da5 = xr.DataArray(sttn_mtrx, coords=[loc_indx, att_sttn],
                           dims=['index', 'station'])
        da6 = xr.DataArray(sttn_depth_vec, coords=[loc_indx], dims=['index'])

        # put into xr.ds
        ds1 = da1.to_dataset(name='coord_mtrx')
        ds1['data_mtrx'] = da2
        ds1['mrkr_mtrx'] = da3
        ds1['ts_mtrx'] = da4
        ds1['sttn_mtrx'] = da5
        ds1['sttn_depth'] = da6

        # create output name and save xarray dataset to netcdf file
        # --check if directory exist, if it doesn't then create
        try:
            os.makedirs(dir_out)
        except OSError:
            if not os.path.isdir(dir_out):
                raise

        fn_out = '{}data_mn5_trnd5_qrtr{}_{}_between_{}_{}m.nc'.format(
            dir_out, j+1, var_wnt[i], z1, z2)

        # --Save Dataset to a netcdf file
        ds1.to_netcdf(fn_out)

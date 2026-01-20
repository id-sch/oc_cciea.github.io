import os
import re
import shutil
import xarray as xr
import pandas as pd
import numpy as np
from numpy import ma

# --Create Netcdf files for the 113 "core" stations of CalCOFI. The data is
# --from the CalCOFI CSV and prelimary CTD data. These files can are
# --downloaded from the CalCOFI data page. This code should be run as often
# --as the data (CSV or prelimCTD) gets updated by CalCOFI.

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# --Notes:
#   1) I'm creating pandas dataframes (pd df) by reading CSV files of
#      cast, bottle, station, CTD downloaded from the CalCOFI website.
#   2) These files have these suffix: _c = cast, -b = bottle,
#      -s = station, _d = CTD
#   3) The pd df have column headers and the code in the "Change These"
#      section defines the columns wanted
#   4) Lists of columns wanted are made to be saved as variables in
#      a NetCDF files. The variables will be 2D (nt x nd) denoted by "var2"
#      or 1D (nt) denoted by "var1".
#   5) cc = consecative cast number, column in both the bottle and cast
#           CSV files
#   6) I now allow for you to set the output directory name

# output directory
dir_out = './data_gha/CalCofiCSV/'

# --max depth of data matrix (time x depth) for each station
max_dpth = 500

# --output file type extension
file_type = 'nc'

# input files
file_s = './csv_database_gha/calcofi_stations66.csv'
file_b = './csv_database_gha/194903-202105_Bottle.csv'
file_c = './csv_database_gha/194903-202105_Cast.csv'

# --Column header names for the 4 files
# ---------------------BOTTLE
cc_b = 'Cst_Cnt'
dpth_b = 'R_Depth'

# NOTE: May 16 2024, The R_Sal values do not make any sense to me, I think they are wrong 
# var2_b_wnt = ['R_TEMP', 'R_Sal', 'STheta', 'R_DYNHT',
#               'O2ml_L', 'O2Sat', 'SiO3uM', 'PO4uM', 'NO3uM', 'NO2uM', 'ChlorA', 'Phaeop']

# NOTE: May 16, 2024, I'm goingto use the 'other' temp and salinity columns and leave everything else the same
var2_b_wnt = ['T_degC', 'Salnty', 'STheta', 'R_DYNHT',
              'O2ml_L', 'O2Sat', 'SiO3uM', 'PO4uM', 'NO3uM', 'NO2uM', 'ChlorA', 'Phaeop']


# the R_TEMP is actually potential temperature, change the label to show this,
# also change the labels so that they can be used by old code
var2_b_lbl = ['R_POTEMP', 'R_SALINITY', 'R_SIGMA', 'R_DYNHT', 'R_O2', 'R_O2Sat', 'R_SIO3', 'R_PO4', 'R_NO3', 'R_NO2', 'R_CHLA', 'R_PHAEO']


# ---------------------CAST
cc_c = 'Cst_Cnt'
line_c = 'Rpt_Line'
sttn_c = 'Rpt_Sta'
lat_c = 'Lat_Dec'
lon_c = 'Lon_Dec'
qrtr_c = 'Quarter'
yr_c = 'Year'
jd_c = 'Julian_Day'
mon_c = 'Month'
var1_c_wnt = [line_c, sttn_c, lat_c, lon_c, qrtr_c, mon_c]

# ---------------------STATION
line_s = 'Line'
sttn_s = 'Station'

# --------------------- Prelim CTD
fn_d = './prelim_ctd.nc'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

# --get the size of input variables
num_var2_b = len(var2_b_wnt)
num_var1_c = len(var1_c_wnt)


# --interpolate variables over these depths
dpth_int = np.arange(0, max_dpth+1)
num_int = np.size(dpth_int)

# --get line/sttn of 113 "core" stations
df_s = pd.read_csv(file_s)
sttn_ps = df_s[line_s].astype(str)+'_'+df_s[sttn_s].astype(str)
sttn_wnt = sttn_ps.values
num_sttn = len(sttn_wnt)

# -- open the prelim CTD data
ds1_d = xr.open_dataset(fn_d)

# --------------------------------------------------------------------
# --Process Bottle data
# --------------------------------------------------------------------
# --Bottle CSV is very large, so use 'nrow' and 'use_cols' feature
#   of pd.read_csv to speed things up

# --column headers
df_hdr_b = pd.read_csv(file_b, nrows=1)
df_hdr_c = pd.read_csv(file_c, nrows=1)

# --just get the pandas dataframes of columns wanted
hdr_b = var2_b_wnt + [cc_b, dpth_b]
hdr_c = var1_c_wnt + [cc_c, yr_c, jd_c]
hdr_c.append(cc_b)
in_clmns_b = np.isin(df_hdr_b.columns.values, np.array(hdr_b)).nonzero()[0]
in_clmns_c = np.isin(df_hdr_c.columns.values, np.array(hdr_c)).nonzero()[0]

# --read cast CSV data, get pd Dataframe of hdr_c columns
df_c = pd.read_csv(file_c, usecols=in_clmns_c)

# --create time variable from df_c
time_c = pd.to_datetime(df_c[yr_c].astype(
    str) + df_c[jd_c].astype(str), format="%Y%j")

# --read Bottle CSV data, but only the consecative cast number
in_cc_b = (df_hdr_b.columns == cc_b).nonzero()[0]
df_cc_b = pd.read_csv(file_b, usecols=in_cc_b)

for i in range(num_sttn):
    print('{}: {}'.format(i, num_sttn))

    # ------------------------------------------------
    # --Subset Cast file by station
    in_sttn = ((df_c[line_c] == df_s[line_s][i]) &
               (df_c[sttn_c] == df_s[sttn_s][i]))
    # indx_sttn = in_sttn.nonzero()[0]
    indx_sttn = in_sttn.values.nonzero()[0]

    # --unique consecative cast number, ie number of CTD casts at the station
    sttn_time_c = time_c.values[indx_sttn]
    sttn_cc_c = np.unique(df_c[cc_c][in_sttn])
    nt_c = sttn_cc_c.size

    # --station var1 created from Cast csv
    sttn_var1_c = np.zeros([num_var1_c, nt_c])
    for ii in range(num_var1_c):
        sttn_var1_c[ii, :] = df_c[var1_c_wnt[ii]][in_sttn].values

    # --station var2 created from Bottle csv, loop over all nt
    #   and interpolate the profiles
    sttn_var2_b = np.zeros([num_var2_b, num_int, nt_c])*np.nan

    for j in range(nt_c):
        # station indices
        in_nt = (df_cc_b[cc_b] == sttn_cc_c[j])
        indx_nt = in_nt.values.nonzero()[0]

        # --create pd df of bottle data for this cc number
        df_b = pd.read_csv(file_b, usecols=in_clmns_b,
                           skiprows=indx_nt[0]+1, nrows=indx_nt.size, header=None)
        df_b.columns = df_hdr_b.columns.values[in_clmns_b]

        # --depths of profile, used to interpolate each variable
        dpth = df_b[dpth_b].values

        for k in range(0, num_var2_b):
            var_k = df_b[var2_b_wnt[k]].values
            # --masked array for missing values
            vard = ma.array(var_k, mask=np.isnan(var_k))
            dpthd = ma.array(dpth, mask=np.isnan(var_k))
            if dpthd.compressed().size > 0:
                data_int = np.interp(
                    dpth_int, dpthd.compressed(), vard.compressed(),
                    left=np.nan, right=np.nan)
            else:
                data_int = np.zeros(num_int)*np.nan
            sttn_var2_b[k, :, j] = data_int

    # --Extend the var1 and var2 of CTD,cast,bottle for
    #   final (_f) var1 an var2
    da_jd_d = ds1_d['jd_d'].sel(station=sttn_wnt[i])
    da_yr_d = ds1_d['yr_d'].sel(station=sttn_wnt[i])
    da_var1_d = ds1_d['var1_mtrx_d'].sel(station=sttn_wnt[i])
    da_var2_d = ds1_d['var2_mtrx_d'].sel(station=sttn_wnt[i])
    ind_d = np.isfinite(da_jd_d.values)
    jdd = da_jd_d.values[ind_d]
    yrd = da_yr_d.values[ind_d]
    nt_d = len(yrd)
    yr_jd_d = np.array([yrd[x].astype('int').astype('str')+jdd[x].astype('int').astype('str') for x in range(0, nt_d)])

    sttn_var1_d = da_var1_d.data[:, ind_d]
    sttn_var2_d = da_var2_d.data[:, :, ind_d]
    sttn_time_d = pd.to_datetime(yr_jd_d, format="%Y%j")

    time_f = np.concatenate((sttn_time_c, sttn_time_d))
    sttn_var1_f = np.concatenate((sttn_var1_c, sttn_var1_d), axis=1)
    sttn_var2_cmb = np.concatenate((sttn_var2_b, sttn_var2_d), axis=2)

    # --Remove depths of no data (or very little, need more than 1% to keep)
    ob_mtrx = np.copy(sttn_var2_cmb)
    in_ob = np.isfinite(ob_mtrx)
    ob_mtrx[in_ob] = 1
    ob_sum1 = np.nansum(ob_mtrx, 2)
    ob_sum2 = np.nansum(ob_sum1, 0)
    [nv, nd, nt] = ob_mtrx.shape
    max_ob = nv*nt
    in_keep = np.where(ob_sum2 > max_ob*0.01)[0]

    sttn_var2_f = sttn_var2_cmb[:, in_keep, :]
    dpth_f = dpth_int[in_keep]

    # -----------------
    # --Now add data_mtrx into a Xarray, can create Netcdf file from these
    # -----------------
    for k in range(0, num_var2_b):
        # put depth interpolated matrix (nd x nt) in a xr.DataArray
        da1 = xr.DataArray(sttn_var2_f[k, :, :], coords=[
                           dpth_f, time_f], dims=['depth', 'time'])
        # sort by time to make sure it is increasing
        da1s = da1.sortby(da1.time)

        # --add the DataArray to a Dataset, but first we need to create it,
        # --otherwise just add additional varialbes to the Dataset
        if k == 0:
            ds1 = da1s.to_dataset(name=var2_b_lbl[k])
        else:
            ds1[var2_b_lbl[k]] = da1s

    # --create DataArray from the var_vec_wnt
    for k in range(0, num_var1_c):
        # put vectors (nt) in a xr.DataArray
        va1 = xr.DataArray(sttn_var1_f[k, :], coords=[
                           time_f], dims=['time'])
        # sort by time to make sure it is increasing
        va1s = va1.sortby(va1.time)

        # --add this new DataArray to the Dataset
        ds1[var1_c_wnt[k]] = va1s

    # -------------------------------------------------------
    # --create output directory and save xarray dataset to netcdf file
    # -------------------------------------------------------
    pwd1 = os.getcwd()
    # dir_out = pwd1[0:12] + 'data_files' + pwd1[11:]
    fn_out = re.sub(r'[^\w]', '', sttn_wnt[i]) + '.' + file_type
    dir_fn = dir_out + '/' + fn_out

    # --check if directory exist, if it doesn't then create
    try:
        os.makedirs(dir_out)
    except OSError:
        if not os.path.isdir(dir_out):
            raise

    # --Save Dataset to a netcdf file
    ds1.to_netcdf(dir_fn)

# remove some large files that can not be commit to github
os.remove(fn_d)

# remove the directory and files that has the downloaded
# CalCOFI bottle and ctd csv files
shutil.rmtree('./csv_database_gha/')

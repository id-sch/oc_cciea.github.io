import xarray as xr
import pandas as pd
import numpy as np
import seawater as sw
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
#      station, CTD downloaded from the CalCOFI website, and a nutrient
#      csv provided by Rasmus to Steven (see work gmail on Sept 3, 2025)
#   2) These files have these suffix: 
#      -s = station, _d = CTD, _n = nutrient
#   3) The pd df have column headers and the code in the "Change These"
#      section defines the columns wanted
#   4) Lists of columns wanted are made to be saved as variables in
#      a NetCDF files. The variables will be 2D (nt x nd) denoted by "var2"
#      or 1D (nt) denoted by "var1".
#   5) cc = consecative cast number, column in both the bottle and cast
#           CSV files
#   6) I now allow for you to set the output directory name

# --max depth of data matrix (time x depth) for each station
max_dpth = 500

# --output file type extension
file_type = 'nc'
file_s = './csv_database_gha/calcofi_stations66.csv'

dir_ctd_d = './csv_database_gha/ctd_prelim/'
fn_d = ['20-2107SR_CTDBTL_001-071D.csv', '20-2111SR_CTDBTL_001-075D.csv',
        '20-2204SH_CTDBTL_001-101D.csv', '20-2208BH_CTDBTL_001-067D.csv', '20-2211SR_CTDBTL_001-073D.csv',
        '20-2301RL_CTDBTL_001-072D.csv', '20-2304SH_CTDBTL_001-114D.csv', '20-2307SR_CTDBTL_001-056D_057-069D_combined_IDS.csv', '20-2311SR_CTDBTL_001-050D_051-075D_combined_IDS.csv',
        '20-2401RL_CTDBTL_001-016D_017-033D_IDS.csv',  '20-2404SH_CTDBTL_001-116D.csv', '20-2408SR_CTDBTL_001-075D.csv', '20-2411SR_CTDBTL_001-048D.csv',
        '20-2501RL_CTDBTL_001-116D.csv', '20-2504SH_CTDBTL_001-116D_IDS.csv', '20-2507_CTD_001-055D.csv', '20-2511_CTD_001-047D.csv']
yr_d = [2021, 2021,
        2022, 2022, 2022,
        2023, 2023, 2023, 2023,
        2024, 2024, 2024, 2024,
        2025, 2025, 2025, 2025]
qrtr_d = [3, 4,
          2, 3, 4,
          1, 2, 3, 4,
          1, 2, 3, 4,
          1, 2, 3, 4]

# set the dataframe variables for the CTD, some CSV files have missing sensors and the 'Ave' columns are not correct, this will set the index of a list of var2_d_wnt
in_var2_d_list = [4, 0,
                  4, 0, 0,
                  0, 4, 0, 2,
                  1, 0, 0, 0,
                  0, 4, 3, 3]
in_var2_d_list = [5, 0,
                  4, 0, 0,
                  5, 4, 0, 2,
                  1, 0, 0, 0,
                  0, 4, 3, 3]

# ---------------------STATION
# line_s = 'Line '
line_s = 'Line'
sttn_s = 'Station'

# ---------------------CTD
line_d = 'Line'
sttn_d = 'Sta'
dpth_d = 'Depth'
cc_d = 'Ord_Occ'
time_d = 'Date_Time_PST'
lon_d = 'Lon_Dec'
lat_d = 'Lat_Dec'
var1_d_wnt = [line_d, sttn_d, lat_d, lon_d]

# can have different columns based on the CTD CSV and the available columns,
# another one can be created if these don't account for columns with data
# NOTE: Rasmus cautioned using the O2 station avg data, suggested using Ox1 or Ox1_StaCorr
var2_d_wnt1 = ['TempAve', 'Salt2_Corr', 'SigThetaTS1', 'DynHt',
               'Ox1_StaCorr', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']

var2_d_wnt2 = ['Temp1', 'Salt1_Corr', 'SigThetaTS1', 'DynHt',
               'Ox1_StaCorr', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']
var2_d_wnt3 = ['Temp2', 'Salt2_Corr', 'SigThetaTS2', 'DynHt',
               'Ox1_StaCorr', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']

var2_d_wnt4 = ['TempAve', 'Salt2', 'SigThetaTS1', 'DynHt',
               'Ox1', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']

var2_d_wnt5 = ['TempAve', 'Salt2_Corr', 'SigThetaTS1', 'DynHt',
               'Ox1', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']

var2_d_wnt6 = ['TempAve', 'Salt2_Corr', 'SigThetaTS1', 'DynHt',
               'Ox1', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']


var2_d_list = [var2_d_wnt1, var2_d_wnt2, var2_d_wnt3, var2_d_wnt4, var2_d_wnt5, var2_d_wnt6]

# ----------------------Nutrient CSV
file_n = './csv_database_gha/nutrient_prelim/compiled_finalized_nutrients_2111_2504_ALE.csv'

# ------------------Nutrient
line_n = 'Line'
sttn_n = 'Sta'
dpth_n = 'Depth'
cc_n = 'Ord_Occ'
time_n = 'Date_Time_PST'
lon_n = 'Lon_Dec'
lat_n = 'Lat_Dec'
var1_n_wnt = [line_n, sttn_n, lat_n, lon_n]

var2_n_wnt = ['NO3', 'NO2', 'PO4', 'SIL']

# -- output NetCDF file
# set the variables with Bottle CSV labels, will make it easier to merge
# the prelim CTD data with the bottle data
var2_b_lbl = ['R_POTEMP', 'R_SALINITY', 'R_SIGMA', 'R_DYNHT', 'R_O2', 'R_O2Sat', 'R_SIO3', 'R_PO4', 'R_NO3', 'R_NO2', 'R_CHLA', 'R_PHAEO']
var1_d_lbl = ['line', 'sttn', 'lat', 'lon', 'quarter', 'month']

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

# --get the size of input variables
num_var1_d = len(var1_d_wnt)
num_var2_d = len(var2_d_list[0])
num_qrtr_d = len(qrtr_d)
num_var2_n = len(var2_n_wnt)

# --interpolate variables over these depths
dpth_int = np.arange(0, max_dpth+1)
num_int = np.size(dpth_int)

# --get line/sttn of 66 "core" stations
df_s = pd.read_csv(file_s)
sttn_ps = df_s[line_s].astype(str)+'_'+df_s[sttn_s].astype(str)
sttn_wnt = sttn_ps.values
num_sttn = len(sttn_wnt)

# --------------------------------------------------------------------
# --Process CTD data
# --------------------------------------------------------------------

# --interpolate variables over these depths
dpth_int = np.arange(0, max_dpth+1)
num_int = np.size(dpth_int)

# --get line/sttn of 113 "core" stations
df_s = pd.read_csv(file_s)
sttn_ps = df_s[line_s].astype(str)+'_'+df_s[sttn_s].astype(str)
sttn_wnt = sttn_ps.values
num_sttn = len(sttn_wnt)

# --------------------------------------------------------------------
# --Process CTD data
# --------------------------------------------------------------------

# -- create CTD data vector (sttns X variables_1d X quarters)
#    note: +2 to num_var1_d is for the quarter and month variable have to input
#          these by hand since there is no column data with this info
var1_mtrx_d = np.zeros([num_sttn, len(var1_d_lbl), num_qrtr_d])*np.nan

# --create ctd data matrix (sttns X variables_2d X depth X quarters)
var2_mtrx_d = np.zeros([num_sttn, num_var2_d, num_int, num_qrtr_d])*np.nan

# -- create Year and Julian date vectors to create a date, do this to
#    match year and jd columns in the pd df of the Cast file
yr_vec_d = np.zeros([num_sttn, num_qrtr_d])*np.nan
jd_vec_d = np.zeros([num_sttn, num_qrtr_d])*np.nan

for ii in range(0, num_qrtr_d):
    # get the column variable names
    var2_d_wnt = var2_d_list[in_var2_d_list[ii]]

    # --create ctd filename and open as pd Dataframe
    d_file = '{}/{}/q{}/db-csvs/{}'.format(dir_ctd_d, yr_d[ii], qrtr_d[ii], fn_d[ii])
    df_d = pd.read_csv(d_file, low_memory=False)

    # --change the "generic" Dataframe index to an datetime index
    # df_d['Datetime'] = pd.to_datetime(df_d[time_d].astype(str), infer_datetime_format=True)
    df_d['Datetime'] = pd.to_datetime(df_d[time_d].astype(str))
    df_d = df_d.set_index('Datetime')

    # --get year and jd from Datetime index, need these if more than one
    # --cast per station, will take average of the dates for final date
    jd = df_d.index.dayofyear
    year = df_d.index.year
    month = df_d.index.month

    # --
    for jj in range(0, num_sttn):
        # --indices of all casts for a give station
        in_d_wnt = ((df_d[line_d] == df_s[line_s][jj]) &
                    (df_d[sttn_d] == df_s[sttn_s][jj]))
        # indx_d_wnt = in_d_wnt.nonzero()[0]
        indx_d_wnt = in_d_wnt.values.nonzero()[0]

        # --unique consecative cast number, ie number of casts at the station
        cc_d_unq = np.unique(df_d[cc_d][indx_d_wnt])
        num_d_unq = cc_d_unq.size

        if num_d_unq > 0:
            # --create data matrix (variables X interp depths X unique dates)
            unq_mtrx_d = np.zeros([num_var2_d, num_int, num_d_unq])*np.nan
            yr_unq = np.zeros([num_d_unq])
            mon_unq = np.zeros([num_d_unq])
            jd_unq = np.zeros([num_d_unq])
            lon_unq = np.zeros([num_d_unq])
            lat_unq = np.zeros([num_d_unq])

            for kk in range(0, num_d_unq):
                # --indices of an individual cast
                in_dq = (df_d[cc_d][in_d_wnt] == cc_d_unq[kk])
                indx_dq = in_dq.values.nonzero()[0]

                # --depths of profile, used to interpolate each variable
                dpth_dq = df_d[dpth_d].iloc[indx_d_wnt[indx_dq]].values
                yr_unq[kk] = year[indx_d_wnt[in_dq]].values[0]
                mon_unq[kk] = month[indx_d_wnt[in_dq]].values[0]
                jd_unq[kk] = jd[indx_d_wnt[in_dq]].values[0]
                lon_unq[kk] = df_d[lon_d].iloc[indx_d_wnt[in_dq]].values[0]
                lat_unq[kk] = df_d[lat_d].iloc[indx_d_wnt[in_dq]].values[0]

                for ll in range(0, num_var2_d):
                    lbl_var_wnt = var2_d_wnt[ll]

                    # Special exceptions to the oxygen data, mostly for "unimportant"
                    # stations. This is very ugly code, but it was constructed
                    # by inspecting each ctd data file using Libre office and identifying
                    # oxygen columns that would work and inspecting plots of all 66
                    # stations over all depths to find stations that had extreme DO values.
                    if fn_d[ii] == '20-2107SR_CTDBTL_001-071D.csv' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 28.0 and ll == 4:
                        lbl_var_wnt = 'Ox2_StaCorr'
                    if fn_d[ii] == '20-2107SR_CTDBTL_001-071D.csv' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 30.0 and ll == 4:
                        lbl_var_wnt = 'Ox2_StaCorr'
                    if fn_d[ii] == '20-2208BH_CTDBTL_001-067D.csv' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 55.0 and ll == 4:
                        lbl_var_wnt = 'Ox2_StaCorr'
                    if fn_d[ii] == '20-2301RL_CTDBTL_001-072D.csv' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 28.0 and ll == 4:
                        lbl_var_wnt = 'Ox2'
                    if fn_d[ii] == '20-2307SR_CTDBTL_001-056D_057-069D_combined_IDS.csv' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 28.0 and ll == 4:
                        lbl_var_wnt = 'Ox2'
                    if fn_d[ii] == '20-2307SR_CTDBTL_001-056D_057-069D_combined_IDS.csv' and df_s['Line'].values[jj] == 90.0 and df_s['Station'].values[jj] == 110.0 and ll == 4:
                        lbl_var_wnt = 'Ox2'
                    if fn_d[ii] == '20-2311SR_CTDBTL_001-050D_051-075D_combined_IDS.csv' and df_s['Line'].values[jj] == 83.3 and df_s['Station'].values[jj] == 55.0 and ll == 4:
                        lbl_var_wnt = 'Ox2'
                    if fn_d[ii] == '20-2311SR_CTDBTL_001-050D_051-075D_combined_IDS.csv' and df_s['Line'].values[jj] == 90.0 and df_s['Station'].values[jj] == 120.0 and ll == 4:
                        lbl_var_wnt = 'Ox2'
                    if fn_d[ii] == '20-2404SH_CTDBTL_001-116D.csv' and df_s['Line'].values[jj] == 90.0 and df_s['Station'].values[jj] == 53.0 and ll == 4:
                        lbl_var_wnt = 'Ox1'
                    if fn_d[ii] == '20-2408SR_CTDBTL_001-075D.csv' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 30.0 and ll == 4:
                        lbl_var_wnt = 'Ox1'

                    # interpolate
                    var_ll = df_d[lbl_var_wnt].iloc[indx_d_wnt[indx_dq]].values
                    # --masked array for missing values
                    vard = ma.array(var_ll, mask=np.isnan(var_ll))
                    dpthd = ma.array(dpth_dq, mask=np.isnan(var_ll))
                    if dpthd.compressed().size > 0:
                        data_int = np.interp(dpth_int, dpthd.compressed(
                        ), vard.compressed(), left=np.nan, right=np.nan)
                    else:
                        data_int = np.zeros(num_int)*np.nan

                    unq_mtrx_d[ll, :, kk] = data_int

                # recalculate potential temperature and replace the one that
                # was in the CSV file, want to do this because sometimes the
                # CSV file doesn't have data in the PoT1 column
                t1 = np.squeeze(unq_mtrx_d[0, :, kk])
                s1 = np.squeeze(unq_mtrx_d[1, :, kk])
                datad = sw.ptmp(s1, t1, dpth_int, 0)
                unq_mtrx_d[0, :, kk] = datad


            # --2D variables, take mean of all casts to get a single profile
            # --for a given station
            var2_mtrx_d[jj, :, :, ii] = np.nanmean(unq_mtrx_d, axis=2)

            # --1D variables, take mean of yr and jd to get single date
            yr_vec_d[jj, ii] = np.mean(yr_unq)
            jd_vec_d[jj, ii] = np.mean(jd_unq)
            var1_mtrx_d[jj, 0, ii] = df_s[line_s][jj]
            var1_mtrx_d[jj, 1, ii] = df_s[sttn_s][jj]
            var1_mtrx_d[jj, 2, ii] = np.mean(lat_unq)
            var1_mtrx_d[jj, 3, ii] = np.mean(lon_unq)
            var1_mtrx_d[jj, 4, ii] = qrtr_d[ii]
            var1_mtrx_d[jj, 5, ii] = np.mean(mon_unq)
# --------------------------------------------------------------------
# --Process Nutrient data
# --------------------------------------------------------------------
# Use the nutrient CSV to add the nutrients (most of the CTD csv files
# do not have nutrient data, but some do -- such as Summer 2021).
df_n = pd.read_csv(file_n)
study_all_n = np.unique(df_n['Study'].values)

study_d = []
for i in range(len(fn_d)):
    study_type1 = fn_d[i].split('-')[1]
    study1 = study_type1.split('_')[0]
    type1 = study_type1.split('_')[1]

    if type1 == 'CTDBTL':
        study_d.append(study1)

# remove the extra Feb 2025 survey
in_study_wnt = np.where(np.array(study_all_n) != '2502RL')[0]
study_n = study_all_n[in_study_wnt]

num_qrtr_n = len(study_n)

for ii in range(num_qrtr_n):
    df_n1 = df_n.loc[df_n['Study']==study_n[ii]]

    # index that the nutrient data matches the quarter of the ctd_prelim
    in_qrtr_mtrx = np.where(np.array(study_d)==study_n[ii])[0]

    # --change the "generic" Dataframe index to an datetime index
    df_n1['Datetime'] = pd.to_datetime(df_n1[time_d].astype(str))
    df_n1 = df_n1.set_index('Datetime')

    # --get year and jd from Datetime index, need these if more than one
    # --cast per station, will take average of the dates for final date
    jd = df_n1.index.dayofyear
    year = df_n1.index.year
    month = df_n1.index.month

    # --
    for jj in range(0, num_sttn):
        # --indices of all casts for a give station
        in_n_wnt = ((df_n1[line_d] == df_s[line_s][jj]) &
                    (df_n1[sttn_d] == df_s[sttn_s][jj]))
        indx_n_wnt = in_n_wnt.values.nonzero()[0]

        # --unique consecative cast number, ie number of casts at the station
        cc_n_unq = np.unique(df_n1[cc_n].iloc[indx_n_wnt])
        num_n_unq = cc_n_unq.size

        if num_n_unq > 0:
            # --create data matrix (variables X interp depths X unique dates)
            unq_mtrx_n = np.zeros([num_var2_n, num_int, num_n_unq])*np.nan

            for kk in range(0, num_n_unq):
                # --indices of an individual cast
                in_nq = (df_n1[cc_n][in_n_wnt] == cc_n_unq[kk])
                indx_nq = in_nq.values.nonzero()[0]

                # --depths of profile, used to interpolate each variable
                dpth_nq = df_n1[dpth_n].iloc[indx_n_wnt[indx_nq]].values

                for ll in range(0, num_var2_n):
                    var_ll = df_n1[var2_n_wnt[ll]].iloc[indx_n_wnt[indx_nq]].values
                    # --masked array for missing values
                    vard = ma.array(var_ll, mask=np.isnan(var_ll))
                    dpthn = ma.array(dpth_nq, mask=np.isnan(var_ll))
                    if dpthn.compressed().size > 0:
                        data_int = np.interp(
                            dpth_int, dpthn.compressed(), vard.compressed(),
                            left=np.nan, right=np.nan)
                    else:
                        data_int = np.zeros(num_int)*np.nan

                    unq_mtrx_n[ll, :, kk] = data_int

            unq_mtrx_n_mn = np.nanmean(unq_mtrx_n, axis=2)
            for ll in range(num_var2_n):
                in_var_mtrx = np.where(np.array(var2_d_wnt) == var2_n_wnt[ll])[0]
                var2_mtrx_d[jj, in_var_mtrx, :, in_qrtr_mtrx] = unq_mtrx_n_mn[ll, :]

# Create a netcdf file with the prelim data
indx_time_d = np.arange(0, num_qrtr_d)
da1_out = xr.DataArray(var1_mtrx_d, coords=[sttn_wnt, var1_d_lbl, indx_time_d], dims=['station', 'var1', 'index_time'])
da2_out = xr.DataArray(var2_mtrx_d, coords=[sttn_wnt, var2_b_lbl, dpth_int, indx_time_d], dims=['station', 'var2', 'depth', 'index_time'])
da3_out = xr.DataArray(yr_vec_d, coords=[sttn_wnt, indx_time_d], dims=['station', 'index_time'])
da4_out = xr.DataArray(jd_vec_d, coords=[sttn_wnt, indx_time_d], dims=['station', 'index_time'])

ds1_out = da1_out.to_dataset(name='var1_mtrx_d')
ds1_out['var2_mtrx_d'] = da2_out
ds1_out['yr_d'] = da3_out
ds1_out['jd_d'] = da4_out

fn_out = 'prelim_ctd.nc'
ds1_out.to_netcdf(fn_out)

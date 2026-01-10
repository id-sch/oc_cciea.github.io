import datetime as DT
import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
from fun_pd_df2csvR_time import fun_pd_df2csvR_time


# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------
# --directory of the 6hr UI
dir_out = './data_gha/bakunUI/'
# dir_out = './data_x13/bakunUI/'
dir_in = dir_out
fn_in = '{}UI_daily.nc'.format(dir_in)

# -- Input variables, change these
lat_wnt = [33, 39, 45]

# rolling mean
roll = 1

# output directory
dir_out = './csv_for_erddap/'

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# len input variables
num_wnt = len(lat_wnt)

# --open the netcdf files as an xarray
ds1 = xr.open_dataset(fn_in)

var_ds1 = list(ds1.keys())
coord_ds1 = list(ds1.coords.keys())

# construct time array to check if any missing dates
dtD = np.arange(ds1.time.data[0].astype('datetime64[D]'),
                ds1.time.data[-1].astype('datetime64[D]')+1,
                dtype='datetime64[D]')

# years
yrs = np.unique(ds1.time.dt.year.data)
num_yrs = len(yrs)

# --loop over all lats wanted, open file and create sti,lusi,tumi
ftiA = np.zeros([num_yrs, num_wnt])
stiA = np.zeros([num_yrs, num_wnt])
lusiA = np.zeros([num_yrs, num_wnt])
tumiA = np.zeros([num_yrs, num_wnt])

for i in range(0, num_wnt):
    # get lat
    da1_lat = ds1[var_ds1[0]].sel(lat=lat_wnt[i])
    ui1D = da1_lat.data
    tt1D = da1_lat.time.data.astype('datetime64[D]')

    # -- ui1D can have missing values, interpolate to remove NaNs
    ui1Dd = ma.array(ui1D, mask=np.isnan(ui1D))  # Use a mask to mark nan
    ttd = ma.array(tt1D, mask=np.isnan(ui1D))  # Use a mask to mark nan
    z1 = ttd.compressed()
    z2 = ui1Dd.compressed()

    ui1Ddi = np.interp(pd.to_datetime(dtD).to_julian_date(),
                       pd.to_datetime(z1).to_julian_date(),
                       z2, left=np.nan, right=np.nan)

    # dataarray for interpolated data
    da1_intrp = xr.DataArray(ui1Ddi, coords=[dtD.astype('datetime64')], dims=['time'])

    # rolling mean
    da1_roll = da1_intrp.rolling(time=roll, center=True).mean()

    # rolling mean removes some of the data at the beggining of the time series, replace
    # with the mean of the missing data
    # ind_bgn = np.isnan(da1_roll.data[0:roll])
    # data_bgn_mn = np.nanmean(da1_intrp.data[0:roll][ind_bgn])
    data_bgn_mn = np.nanmean(da1_intrp.data[0:roll])
    da1_roll[0:roll] = data_bgn_mn

    # --create sti,lusi,tumi
    for j in range(0, num_yrs):
        da_yr = da1_roll.sel(time=da1_intrp.time.dt.year.isin(yrs[j]))
        dt_yr = da_yr.time.data
        ui_yr = da_yr.data


        # check for nan, can occur for last year
        ind = np.isfinite(ui_yr)
        ui_yr = ui_yr[ind]
        dt_yr = dt_yr[ind]

        # --check size of dt_yr, can be 365, 366 or
        # --less (depending on mon_wnt1,mon_wnt2)
        num_in = len(dt_yr)
        in_end = 365
 
        if num_in < in_end:
            in_end = num_in-1
            

        print('{}: {}, {}: {}, num_in={}, in_end={}'.format(i, num_wnt, j, num_yrs, num_in, in_end))
        
        ui_365 = np.zeros(365)*np.nan
        ui_365[0:in_end] = ui_yr[0:in_end]
        cui = np.nancumsum(ui_365)

        # --find the index of the minimum cui value in the spring,
        # --this will be the day of the Spring Transition (ST)
        in_sprng = np.where(dt_yr <= np.datetime64(DT.date(yrs[j], 6, 30)))
        cui_sprng = cui[in_sprng]
        in_sti = np.argmin(cui_sprng)

        # --find the maximum cui value after the ST, this will
        # --give you LUSI, the length of the upwelling season
        cui_from_sti = cui[in_sti:]
        in_lusi = np.nanargmax(cui_from_sti)

        # --the Fall Transition (FT) will be the date from the
        # --spring transition plus lusi
        in_fti = in_sti + in_lusi

        # --calculate TUMI, ie the cumsum of ui between ST and FT
        ui_for_tumi = ui_365[in_sti:in_fti]
        tumi = np.sum(ui_for_tumi)

        # --place in vectors, need to add one to convert the index to an actual day value
        # 1) first check to see if year has data upto june 30
        jd_end = pd.to_datetime(dt_yr[-1]).to_julian_date()
        jd_end_sprng = pd.to_datetime(DT.date(yrs[j], 6, 30)).to_julian_date()
        if jd_end < jd_end_sprng:
            sti = np.nan
        else:
            sti = in_sti + 1

        # 2) check to see if there are 365 or 366 dates in data
        nt = len(dt_yr)
        if nt < 365:
            fti = np.nan
            lusi = np.nan
            tumi = np.nan
        else:
            fti = in_fti + 1
            lusi = in_lusi + 1

        ftiA[j, i] = fti
        stiA[j, i] = sti
        lusiA[j, i] = lusi
        tumiA[j, i] = tumi

# --------------------------------------
# --Write R-style (IEA format) CSV files
# --------------------------------------

lon = np.zeros(num_wnt)
depth = np.zeros(num_wnt)

# --create ndarrays of yrs, data and std
yrsA = np.dot(np.transpose(np.ones([1, num_yrs])*yrs), np.ones([1, num_wnt]))
stdA = np.zeros([num_yrs, num_wnt])

# R style csv will have these columns
# clmns_iea = ['year', 'time', 'data', 'error', 'SElo', 'SEup',
#              'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']
clmns_iea = ['year', 'time', 'index', 'error', 'SElo', 'SEup',
             'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

# --STI
file_pre = 'oc_sti'
metric_lbl = 'STI (yearday)'
ts_lbl = []
df_list = list()
lat_list = list()
lon_list = list()
depth_list = list()
for i in range(num_wnt):
    # time series label
    dd1 = 'Spring Transition Index @ {}N'.format(lat_wnt[i])
    ts_lbl.append(dd1)

    # last year can be NaN, remove if it is
    if np.isnan(stiA[-1, i]):
        y1 = yrsA[0:-1, i]
        d1 = stiA[0:-1, i]
        sd1 = stdA[0:-1, i]
    else:
        y1 = yrsA[:, i]
        d1 = stiA[:, i]
        sd1 = stdA[:, i]
    # dataframe
    date1 = pd.to_datetime({'year': y1, 'month': np.ones(len(y1)), 'day': np.ones(len(y1))})
    df1 = pd.DataFrame({'Datetime': date1, 'data': d1, 'sd': sd1})
    df1 = df1.set_index('Datetime')
    df_list.append(df1)

    # yrs clim
    yrs_csv_bgn = int(df1.index.year[0])
    yrs_csv_end = int(df1.index.year[-1])

    # lat, lon, depth lists
    lat_list.append(lat_wnt[i])
    lon_list.append(lon[i])
    depth_list.append(depth[i])

fn_out = '{}_A.csv'.format(file_pre)
fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_list, lat_list, lon_list, depth_list,
                            metric_lbl, ts_lbl, dir_out, fn_out,
                            yrs_csv_bgn, yrs_csv_end)

# --LUSI
file_pre = 'oc_lusi'
metric_lbl = 'LUSI (days)'
ts_lbl = []
df_list = list()
lat_list = list()
lon_list = list()
depth_list = list()
for i in range(num_wnt):
    # time series label
    dd1 = 'Length of Upwelling Season Index @ {}N'.format(lat_wnt[i])
    ts_lbl.append(dd1)

    # last year can be NaN, remove if it is
    if np.isnan(lusiA[-1, i]):
        y1 = yrsA[0:-1, i]
        d1 = lusiA[0:-1, i]
        sd1 = stdA[0:-1, i]
    else:
        y1 = yrsA[:, i]
        d1 = lusiA[:, i]
        sd1 = stdA[:, i]
    # dataframe
    date1 = pd.to_datetime({'year': y1, 'month': np.ones(len(y1)), 'day': np.ones(len(y1))})
    df1 = pd.DataFrame({'Datetime': date1, 'data': d1, 'sd': sd1})
    df1 = df1.set_index('Datetime')
    df_list.append(df1)

    # yrs clim
    yrs_csv_bgn = int(df1.index.year[0])
    yrs_csv_end = int(df1.index.year[-1])

    # lat, lon, depth lists
    lat_list.append(lat_wnt[i])
    lon_list.append(lon[i])
    depth_list.append(depth[i])

fn_out = '{}_A.csv'.format(file_pre)
fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_list, lat_list, lon_list, depth_list,
                            metric_lbl, ts_lbl, dir_out, fn_out,
                            yrs_csv_bgn, yrs_csv_end)


# --TUMI
file_pre = 'oc_tumi'
metric_lbl = 'TUMI (m$^3$/s/100m coastline)'
ts_lbl = []
df_list = list()
lat_list = list()
lon_list = list()
depth_list = list()
for i in range(num_wnt):
    # time series label
    dd1 = 'Total Upwelling Magnitude Index @ {}N'.format(lat_wnt[i])
    ts_lbl.append(dd1)

    # last year can be NaN, remove if it is
    if np.isnan(tumiA[-1, i]):
        y1 = yrsA[0:-1, i]
        d1 = tumiA[0:-1, i]
        sd1 = stdA[0:-1, i]
    else:
        y1 = yrsA[:, i]
        d1 = tumiA[:, i]
        sd1 = stdA[:, i]
    # dataframe
    date1 = pd.to_datetime({'year': y1, 'month': np.ones(len(y1)), 'day': np.ones(len(y1))})
    df1 = pd.DataFrame({'Datetime': date1, 'data': d1, 'sd': sd1})
    df1 = df1.set_index('Datetime')
    df_list.append(df1)

    # yrs clim
    yrs_csv_bgn = int(df1.index.year[0])
    yrs_csv_end = int(df1.index.year[-1])

    # lat, lon, depth lists
    lat_list.append(lat_wnt[i])
    lon_list.append(lon[i])
    depth_list.append(depth[i])

fn_out = '{}_A.csv'.format(file_pre)
fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_list, lat_list, lon_list, depth_list,
                            metric_lbl, ts_lbl, dir_out, fn_out,
                            yrs_csv_bgn, yrs_csv_end)

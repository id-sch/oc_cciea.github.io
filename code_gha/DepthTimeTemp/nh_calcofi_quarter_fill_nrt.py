import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
from matplotlib import rcParams


def sttn_create_lbl(ds, sttn_pre, sttn_lst1):
    ttl = sttn_pre
    for i in range(0, len(sttn_lst1)):
        ttl = ttl + str(ds[sttn_lst1[i]].data[0]) + ' '
    return ttl

# plot paramaters
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'RdYlBu_r',
    'axes.grid': False,
    'savefig.dpi': 300,  # to adjust notebook inline plot size
    'xtick.top':        False,  # shold the top and bottom have tick marks
    'xtick.bottom':     True,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'ytick.direction': 'out',
    'xtick.direction': 'out',
    'axes.labelsize': 9,  # fontsize for x and y labels 
    'axes.titlesize': 9,
    'font.size': 9,  # was 10
    'legend.fontsize': 9,  # was 10
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.figsize': [11, 8.5],
    'font.family': 'STIXGeneral',
    'mathtext.fontset': 'stix',
    'toolbar': 'None',
    'savefig.bbox': 'tight'
}
rcParams.update(params)

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2025

# --station wanted, two stations newport hydrographic at 5 and 25 km
sttn_wnt = ['NH25_CTD', '933_300']
# sttn_wnt = ['NH05_CTD', '933_267']
num_sttn = len(sttn_wnt)

# --input directory
dir_in = ['./',
          './data_gha/CalCofiCSV/']

# --input file type extension
file_type_in = 'nc'

# --names of salinity and temperature in the netcdf file
var_wnt = ['Temperature', 'R_POTEMP']
var_name = ['Temperature', 'Temperature']
sttn_pre = ['Newport Hydrographic Station ', 'CalCOFI Station ']

# --station can have multiple names, create list of station names
sttn_lst = [['Station'], ['Rpt_Line', 'Rpt_Sta']]

# change the months of CalCOFI data to a uniform 'seasonal' spacing
monS = [1, 4, 7, 10]


# clim quarter months
mon_qrtr = [1, 4, 7, 10]

ylbl = 'Depth (m)'
# can have seperate label for each ds
clrbr_lbl = [var_name[0]+u' Anom. (\N{DEGREE SIGN}C)']
yr_clim_bgn = 1997
yr_clim_end = iea_yr
dpth_bgn = 0
dpth_end = 250

# min/max contour limit
nlvl = 3

# Mark missing samples on top of the x-axis
mark_missing_NH = 1
mark_missing_CC = 1

# directory out
dir_plot_out = './figures_gha/DepthTimeTemp/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# define fireice colormap
clr_map = mpl.colors.LinearSegmentedColormap.from_list(
    "", ["mediumblue", "dodgerblue", "cyan", "whitesmoke", "yellow", "tomato", "red"])

# months over clim time period
time_bgn = np.datetime64('{}-01'.format(yr_clim_bgn), 'M')
time_end = np.datetime64('{}-12'.format(yr_clim_end), 'M')
dtM = np.arange(time_bgn, time_end+1, dtype='datetime64[M]')
ntM = len(dtM)

# ----------------------------------------------------------------------------- Newport DATA
fn_in = dir_in[0] + '/' + sttn_wnt[0] + '.' + file_type_in
ds0 = xr.open_dataset(fn_in)
ds0M = ds0.resample(time='1MS').mean()

# nt, nd = ds0M[var_wnt[0]].shape
nt = ds0M.time.shape[0]
nd = ds0M.depth.shape[0]


# find missing dates
ia = np.in1d(dtM, ds0M.time.data).nonzero()[0]
ib = np.in1d(ds0M.time.data, dtM).nonzero()[0]

# create new data array on same time scale as dtM, with missing months marked with np.nan
data_new = np.zeros([nd, ntM])*np.nan
# data_old = ds0M[var_wnt[0]].data
data_old = ds0M[var_wnt[0]].data.T
data_new[:, ia] = data_old[:, ib]

da1NH = xr.DataArray(data_new, coords=[ds0M.depth.data, dtM.astype('datetime64')], dims=['depth', 'time'])


da0 = da1NH

# create anoms
clm0 = da0.groupby('time.month').mean('time')
da0A = da0.groupby('time.month') - clm0
sl = sttn_lst[0]
ttl0 = sttn_create_lbl(ds0, sttn_pre[0], sl)


# ----------------------------------------------------------------------------- CalCOFI data
fn_in = dir_in[1] + '/' + sttn_wnt[1] + '.' + file_type_in
ds1CC = xr.open_dataset(fn_in)

# get data on clim time period
da1CC = ds1CC[var_wnt[1]].sel(
    time=ds1CC[var_wnt[1]].time.dt.year.isin(np.arange(yr_clim_bgn, yr_clim_end+1)))

# summer 2020 has very large value > 40, remove and set to nan
time_correct = np.array('2020-07-13T00:00:00.000000000', dtype='datetime64[ns]')

in_time_correct = np.where(da1CC.time.data == time_correct)[0]
in_z_correct = np.logical_and(da1CC.depth.data >= 42, da1CC.depth.data <= 51).nonzero()[0]
data_correct = da1CC[:, in_time_correct].data.squeeze()
data_correct[in_z_correct] = np.nan
ind_correct = np.isfinite(data_correct)
data_int = np.interp(da1CC.depth.data, da1CC.depth.data[ind_correct],
                     data_correct[ind_correct], left=np.nan, right=np.nan)
da1CC[:, in_time_correct[0]] = data_int

#
da1CC_qrtr = ds1CC['Quarter'].sel(
    time=ds1CC['Quarter'].time.dt.year.isin(np.arange(yr_clim_bgn, yr_clim_end+1)))

qrtr = np.unique(da1CC_qrtr.data)
num_qrtr = len(qrtr)

yrs = np.unique(da1CC.time.dt.year)
num_yrs = len(yrs)

depthCC = da1CC.depth.data
nd = len(depthCC)

dataCC_qrtr = np.zeros([nd, num_yrs, num_qrtr])*np.nan

for i in range(num_yrs):
    da1CC_yr = da1CC.sel(time=da1CC.time.dt.year.isin(yrs[i]))
    da1CC_qrtr_yr = da1CC_qrtr.sel(time=da1CC_qrtr.time.dt.year.isin(yrs[i]))
    for j in range(num_qrtr):
        in_qrtr = np.where(da1CC_qrtr_yr == qrtr[j])[0]
        data_qrtr = da1CC_yr.data[:, in_qrtr]
        data_qrtr_mn = np.nanmean(data_qrtr, axis=1)
        dataCC_qrtr[:, i, j] = data_qrtr_mn

#
anomCC_qrtr = np.zeros([nd, num_yrs, num_qrtr])
for i in range(num_qrtr):
    clmn1 = np.nanmean(dataCC_qrtr[:, :, i], axis=1)
    anom1 = np.squeeze(dataCC_qrtr[:, :, i]) - np.expand_dims(clmn1, axis=1)
    anomCC_qrtr[:, :, i] = anom1

# reshape
anomCC_mtrx = np.zeros([nd, num_yrs*num_qrtr])*np.nan
dateCC_vec = np.zeros(num_yrs*num_qrtr, dtype=('datetime64[D]'))
mrkr = 0
for i in range(num_yrs):
    for j in range(num_qrtr):
        anomCC_mtrx[:, mrkr] = anomCC_qrtr[:, i, j]
        dateCC_vec[mrkr] = '{}-{:02d}-15'.format(yrs[i], mon_qrtr[j])
        mrkr = mrkr + 1

da1CC_qrtr = xr.DataArray(anomCC_mtrx, coords=[depthCC, dateCC_vec.astype('datetime64')], dims=['depth', 'time'])


# ------------------------------------------------------------------------ PLOT
# --make figure
sl = sttn_lst[1]
ttl1 = sttn_create_lbl(ds1CC, sttn_pre[1], sl)

# save dataarray in a list
da = [da0A, da1CC_qrtr]

# figure subplots, titles, and limits
nr_nc = [4, 1]
ttl = [ttl0, ttl1]
ylm = [dpth_bgn, dpth_end]
xlm = [yr_clim_bgn, yr_clim_end]

# contour levels
nlvl1 = np.arange(-30, 30.25, 0.25)
nlvl2 = np.ones(3)*np.nan
nlvl2 = np.array([-0.5, 0.5])

# open new figure
plt.close()
fig = plt.figure()

# number of xtick marks, based on IEA window size in years
wndw = 5

# subset da by these x and y limits
dpth_bgn = ylm[0]
dpth_end = ylm[1]
yr_clim_bgn = xlm[0]
yr_clim_end = xlm[1]


# mark missing data with symbol on top axis
mark_missing_list = [mark_missing_NH, mark_missing_CC]

# --number of plots to make, depends on the number of xr.da
order_list = range(0, len(da))
num_order_list = len(order_list)

# --set up plot grid
nr = nr_nc[0]
nc = nr_nc[1]

gs1 = gridspec.GridSpec(nr, nc)
gs1.update(left=0.2, right=0.8, bottom=0.05,
           top=0.9,  wspace=0.3, hspace=0.35)


for iii in range(0, num_order_list):
    print([iii, num_order_list])
    ax = plt.subplot(gs1[iii])

    # --get pd.df unique to the order number and create datetime index
    dao = da[iii]
    [dim1, dim2] = np.shape(dao)

    # --shorten by depth intreval
    in_dpth = ((dao['depth'] >= dpth_bgn) & (dao['depth'] <= dpth_end)).data

    # time dimension
    tt = dao['time'].to_dataframe()
    [dim_time, num_col] = tt.shape

    # check the shape of dao, should be depth x time
    if dim1 == dim_time:
        dao = dao.transpose()

    # --shorten by yrs_clim
    # in_yr = ((tt.index.year >= yr_clim_bgn) &
    #          (tt.index.year <= yr_clim_end))

    # --subset by depth and time
    # dao = dao[in_dpth, in_yr]

    # remove months with no samples
    datao = dao.data
    datao_mn = np.nanmean(datao, axis=0)

    in_data = np.where(datao_mn > -6)[0]

    # force the day to be one
    dt = dao.time[in_data]
    dt_dic = {'year': dt.time.dt.year.data,
              'month': dt.time.dt.month.data, 'day': np.ones(dt.shape)}
    tt = pd.to_datetime(dt_dic).values

    plt.contourf(tt, dao.depth.data, datao[:, in_data],
                 nlvl1, cmap=clr_map, vmin=-1*nlvl, vmax=nlvl)

    # x-limits set to clim interval
    xlm1 = np.datetime64(str(yr_clim_bgn)+'-01')
    xlm2 = np.datetime64(str(yr_clim_end)+'-12')
    xlm2 = np.datetime64(str(yr_clim_end+1)+'-01-01')
    # xlm2 = tt[-1]
    plt.xlim([xlm1, xlm2])

    # y-limits to depth interval
    plt.ylim([dpth_bgn, dpth_end])

    # --reverse the y-axis, so surface at the top
    plt.gca().invert_yaxis()

    # set xticks, xticklabels at start and end of the window
    yr_all = yr_clim_end - yr_clim_bgn
    if yr_all < 38:
        dt_intrvl = -1*(wndw-1)
    else:
        dt_intrvl = -1*(2*wndw-1)
    xt1 = np.flip(np.arange(yr_clim_end, yr_clim_bgn, dt_intrvl), 0)
    xt1_dt = [np.datetime64(xt1[j].astype('str')+'-01')
              for j in range(0, np.size(xt1))]

    # --xticks and x-label
    if (iii == num_order_list-1):
        plt.xticks(xt1_dt)
        plt.xlabel('Year')
    else:
        plt.xticks(xt1_dt, '')
    minor_locator = AutoMinorLocator(4)
    ax.xaxis.set_minor_locator(minor_locator)

    # make sure that only the year gets written as xticklabels
    monthyearFmt = mdates.DateFormatter('%Y')
    ax.xaxis.set_major_formatter(monthyearFmt)

    # add y-label
    plt.ylabel(ylbl)

    # add title
    plt.title(ttl[iii])

    # add visual indicator of dates that were filled
    if mark_missing_list[iii] == 1:
        data_chck = np.copy(dao.data)
        ind = np.isfinite(data_chck)
        data_chck[ind] = 1
        data_chck_sum = np.nansum(data_chck, axis=0)
        ind_miss = np.where(data_chck_sum == 0)[0]
        for k in range(len(ind_miss)):
            plt.text(dao.time.data[ind_miss][k], -3, '\u25CF',
                     ha='center', va='bottom', fontsize=2)

    # colorbar
    # m = plt.cm.ScalarMappable(cmap='RdYlBu_r')
    m = plt.cm.ScalarMappable(cmap=clr_map)
    m.set_array(datao[:, in_data])
    m.set_clim(-1*nlvl, nlvl)
    ax_pos = ax.get_position()
    if (len(clrbr_lbl) == num_order_list):
        x_cb = ax_pos.x0 + ax_pos.width + (ax_pos.width)/20
        y_cb = ax_pos.y0
        y_hght = ax_pos.height
        cbaxes = plt.gcf().add_axes([x_cb, y_cb, 0.006, y_hght])
        plt.colorbar(m, cax=cbaxes, extend='both', label=clrbr_lbl[iii])
    elif (iii == num_order_list-1):
        ax_pos = ax.get_position()
        x_cb = ax_pos.x0 + ax_pos.width + (ax_pos.width)/20
        y_cb = ax_pos.y0 + ax_pos.height*(1/3.0)
        y_hght = (gs1.top - y_cb)*(2/3.0)
        cbaxes = plt.gcf().add_axes([x_cb, y_cb, 0.006, y_hght])
        plt.colorbar(m, cax=cbaxes, extend='both', label=clrbr_lbl[0])

# create plot output directory
dir_plots = '{}/{}/'.format(dir_plot_out, iea_yr)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_plots)
except OSError:
    if not os.path.isdir(dir_plots):
        raise

# --Save the figure
fn_fig = '{}oc_z_{}_{}_{}_contour_fill_NRT.png'.format(
    dir_plots, var_name[0], sttn_wnt[0], sttn_wnt[1])
plt.savefig(fn_fig, dpi=300, bbox_inches='tight')

# remove NH25 netcdf, keep this private
os.remove("./NH25_CTD.nc")

import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt



# turn of toolbar on fig, set font
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 11})

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------

iea_yr = 2024

# Bakun UI
dir_out = './data_gha/bakunUI/'
dir_in = dir_out
fn_ui = ['{}UI_daily.nc'.format(dir_in),
         '{}UI_daily.nc'.format(dir_in)]

# dataset variables
ds_var = ['ui_day']
ds_coord = ['lat', 'year', 'days']

# lats to plot CUI
lat_list = [[45, 39, 33], [48, 42, 36]]
lat_list = [[48, 42, 36], [45, 39, 33]]

# IEA stuff
wndw = 5
yy_end = iea_yr
yr_clim_bgn = 1988
yr_clim_end = yy_end

# rolling mean size
roll = 30

# vertical lines at these dates, use 1981 as it is non-leapyear
dates_vline = ['1981-01-30', '1981-04-30', '1981-07-31', '1981-10-31']

xlbl = 'Yearday'

xtcks = np.arange(0, 390, 30)
xtcks[0] = 1

ylbl = ['UI (m^3/s/100 m coastline)', 'UI (m^3/s/100 m coastline)']

# --plot directory
dir_plot_out = './figures_gha/bakunUI_ts/'
fig_type = '.png'

# figure size
fig_wdth = 11
fig_hght = 8.5
# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

# create plot output directory
dir_plots = '{}{}/'.format(dir_plot_out, iea_yr)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_plots)
except OSError:
    if not os.path.isdir(dir_plots):
        raise

# size of input variables
num_data = len(fn_ui)
num_lat_wnt = len(lat_list[0])

# get yearday of vline dates
yearday_vline = (pd.to_datetime(dates_vline) -
                 pd.to_datetime('1981-01-01')).days + 1

# setup year vectors
yrs_clim = np.arange(yr_clim_bgn, yr_clim_end + 1)
yrs_wndw = np.arange(yr_clim_end - wndw + 1, yr_clim_end+1)
num_yrs_clim = len(yrs_clim)

# find the index of yy_end
in_yr_wnt = np.where(yrs_clim == yy_end)


# create daily matrix
dataD_mtrx = np.zeros([num_data, num_lat_wnt, num_yrs_clim, 366])*np.nan
dataD_run = np.zeros([num_data, num_lat_wnt, num_yrs_clim, 366])*np.nan
for i in range(0, num_data):
    print('i={}, fn={}'.format(i, fn_ui[i]))

    lat_wnt = lat_list[i]
    ds1 = xr.open_dataset(fn_ui[i])
    da1 = ds1[ds_var[0]]

    for j in range(num_lat_wnt):
        in_lat = np.where(ds1[ds_coord[0]] == lat_wnt[j])[0]
        da2 = np.squeeze(da1[in_lat, :])
        da2_run = da2.fillna(np.nanmean(da2.data)).rolling(time=roll, center=True).mean()

        for k in range(num_yrs_clim):
            in_yr = np.where(da2.time.dt.year == yrs_clim[k])[0]
            da3 = da2[in_yr]
            # da3_run = da2[in_yr].rolling(time=14, center=True).mean()
            da3_run = da2_run[in_yr]
            nt3 = da3.shape[0]
            dataD_mtrx[i, j, k, 0:nt3] = da3.data[0:366]
            dataD_run[i, j, k, 0:nt3] = da3_run.data[0:366]


climD = np.nanmean(dataD_run, axis=2)
stdD = np.nanstd(dataD_run, axis=2)

# subplots
num_clmn = num_data
num_row = len(lat_wnt)
gs1 = gridspec.GridSpec(num_row, num_clmn)


# change spacing
gs1.update(left=0.1, right=0.9, bottom=0.05,
           top=0.9, wspace=0.4, hspace=0.25)

plt.close()
fig = plt.figure(figsize=(fig_wdth, fig_hght))


# find common ylimit
ylm_mtrx = np.zeros([num_data, num_lat_wnt, 2])
for i in range(0, num_data):
    lat_wnt = lat_list[i]
    for j in range(num_lat_wnt):
        # var ending in "1" are 366 days long
        # var ending in "f" are 365 days long
        # data1 = np.squeeze(dataD_mtrx[i, j, in_yr_wnt, :])
        data1 = np.squeeze(dataD_run[i, j, in_yr_wnt, :])
        dataf = data1[0:365]
        ylm_mtrx[i, j, 0] = np.nanmin(dataf)
        ylm_mtrx[i, j, 1] = np.nanmax(dataf)


# make the time series plot
for i in range(0, num_data):
    lat_wnt = lat_list[i]
    for j in range(num_lat_wnt):
        # var ending in "1" are 366 days long
        # var ending in "f" are 365 days long
        # data1 = np.squeeze(dataD_mtrx[i, j, in_yr_wnt, :])
        data1 = np.squeeze(dataD_run[i, j, in_yr_wnt, :])
        dataf = data1[0:365]

        clim1 = np.squeeze(climD[i, j, :])
        climf = clim1[0:365]

        std1 = np.squeeze(stdD[i, j, :])
        stdf = std1[0:365]

        x1 = np.arange(1, 367)
        xf = np.arange(1, 366)

        # need one label for the std envelope, do this by combining them into
        # one line, however this puts an ugly vertical line at the end, so
        # for the std window use the whole 366 days and remove the vertical
        # line with change xlim to [1,365]
        x_std = np.append(x1, np.flipud(x1))
        y_std = np.append(clim1-std1, np.flipud(clim1+std1))

        # plot the three lines
        ax1 = plt.subplot(gs1[j, i])

        # data
        plt.plot(xf, dataf, '-k', label=yrs_clim[in_yr_wnt][0], linewidth=2)

        # clim
        plt.plot(xf, climf, '--', color='dodgerblue', linewidth=1.0, label='Clim')

        # std envelope
        plt.plot(x_std, y_std, '-', color='dodgerblue', linewidth=0.75, label='±1 s.d.', alpha=0.3)

        plt.fill_between(x1, clim1-std1, clim1+std1,
                         facecolor=np.ones(3)*0.8, alpha=0.1)
        # horizontal legend
        if j == 0 and i == 1:
            lgnd = ax1.legend(loc='upper center', ncol=int((wndw+1)/2),
                              bbox_to_anchor=(-0.25, 1.4), edgecolor='none',
                              framealpha=0, fontsize=12)

        # xticks and xlabel
        ax1.set_xticks(xtcks)
        if j < num_lat_wnt-1:
            ax1.set_xticklabels('')
        else:
            plt.xlabel(xlbl)

        # xlim
        plt.xlim([1, 365])
        xrng = plt.xlim()
        dx = np.diff(xrng)/20.0

        # ylim
        ymin = np.min(ylm_mtrx)
        ymax = np.max(ylm_mtrx)

        yrng = [ymin, ymax]
        dy = np.diff(yrng)/10
        ylm = [yrng[0]-dy, yrng[1]+dy]
        plt.ylim(ylm)
        plt.ylabel(ylbl[i])

        # --remove top and right axes lines
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        # vlines for set dates
        ax1.set_axisbelow(False)
        plt.vlines(yearday_vline, plt.ylim()[0], plt.ylim()[1],
                   colors='black', linestyles=':', linewidth=0.5, zorder=30)
        for jjj in range(0, len(dates_vline)):
            pddt = pd.to_datetime(dates_vline[jjj])
            txt = '{}/{}'.format(pddt.month, pddt.day)
            plt.text(yearday_vline[jjj], ylm[1]+dy, txt,
                     horizontalalignment='center',
                     verticalalignment='top',
                     color='black', fontsize=8)

        # Latitude Text
        xtxt = plt.xlim()[1]+0.5*dx
        ytxt = plt.ylim()[1]-1.5*dy
        tbox = plt.text(xtxt, ytxt, '{} °N'.format(lat_wnt[j]),
                        fontsize=15, horizontalalignment='right',
                        fontweight='bold',
                        bbox={'facecolor': 'none', 'edgecolor': 'none'})

# create output name and save xarray dataset to netcdf file
fn_out = '{}bakunUI_daily_{}roll_same_ylim.png'.format(
    dir_plots, roll)

plt.savefig(fn_out, dpi=300, bbox_inches='tight')

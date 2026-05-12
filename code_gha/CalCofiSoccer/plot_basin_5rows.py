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

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# end year
iea_yr = 2026

# input file
dir_out = './data_gha/CalCofiSoccer/'
dir_in = dir_out
fn_in = '{}/CalCofiSoccer/basin_5rows.nc'.format(dir_in)

# directory out
dir_plot_out = './figures_gha/CalCofiSoccer/'
# dir_plot_out = './figures_x13/CalCofiSoccer'


# plot labels for the different variables
plt_lbl = ['a', 'b', 'c', 'd', 'e']
ttl_lbl = ['ONI', 'PDO', 'NPGO', 'North Pacific\nSST anom (°C)', 'CCLME\nSST anom (°C)']

var_wnt = ['ONI', 'PDO', 'NPGO', 'NEP', '150km']
roll_vec = [1, 5, 15, 31]

# ylim
ylm_wnt = [[-2, 3], [-3, 4], [-4, 3], [-1, 1], [-2, 2]]

# last x years
x_yrs = 4

# figure size
fig_wdth = 8.5
fig_hght = 6

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# len input vars
num_var_wnt = len(var_wnt)
num_roll_wnt = len(roll_vec)


for iii in range(num_roll_wnt):
    roll_wnt = roll_vec[iii]

    ds1 = xr.open_dataset(fn_in).sel(day_roll=roll_wnt)
    yr_bgn = ds1.time.dt.year[0].data
    yr_end = ds1.time.dt.year[-1].data


    # CalCOFI reports usually focus on the last x years, get these year
    yrs_last = np.arange(yr_end-x_yrs+1, yr_end+1)
    # setup subplots, spacing and figure size
    plt.close()
    gs1 = gridspec.GridSpec(num_var_wnt, 1)
    gs1.update(left=0.1, right=0.8, bottom=0.05, top=0.9, wspace=0.1, hspace=0.1)
    fig = plt.figure(figsize=(fig_wdth, fig_hght))
    rcParams.update({'font.size': 7})

    # subplots, rows = erddap variables, columns = 1
    gs1 = gridspec.GridSpec(num_var_wnt, 1)
    gs1.update(left=0.1, right=0.65, bottom=0.05,
               top=0.95, wspace=0.0, hspace=0.2)

    # loop over variables

    for i in range(0, num_var_wnt):
        da1 = ds1[var_wnt[i]]


        # print the last dates
        print(da1.time.data[-1])


        dates = da1.time.astype('datetime64[M]')
        Y = da1.data


        ax = plt.subplot(gs1[i])
        plt.plot(dates, Y, linewidth=0.2, color='black')
        ax.fill_between(dates, 0, Y, where=Y >= 0, color='red',
                        linewidth=0.1, interpolate=True)
        ax.fill_between(dates, 0, Y, where=Y <= 0, color='blue',
                        linewidth=0.1, interpolate=True)

        # xlimits
        # plt.xlim(dates[0], dates[-1])
        xlm1 = np.datetime64('{}-01-01'.format(yr_bgn))
        xlm2 = np.datetime64('{}-02-01'.format(yr_end+1))
        plt.xlim([xlm1, xlm2])

        # --get rid of the top and right
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # y labels
        plt.ylabel(ttl_lbl[i])

        # --set format of the dates for the x-axis
        if i < num_var_wnt-1:
            xfmt = mdates.DateFormatter('')
        else:
            # xfmt = mdates.DateFormatter('%Y-%m')
            xfmt = mdates.DateFormatter('%Y')
        ax.xaxis.set_major_formatter(xfmt)
        # ax.set_xticks(dates[0::12 * 5])
        # plt.xticks(rotation=0)

        # --add x-minor tick marks at intervals of 5, this will produce
        # --a tick mark at Jan of every year
        minorLocator = AutoMinorLocator(5)
        ax.xaxis.set_minor_locator(minorLocator)

        minorLocator = AutoMinorLocator(2)
        ax.yaxis.set_minor_locator(minorLocator)

        # --x&y tick fontsizes
        # mpl.rc('xtick', labelsize=8)
        # mpl.rc('ytick', labelsize=8)

        # --vertical lines for the last x years
        plt.ylim(ylm_wnt[i])
        ylm = ax.get_ylim()
        for j in range(0, len(yrs_last)):
            datej = '{}-01-01'.format(yrs_last[j])
            plt.vlines(np.datetime64(datej),
                       ylm[0], ylm[1], color='gray',
                       alpha=0.3, linewidth=0.5)

        # plot labels
        xlm = ax.get_xlim()
        plt.text(xlm[0]+np.diff(xlm)/25, ylm[1] -
                 np.diff(ylm)/25, plt_lbl[i], fontsize=12)

    # create plot output directory
    dir_plots = '{}/{}/'.format(dir_plot_out, iea_yr)

    # check if directory exist, if it doesn't then create
    try:
        os.makedirs(dir_plots)
    except OSError:
        if not os.path.isdir(dir_plots):
            raise

    # --Save the figure
    fn_fig = '{}blue_red_basin_5rows_roll{}.png'.format(
        dir_plots, roll_wnt)
    plt.savefig(fn_fig, dpi=300, bbox_inches='tight')

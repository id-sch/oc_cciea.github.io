import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import xarray as xr
import matplotlib.gridspec as gridspec


# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Times New Roman']})

# turn of toolbar on fig, set font
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 12})

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2024
# iea_yr = 2025

# --basin index wanted
# set lat range
lat_bgn = 33
lat_end = 48
dlat = 3


pheno_wnt = ['sti', 'lusi', 'tumi']
num_basin = len(pheno_wnt)

# --input directory
dir_out = './data_gha/bakunUI/'
# dir_out = './data_x13/bakunUI/'
dir_in = dir_out

# variable name in the xr.ds
var_lbl = ['cui_mtrx']

# dimensions labels
dim_lbl = ['lat', 'year', 'days']

# --IEA file names
file_pre = 'oc'
num_pre = len(file_pre)

# --IEA year clim
yr_bgn = 1967
yr_end = iea_yr

# --IEA window size
num_xyrs = 5

# vertical lines at these dates, use 1981 as it is non-leapyear
dates_vline = ['1981-01-30', '1981-04-30', '1981-07-31', '1981-10-31']

# figure usually has 2 columns
num_clmn = 2

# xticks
xtck = np.arange(0, 365, 30)

# get rid of big numbers in ylabels, set scale factor
sf = 10000

# color for the xyrs
color_xyrs = ['limegreen', 'orange', 'dodgerblue', 'red', 'black']

# xlabel
xlbl = 'Yearday'

# ylabels
ylbl = 'CUI x10⁻⁴\n(m³/s/100 m coastline)'

# figure size
fig_wdth = 8.5
fig_hght = 11

# --plot directory
dir_plot_out = './figures_gha/phenologyUI/'
# dir_plot_out = './figures_x13/bakunUI/'
fig_type = '.png'


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

# lat position, start with northern lats at top of figure
lat_pos = np.arange(lat_end, lat_bgn-dlat, -1*dlat)
num_lat_pos = len(lat_pos)

# get yearday of vline dates
yearday_vline = (pd.to_datetime(dates_vline) -
                 pd.to_datetime('1981-01-01')).days + 1

# open the UI data
fn_data = '{}/UI_cui_mtrx.nc'.format(dir_in)
ds = xr.open_dataset(fn_data)

# get years clim
in_yr = np.logical_and(ds[dim_lbl[1]] >= yr_bgn, ds[dim_lbl[1]] <= yr_end)

# get the cui_mtrx as a xr.da
da = ds[var_lbl[0]][:, in_yr, :]

# da dimensions
lat = da[dim_lbl[0]]
year = da[dim_lbl[1]]
days = da[dim_lbl[2]]

# get the last x_yrs
xyrs = np.arange(year[-1]+1-num_xyrs, year[-1]+1, 1)

# setup subplots
num_row = int(num_lat_pos/num_clmn)
gs1 = gridspec.GridSpec(num_row, num_clmn)

# get the latitudes that appear in a row, use to find ylim
lat_row = lat_pos.reshape([num_row, num_clmn])

# change spacing
plt.close()
gs1.update(left=0.01, right=0.99, bottom=0.12, top=0.9, wspace=0.15, hspace=0.2)

# loop over all lat positions


fig = plt.figure(figsize=(fig_wdth, fig_hght))
for i in range(0, num_lat_pos):
    # get lat wanted
    in_lat = np.where(lat == lat_pos[i])[0]
    cui = np.squeeze(ds[var_lbl[0]][in_lat, :, :])/sf
    cui_mn = np.nanmean(cui, axis=0)

    # get the row that the lat appears in
    in_row = np.where(lat_pos[i] == lat_row)
    lats_in_row = np.squeeze(lat_row[in_row[0], :])
    max_vec = np.zeros(num_clmn)
    min_vec = np.zeros(num_clmn)
    for j in range(0, num_clmn):
        inj = np.where(lat == lats_in_row[j])[0]
        cuij = np.squeeze(ds[var_lbl[0]][inj, :, :])/sf
        maxj = np.nanmax(cuij)
        minj = np.nanmin(cuij)
        max_vec[j] = maxj
        min_vec[j] = minj

    yrng = [min(min_vec), max(max_vec)]
    dyrng = np.diff(yrng)/20
    ylm = [yrng[0]-dyrng, yrng[1]+dyrng]

    # ax = plt.subplot(gs1[i])
    ax = fig.add_subplot(gs1[i])

    # plot all years in gray
    plt.plot(days, cui.data.T, color=np.ones(3)*0.9, linewidth=0.5)

    # plot mean cui in black
    plt.plot(days, cui_mn, color='black', label='mean', linestyle='--')

    # loop and plot the last x_yrs
    for j in range(0, num_xyrs):
        in_xyr = np.where(year == xyrs[j])[0]
        plt.plot(days, np.squeeze(cui[in_xyr, :]),
                 color=color_xyrs[j], linewidth=1.5, label=xyrs[j].astype('str'))

    # --remove top and right axes lines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xticks(xtck)

    # remove xtick labels
    if (in_row[0][0] < num_row-1):
        # plt.xticks(xtck, ' ')
        ax.set_xticklabels('')

    # remove ytick labels
    if (in_row[1][0] == num_clmn-1):
        plt.yticks(plt.yticks()[0], '')

    # xlabels
    if (in_row[0][0] == num_row-1):
        plt.xlabel(xlbl)

    # ylabels
    if (in_row[1][0] != num_clmn-1):
        plt.ylabel(ylbl)

    # xlim
    plt.xlim([1, 365])

    # ylim
    plt.ylim(ylm)
    print('i={}, ylm0 {} ylm1{}'.format(i, ylm[0], ylm[1]))

    # vlines for set dates
    ax.set_axisbelow(False)
    plt.vlines(yearday_vline, ylm[0], ylm[1],
               colors='red', linestyles=':', linewidth=0.5, zorder=30)
    for j in range(0, len(dates_vline)):
        pddt = pd.to_datetime(dates_vline[j])
        txt = '{}/{}'.format(pddt.month, pddt.day)
        plt.text(yearday_vline[j], ylm[1], txt,
                 horizontalalignment='center', color='red', fontsize=5, verticalalignment='bottom')

    # Latitude Text
    xtxt = plt.xlim()[0]
    ytxt = plt.ylim()[1]
    tbox = plt.text(xtxt, ytxt, '{}$\degree$N'.format(lat_pos[i]),
                    fontsize=9, horizontalalignment='center',
                    bbox={'facecolor': 'white', 'edgecolor': 'white'})

# add legend on top
plt.legend(bbox_to_anchor=(-0.8, 3.52, 1.5, .102), loc=3,
           ncol=num_xyrs+1, mode="expand", borderaxespad=0.,
           frameon=None)

# save
fn_out = '{}oc_CUI_bakun.png'.format(dir_plots)
plt.savefig(fn_out, dpi=300, bbox_inches='tight')

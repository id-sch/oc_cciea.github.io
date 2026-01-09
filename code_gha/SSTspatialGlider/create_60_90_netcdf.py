import os
import shutil
import requests
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec


mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 11})
mpl.rcParams.update({'savefig.bbox': 'tight'})


# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------
# iea year
iea_yr = 2025

# years to plot
yr_bgn = 2007
yr_end = iea_yr

# download and plot thes
line_wnt = [66, 80, 90]
fn1_tmplt  = './tmp/anomaly_z_{}.nc'

var_wnt = 'temperature'
dis_wnt = np.arange(0, 205, 5)

nlvl = np.arange(-3.0, 3.25, 0.25)

clr_bar_lbl = 'Anomaly of Temperature (°C)'

# --plot directory
dir_plot_out = './figures_gha/SSTspatialGlider/'
fig_type = '.png'

# figure size
fig_wdth = 11
fig_hght = 8.5

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# len input variables
num_line_wnt = len(line_wnt)

# create tmp directory, will save the download data here
dir_tmp = 'tmp'
# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_tmp)
except OSError:
    if not os.path.isdir(dir_plots):
        raise

# create plot output directory
dir_plots = '{}{}/'.format(dir_plot_out, iea_yr)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_plots)
except OSError:
    if not os.path.isdir(dir_plots):
        raise

# years wanted
yrs_wnt = np.arange(yr_bgn, yr_end+1)
xlm = [np.datetime64('{}-01-01'.format(yr_bgn)), np.datetime64('{}-12-31'.format(yr_end))]

# download and put dataset in a list, find common years
ds_list = []
for i in range(num_line_wnt):
    url = 'https://spraydata.ucsd.edu/data/cugn/clim/netcdf/anomaly_z_{}.nc'.format(line_wnt[i])

    fn1 = fn1_tmplt.format(line_wnt[i])
    ddd = requests.get(url, timeout=20)
    open(fn1, 'wb').write(ddd.content)

    ds1 = xr.open_dataset(fn1)
    ds_list.append(ds1)

da1_list = []
for i in range(num_line_wnt):
    da1 = ds_list[i][var_wnt].sel(distance=dis_wnt).mean('distance').sel(time=ds_list[i].time.dt.year.isin(yrs_wnt))    
    da1_list.append(da1)


# end year is even
if np.mod(yr_end, 2) == 0:
    xend = yr_end
    xbgn = yr_bgn + 1

# end year is odd
if np.mod(yr_end, 2) == 1:
    xend = yr_end
    xbgn = yr_bgn


xtcks = np.arange(xbgn, xend + 1, 2)
xtcks_dt = []
for i in range(len(xtcks)):
    xtcks_dt.append(np.datetime64('{}-01-01'.format(xtcks[i])))

# open new figure of set size
plt.close()
fig = plt.figure(figsize=(fig_wdth, fig_hght))


# subplots using gridspec
num_row = num_line_wnt
num_clmn = 1
gs1 = gridspec.GridSpec(num_row, num_clmn)
gs1.update(left=0.1, right=0.9, bottom=0.1,
           top=0.9, wspace=0.1, hspace=0.4)


for i in range(num_line_wnt):
    # subplot axes
    ax1 = plt.subplot(gs1[i, 0])

    tt = da1_list[i].time.data
    z = da1_list[i].depth.data
    data1 = da1_list[i].data
    plt.contourf(tt, z, data1, nlvl, extend='max', cmap='RdBu_r')

    ax1.invert_yaxis()

    plt.title('Line {}, {} - {} km'.format(line_wnt[i], dis_wnt[0], dis_wnt[-1]))
    plt.ylabel('Depth (m)')

    if i < num_line_wnt - 1:
        ax1.set_xticklabels('')
    else:
        plt.xlabel('Time')
    cb1 = plt.colorbar()
    cb1.set_label(clr_bar_lbl)

    yrs = np.unique(da1_list[i].time.dt.year.data)
    yr_bgn = yrs[0]
    yr_end = yrs[-1]

    ax1.set_xticks(xtcks_dt)
    ax1.set_xticklabels(xtcks)
    plt.minorticks_on()
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
    plt.xlim(xlm)


line_str = '_'.join(list(map(str,line_wnt)))
fn_plot = '{}oc_z_{}_gldr_lines_{}_distance_{}_{}km.png'.format(dir_plots, var_wnt, line_str, dis_wnt[0], dis_wnt[-1])
plt.savefig(fn_plot)

# remove the directory and files that has the downloaded
shutil.rmtree(dir_tmp)

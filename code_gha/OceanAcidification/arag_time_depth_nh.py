import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
from scipy.interpolate import griddata


mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 12})


# ----------------------------------------------------------------------
# BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2025

file_pre = 'oc_arg_Newport'

# --station wanted, two stations newport hydrographic at 5 and 25 km
sttn_wnt = ['NH25_CTD', 'NH05_CTD']
num_sttn = len(sttn_wnt)

# --input directory
dir_in = ['./',
          './']

nlvl1 = np.arange(0, 3.1, 0.1)

ylm_list = [[-180, -2], [-50, -2]]

yr_bgn = 1998
yr_end = iea_yr
time_bgn = np.datetime64('{}'.format(yr_bgn, 'Y'))
time_end = np.datetime64('{}'.format(yr_end, 'Y'))
xtck = np.arange(time_bgn, time_end+2, 2, dtype='datetime64[Y]')
ttl_lbl = ['NH25', 'NH05']

# --plot directory
dir_plot_out = './figures_gha/OceanAcidification/'
fig_type = '.png'

fig_wdth = 11
fig_hght = 8.5

# ----------------------------------------------------------------------
# END: Change These
# ----------------------------------------------------------------------
num_sttn = len(sttn_wnt)

# Aragonite paramaters
do_ref = 140
temp_ref = 8
alpha0 = 0.9242
alpha1 = 0.004492
alpha2 = 0.00094

# create plot output directory
dir_plots = '{}{}/'.format(dir_plot_out, iea_yr)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_plots)
except OSError:
    if not os.path.isdir(dir_plots):
        raise

# subplots using gridspec
num_row = 4
num_clmn = 1
gs1 = gridspec.GridSpec(num_row, num_clmn)
gs1.update(left=0.1, right=0.9, bottom=0.1,
           top=0.9, wspace=0.1, hspace=0.1)

# --Newport Data
for i in range(num_sttn):
    # create new figure and close old one
    plt.close()
    fig = plt.figure(figsize=(fig_wdth, fig_hght))

    # open data
    fn_in = dir_in[0] + '/' + sttn_wnt[i] + '.nc'
    ds0 = xr.open_dataset(fn_in)
    ds0M = ds0.resample(time='1MS').mean()

    nt, nd = ds0M['Temperature'].shape

    # calc arag1
    do2 = ds0M['Oxygen']*44.661
    temp1 = ds0M['Temperature']
    arag_c = alpha0+alpha1*(do2-do_ref)+alpha2*(temp1-temp_ref)*(do2-do_ref)

    # griddata
    d1 = arag_c.data.T
    nz = arag_c.depth.shape[0]
    nt = arag_c.time.shape[0]

    t1 = np.arange(0, nt)
    z1 = arag_c.depth.data
    t1g, z1g = np.meshgrid(t1, z1)

    ind = np.isfinite(d1)
    d1g = griddata((t1g[ind], z1g[ind]), d1[ind], (t1g,z1g), method='linear')

    inm = np.where(d1g > nlvl1[-1])
    d1g[inm] = nlvl1[-1]

    # contour
    # subplot axes
    ax1 = plt.subplot(gs1[0, 0])

    cs = plt.contourf(arag_c.time, arag_c.depth*-1, d1g, nlvl1, cmap='jet')
    plt.colorbar(cs)
    plt.contour(arag_c.time, arag_c.depth*-1, d1g, [1], colors='black')

    plt.ylim(ylm_list[i])

    ax1.set_xticks(xtck)
    plt.xlim([time_bgn, time_end+1])

    xfmt = mdates.DateFormatter('%Y')
    ax1.xaxis.set_major_formatter(xfmt)
    ax1.xaxis.set_minor_locator(AutoMinorLocator(2))

    plt.title('{} Aragonite Saturation'.format(ttl_lbl[i]))

    fn_fig = '{}{}_time_depth_{}.png'.format(dir_plots, file_pre, ttl_lbl[i])
    plt.savefig(fn_fig)

# remove the directory and files that has the downloaded
os.remove("./NH05_CTD.nc")
os.remove("./NH25_CTD.nc")

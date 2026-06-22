import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import calendar as clndr
import pandas as pd
from matplotlib import gridspec
from matplotlib import rcParams
from matplotlib import interactive
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
interactive(True)


class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()

# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------
# .1) sst time series
iea_yr = 2025
# iea_yr = 2026

# clim year period
yr_clim_bgn = 1982
yr_clim_end = 2010


# lats of the sst data
lat_wnt = [31, 47]

# file name of NOAA OI sst
# fn_sst = './oc_cciea.github.io/data_gha/HCI/sst_oi_lat_30_48_xdis_150km.nc'
fn_sst = './data_gha/HCI/sst_oi_lat_30_48_xdis_150km.nc'

# file name of the BEUTI an CUTI 
# fn_beuti = './oc_cciea.github.io/data_gha/newUI/BEUTI_monthly.nc'
# fn_cuti = './oc_cciea.github.io/data_gha/newUI/CUTI_monthly.nc'
fn_beuti = './data_gha/newUI/BEUTI_monthly.nc'
fn_cuti = './data_gha/newUI/CUTI_monthly.nc'


# files list
file_list = [fn_sst, fn_beuti, fn_cuti]

# roms variable name in the xr.ds
var_sst = ['sst_oi']

# beuti and cuti variable name in the xr.ds
var_ui = ['ui_mon']

# var list
var_list = [var_sst, var_ui, var_ui]
var_lbl = ['SST', 'BEUTI', 'CUTI']

# sst_oi dimensions labels
dim_sst = ['time', 'latitude', 'longitude']

# beuti and cuti dimensions labels
dim_ui = ['time', 'lat']

# dims list
dim_list = [dim_sst, dim_ui, dim_ui]

# figure usually has 1 columns, 2 rows
num_clmn = 1
num_row = 2

# ylim and ticks
ylm = [30, 48]
y_tck = np.arange(ylm[0], ylm[1]+3, 3)

# CalCOFI figures usually focus discussion over the last 3 years, change
# if more years are wanted
num_xyrs = 3

# figure size
fig_wdth = 8.5
fig_hght = 6

# nlevels
d_sst = 0.5
dmin_sst = -2
dmax_sst = 2
nlvl1_sst = np.arange(dmin_sst-3*d_sst, dmax_sst+4*d_sst, d_sst)
nlvl2_sst = np.arange(dmin_sst-1*d_sst, dmax_sst+2*d_sst, 2*d_sst)

d_beuti = 3
dmin_beuti = -18
dmax_beuti = 18
nlvl1_beuti = np.arange(dmin_beuti-6*d_beuti, dmax_beuti+7*d_beuti, d_beuti)
nlvl2_beuti = np.arange(dmin_beuti-1*d_beuti, dmax_beuti+2*d_beuti, 2*d_beuti)

d_cuti = 0.2
dmin_cuti = -1.2
dmax_cuti = 1.2
nlvl1_cuti = np.arange(dmin_cuti-6*d_cuti, dmax_cuti+7*d_cuti, d_cuti)
nlvl2_cuti = np.arange(dmin_cuti-1*d_cuti, dmax_cuti+2*d_cuti, 2*d_cuti)

min_max_list = [[dmin_sst, dmax_sst], [dmin_beuti, dmax_beuti], [dmin_cuti, dmax_cuti]]
nlvl1_list = [nlvl1_sst, nlvl1_beuti, nlvl1_cuti]
nlvl2_list = [nlvl2_sst, nlvl2_beuti, nlvl2_cuti]

# colorbar labelsb
clrbr_lbl_sst = 'SST Anoms ($\degree$C)'
clrbr_lbl_beuti = 'BEUTI Anoms (mmol s$^{-1}$ m$^{-1}$)'
clrbr_lbl_cuti = 'CUTI (m$^2$ s$^{-1}$)'

clrbr_list = [clrbr_lbl_sst, clrbr_lbl_beuti, clrbr_lbl_cuti]

# .5) plot directory
dir_plot_out = './figures_gha/miniESR/'
# dir_plot_out = './figures_x13/miniESR/'

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# input variable size
num_file_list = len(file_list)


yrs_clim = np.arange(yr_clim_bgn, yr_clim_end+1)
yr_bgn = iea_yr - num_xyrs + 1
yrs_wnt = np.arange(yr_bgn, iea_yr + 1)

# vlines
vln = np.zeros(num_xyrs - 1, dtype='datetime64[M]')
for i in range(0, num_xyrs-1):
    vln[i] = '{}-01'.format(yr_bgn+1+i)

# setup subplots, spacing and figure size
gs1 = gridspec.GridSpec(num_row, num_clmn)
gs1.update(left=0.05, right=0.85, bottom=0.05, top=0.9, wspace=0.1, hspace=0.1)

for iii in range(num_file_list):
    ds1 = xr.open_dataset(file_list[iii])

    if iii == 0:
        # need to area average the sst
        da1 = ds1[var_list[iii][0]].mean(dim_list[iii][2]).T
    else:
        da1 = ds1[var_list[iii][0]]

    da1_clim = da1.sel(time=da1.time.dt.year.isin(yrs_clim)).groupby('time.month').mean('time')
    da1_anom = da1.groupby('time.month') - da1_clim

    da1f = da1_anom.sel(time=da1_anom.time.dt.year.isin(yrs_wnt))

    # new figure
    plt.close()
    fig = plt.figure(figsize=(fig_wdth, fig_hght))

    # contour sst anom
    ax = fig.add_subplot(gs1[0])
    x = da1f[dim_list[iii][0]].data.astype('datetime64[M]')
    y = da1f[dim_list[iii][1]].data
    CS = plt.contour(x, y, da1f.data, nlvl2_list[iii], vmin=min_max_list[iii][0], vmax=min_max_list[iii][1], colors='gray', linewidths=0.5)

    # Recast levels to new class
    CS.levels = [nf(val) for val in CS.levels]

    labels1 = plt.clabel(CS, CS.levels, inline=False, fmt='%r', fontsize=7, colors='k', rightside_up=True)

    dmin = np.ceil(np.nanmin(da1f.data))
    dmax = np.ceil(np.nanmax(da1f.data))
    chck_min = 0
    chck_max = 0
    if dmin < min_max_list[iii][0]:
        chck_min = 1
        extnd1 = 'min'
    if dmax > min_max_list[iii][1]:
        chck_max = 1
        extnd1 = 'max'
    if chck_max+chck_min==2:
        extnd1 = 'both'

    plt.contourf(x, y, da1f.data, nlvl1_list[iii], vmin=min_max_list[iii][0], vmax=min_max_list[iii][1], cmap='bwr')

    # xtick labels
    mon_lbl = list()
    for i in x:
        dti = pd.to_datetime(i)
        moni = clndr.month_name[dti.month]
        if moni == 'January':
            mon_lbl.append('{}\n      {}'.format(moni[0], dti.year))
        else:
            mon_lbl.append(moni[0])
    plt.xticks(x, mon_lbl, fontsize=7)

    # ytick labesl
    lat_lbl = list()
    for i in y_tck:
        lat_lbl.append('{}$\degree$N'.format(i))
    plt.yticks(y_tck, lat_lbl, fontsize=7)
    plt.ylim([y[0], y[-1]])

    # vlines
    plt.vlines(vln, ylm[0], ylm[1], colors=np.array([1, 1, 1])*0.4, linestyles='dashed', linewidth=2)

    # colorbar
    ax_pos = ax.get_position()
    x_cb = ax_pos.x0 + ax_pos.width + ax_pos.width/40.0
    y_cb = ax_pos.y0
    y_hght = ax_pos.height
    cbaxes = plt.gcf().add_axes([x_cb, y_cb+(y_hght*0.1)/2.0, 0.005, y_hght*0.9])

    m = plt.cm.ScalarMappable(cmap='bwr')
    m.set_array(da1f.data)
    m.set_clim(min_max_list[iii][0], min_max_list[iii][1])
    dlvl = np.diff(nlvl1_list[iii])[0]
    bound = np.arange(min_max_list[iii][0], min_max_list[iii][1]+dlvl, dlvl)
    plt.colorbar(m, cax=cbaxes, label=clrbr_list[iii],
                 format='%4.1f', extend=extnd1,
                 boundaries=bound)

    # create plot output directory
    dir_year = '{}/{}/'.format(dir_plot_out, iea_yr)
    dir_plots = dir_year + 'hovmoller/'

    # check if directory exist, if it doesn't then create
    try:
        os.makedirs(dir_plots)
    except OSError:
        if not os.path.isdir(dir_plots):
            raise

    # --Save the figure
    fn_fig = '{}oc_miniESR_{}_time_depth_yrs_{}_{}_contour.png'.format(dir_plots, var_lbl[iii], yr_bgn, iea_yr)
    plt.savefig(fn_fig, dpi=300, bbox_inches='tight')

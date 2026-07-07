import os
import calendar
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import calendar as clndr
from matplotlib import rcParams
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


# plot paramaters
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'RdYlBu_r',
    'axes.grid': False,
    'savefig.dpi': 300,  # to adjust notebook inline plot size
    'xtick.top':        False,  # shold the top and bottom have tick marks
    'xtick.bottom':     True,
    'xtick.major.size': 0,
    'ytick.major.size': 0,
    'ytick.direction': 'out',
    'xtick.direction': 'out',
    'axes.labelsize': 10,  # fontsize for x and y labels 
    'axes.titlesize': 10,
    'font.size': 10,  # was 10
    'legend.fontsize': 10,  # was 10
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'figure.figsize': [8.0, 8.5],
    'font.family': 'STIXGeneral',
    'mathtext.fontset': 'stix',
    'toolbar': 'None',
    'savefig.bbox': 'tight'
}
rcParams.update(params)

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# end year

# regions
rgn = np.arange(1, 5)

# --input directory
dir_out = './data_gha/HCI_ROMS/'
# dir_out = './data_x13/HCI_ROMS/'
dir_in = dir_out

# depth level
z_wnt = -2

# variable wanted
var_wnt = ['temp']

# dimensions
dim_wnt = ['latitude', 'longitude']

# File1 (ROM2 temp): x distance and lat range list
var_roms = 'temp_mtrx'
var_level = 'lvlM_mtrx'    # 1=below threshold, 0=above threshold
threshold_wnt = 'mn_75kmM'  # M=monthly threshold, S=seasonal threshold
xdis = 150
lat_rgn = [[43.5, 48], [40, 43.5], [35.5, 40], [30, 35.5]]
rgn_wnt = [1, 2, 3, 4]

# File2 (seasonal thresholds): x distance threshold
xdis_thrsh = 75

# contour levels
d_min = 7
d_max = 22
dx = 1
nlvl1 = np.arange(d_min-4*dx, d_max+8*dx, dx)

cbr_tcks = np.arange(d_min, d_max+3, 3)

# lvl_wnt dimension name of the xr.ds
var_dim = 'temp cutoff'

# time series name
var_vec = 'num_below_cutoff'

# map grid
xlm_list = [[-126.5, -122.5],
            [-126.5, -123.0],
            [-126.5, -120.5],
            [-123.5, -115.3]]
xg1_list = [np.arange(-127, -120, 1),
            np.arange(-127.0, -119.0, 1),
            np.arange(-127.0, -119.0, 1),
            np.arange(-126.0, -114.0, 1)]
yg1_list = [np.arange(43, 49, 1),
            np.arange(39.5, 44, 1),
            np.arange(35, 42, 1),
            np.arange(30, 37, 1)]

xg2_list = [xg1_list[0][0::2], xg1_list[1][1::2], xg1_list[2][2::3], xg1_list[3][1::3]]
yg2_list = [yg1_list[0][1::2], yg1_list[1][1::2], yg1_list[2][1::2], yg1_list[3][1::2]]

# colorbar label
clrbr_lbl = 'Temp ($\degree$C)'

# rows and columns of figure
num_rows = 6
num_clmn = 8

# months
mons = np.arange(1, 13)

# lat and lon position of text labels 
yr_lbl = [[-122.55, 47.3],[-123.1, 42.7],[-120.9, 39.2],[-115.4, 34.4]]
pnts_lbl = [[-123.15, 46.0], [-123.6, 41.7], [-121.4, 37.6], [-116.3, 32.6]]
lgnd_pad = [0.35, 0.35, 0.15, 0.25]

# figure size
fg_sz = [[8.0, 10], [9.75, 9.75], [9.25, 8.5], [11.25, 9.0]]

# annotate hci based on sd
lbl_sd = ['Low Compression', 'Medium Compression', 'High Compression']
clr_sd = ['springgreen', 'yellow', 'tomato']


# --IEA file names
file_pre = 'hci'

# --plot directory
dir_plot_out = './figures_gha/HCI_ROMS/'
# dir_plot_out = './figures_x13/HCI_ROMS/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# len of input variables
num_data = len(rgn)

# size of input variables
num_mons = len(mons)
num_rgn = len(lat_rgn)

# xr.ds variable name
data_var = '{}M_{}m'.format(var_wnt[0], abs(z_wnt))

for iii in np.arange(0, num_rgn):
    lat_rgn1 = lat_rgn[iii]
    lat_bgn = lat_rgn1[0]
    lat_end = lat_rgn1[1]

    # grid for the different regions
    xg1 = xg1_list[iii]
    yg1 = yg1_list[iii]

    # lat and lon labels
    xg2 = xg2_list[iii]
    yg2 = yg2_list[iii]

    # open file as xr.ds
    fn_in = '{}{}_{}_lat_{}_{}_xdis_{}km_thresh_xdis_{}km.nc'.format(dir_in, file_pre, var_roms, lat_bgn, lat_end, xdis, xdis_thrsh)
    ds1 = xr.open_dataset(fn_in)

    # get the 2M temperature data, put into xr.da
    da1M = ds1[var_roms]

    # get the theshold level matrix
    in_thrsh_wnt = np.where(ds1.threshold == threshold_wnt)[0]
    da1M_lvl = np.squeeze(ds1[var_level][:, :, :, in_thrsh_wnt])

    # get the seasonal theshold values as a monthly vector
    da1M_ssn_thrsh = np.squeeze(ds1['thresholdM_vec'][:, in_thrsh_wnt])

    # get the hci
    da1M_hci = np.squeeze(ds1['hci'][:, in_thrsh_wnt])

    # years
    yy = da1M.time.dt.year
    yrs = np.unique(yy)
    num_yrs = len(yrs)

    # mean and std for HCI
    hci_mn = da1M_hci.mean().data
    hci_sd = da1M_hci.std().data
    hci_p1sd = hci_mn + hci_sd
    hci_m1sd = hci_mn - hci_sd

    # lon and lat matrix
    lat = ds1[dim_wnt[0]].data
    lon = ds1[dim_wnt[1]].data
    ny = len(lat)
    nx = len(lon)

    # create new figure and close old one
    params = {'figure.figsize': fg_sz[iii]}
    rcParams.update(params)

    # subplots, rows = water mass + delt_p + error, columns = quarters
    gs1 = gridspec.GridSpec(num_rows, num_clmn)
    gs1.update(left=0.05, right=0.95, bottom=0.04, top=0.7,
               wspace=0.1, hspace=0.1)

    # get coastline
    land1 = cfeature.NaturalEarthFeature(
        'physical', 'land', '50m', edgecolor='face', facecolor='lightgray')

    for i in range(num_mons):
        # get the month wanted for the temp, level matrices, seasonal threshold, hci
        da1Mi = da1M.sel(time=da1M.time.dt.month.isin(mons[i]))
        nti = len(da1Mi.time)

        da1M_lvli = da1M_lvl.sel(time=da1M.time.dt.month.isin(mons[i]))
        da1M_ssn_thrshi = da1M_ssn_thrsh.sel(time=da1M.time.dt.month.isin(mons[i]))
        da1M_hcii = da1M_hci.sel(time=da1M.time.dt.month.isin(mons[i]))

        # create new figure and close old one
        plt.close()
        fig = plt.figure()

        # change the d_max but keep d_min
        # d_max = np.ceil(da1Mi.max())
        # nlvl1 = np.arange(d_min, d_max+dx, dx)

        for j in range(nti):
            print('i{}:{}, j{}:{}'.format(i, num_mons, j, nti))

            # get the temp, level, seasonal threshold for the given year
            data1 = np.squeeze(da1Mi.data[j, :, :])
            lvl1 = np.squeeze(da1M_lvli.data[j, :, :])
            lvl_wnt = da1M_ssn_thrshi.data[j]
            hci1 = da1M_hcii.data[j]

            # plus/minus 1 sd labels
            sd_lnd_lbl = lbl_sd[1]
            sd_lnd_clr = clr_sd[1]
            if hci1 < hci_mn:
                sd_lnd_lbl = lbl_sd[2]
                sd_lnd_clr = clr_sd[2]
            if hci1 > hci_p1sd:
                sd_lnd_lbl = lbl_sd[0]
                sd_lnd_clr = clr_sd[0]

            # contour
            ax = plt.subplot(gs1[j], projection=ccrs.PlateCarree(), frameon=False)
            if len(data1) == 0:
                data1 = np.zeros([ny, nx])
                plt.contourf(lon, lat, data1, cmap='RdYlBu_r')
            else:
                plt.contourf(lon, lat, data1, nlvl1, cmap='RdYlBu_r',
                             vmin=d_min, vmax=d_max)

            # contour lvl_wnt contour
            cs = plt.contour(lon, lat, data1,
                             np.arange(lvl_wnt, lvl_wnt+0.01, 0.01),
                             colors='k', linewidths=1, zorder=100)

            # lon and lat grid
            long, latg = np.meshgrid(lon, lat)

            # get locations below threshold
            in_lvl_thrsh = np.where(lvl1 == 1)
            xp = long[in_lvl_thrsh]
            yp = latg[in_lvl_thrsh]

            # plot points that are below lvl_wnt
            # plt.plot(xp, yp, ',', zorder=50)
            plt.plot(xp, yp, 'o', zorder=50, markersize=0.6, markerfacecolor=sd_lnd_clr,color='red',markeredgewidth=0)

            # Text on land: year
            plt.text(yr_lbl[iii][0], yr_lbl[iii][1], yrs[j], fontsize=12, ha='right', bbox=dict(boxstyle='square,pad=0.0', fc='lightgray', ec='none'))

            # Text on land: number of points
            num1 = hci1

            plt.text(pnts_lbl[iii][0], pnts_lbl[iii][1],
                     'HCI \n{:3.2f}'.format(num1), fontsize=7, ha='center',
                     bbox=dict(boxstyle='square,pad=0.1', fc=sd_lnd_clr, ec='none'))

            # coastlines
            ax.add_feature(land1)

            # grid lines
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                              linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
            gl.xlocator = mticker.FixedLocator(xg1)
            gl.ylocator = mticker.FixedLocator(yg1)

            # tick labels
            ax.set_xticks(xg2, crs=ccrs.PlateCarree())
            ax.set_yticks(yg2, crs=ccrs.PlateCarree())

            ax.set_xticklabels(xg2, color='black')
            ax.set_yticklabels(yg2, color='black')
            # if j < (num_rows-1)*num_clmn:
            if j < (nti-num_clmn):
                ax.set_xticklabels(xg2, color='None')

            if np.mod(j, num_clmn) > 0:
                ax.set_yticklabels(yg2, color='None')

            # fancy grid labels
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)

            # --remove grid box
            # ax.outline_patch.set_edgecolor('white')

            # xlimit
            plt.xlim(xlm_list[iii])

            # ylimit (extend the northern limit a little to get a grid line)
            plt.ylim([np.min(lat), np.max(lat)+0.05])

            # one Month Name title
            num_clmn_half = int(num_clmn/2) - 1
            if j == num_clmn_half:
                xlm = plt.xlim()
                ylm = plt.ylim()
                ttl_thrsh = '{}, Area Below {}$\degree$C Threshold'.format(clndr.month_name[i+1], lvl_wnt)
                plt.text(xlm[1], ylm[1], ttl_thrsh, ha='center', va='bottom', fontsize=14)

        # boundaries for colorbar, this will make discreate color intervals
        boundaries1 = np.arange(d_min, d_max+dx, dx)

        # colorbar
        ax_pos = ax.get_position()
        x_cb = ax_pos.x0 + ax_pos.width + (ax_pos.width)/5
        y_cb = ax_pos.y0
        y_hght = ax_pos.height
        cbaxes = plt.gcf().add_axes([x_cb, y_cb-(y_hght*0.1), 0.005, y_hght*0.8])

        m = plt.cm.ScalarMappable(cmap='RdYlBu_r')
        m.set_array(data1)
        m.set_clim(d_min, d_max)
        cbar = plt.colorbar(m, cax=cbaxes, label=clrbr_lbl, boundaries=boundaries1)

        # set tick on colorbar
        cbar.set_ticks(cbr_tcks)

        # legend for HCI SD
        ax_pos1 = ax.get_position()
        x_cb = ax_pos1.x0 + ax_pos1.width*2
        y_cb = ax_pos1.y0
        y_hght = ax_pos1.height
        cbaxes = plt.gcf().add_axes([x_cb, y_cb-(y_hght*0.1), ax_pos1.width*4, y_hght*0.8])
        cbaxes.axis('off')

        sd_ttl = 'HCI colors based on standard deviation (SD) and\nmean (MN) of all values over {} to {}'.format(da1M_hci.time.data[0].astype('datetime64[M]'), da1M_hci.time.data[-1].astype('datetime64[M]'))
        plt.text(0, 0.9, sd_ttl, fontsize=11)

        plt.text(0.05, 0.6, '{} (HCI > 1SD)'.format(lbl_sd[0]),
                 bbox=dict(boxstyle='round,pad={}'.format(lgnd_pad[iii]), fc=clr_sd[0], ec='none'))
        plt.text(0.05, 0.3, '{} (HCI > MN and HCI < 1SD)'.format(lbl_sd[1]),
                 bbox=dict(boxstyle='round,pad={}'.format(lgnd_pad[iii]), fc=clr_sd[1], ec='none'))
        plt.text(0.05, 0.0, '{} (HCI < MN)'.format(lbl_sd[2]),
                 bbox=dict(boxstyle='round,pad={}'.format(lgnd_pad[iii]), fc=clr_sd[2], ec='none'))

        # save figure
        # create plot output directory
        dir_plots = '{}/all_years/{:02d}/'.format(dir_plot_out, mons[i])

        # check if directory exist, if it doesn't then create
        try:
            os.makedirs(dir_plots)
        except OSError:
            if not os.path.isdir(dir_plots):
                raise
        fn_fig = '{}roms_temp{}m_rgn{}_{:02d}.png'.format(dir_plots, np.abs(z_wnt), rgn_wnt[iii], mons[i])
        plt.savefig(fn_fig)

import os
import calendar
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import calendar as clndr
from matplotlib import rcParams
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


# plot paramaters
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'RdYlBu_r',
    'axes.grid': False,
    'savefig.dpi': 300,  # to adjust notebook inline plot size
    'xtick.top':        False,  # shold the top and bottom have tick marks
    'xtick.bottom':     True,
    'xtick.major.size': 0,
    'ytick.major.size': 0,
    'ytick.direction': 'out',
    'xtick.direction': 'out',
    'axes.labelsize': 10,  # fontsize for x and y labels 
    'axes.titlesize': 10,
    'font.size': 10,  # was 10
    'legend.fontsize': 10,  # was 10
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'figure.figsize': [8.0, 8.5],
    'font.family': 'STIXGeneral',
    'mathtext.fontset': 'stix',
    'toolbar': 'None',
    'savefig.bbox': 'tight'
}
rcParams.update(params)

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# end year

# regions
rgn = np.arange(1, 5)

# --input directory
dir_out = './data_gha/HCI_ROMS/'
# dir_out = './data_x13/HCI_ROMS/'
dir_in = dir_out

# depth level
z_wnt = -2

# variable wanted
var_wnt = ['temp']

# dimensions
dim_wnt = ['latitude', 'longitude']

# File1 (ROM2 temp): x distance and lat range list
var_roms = 'temp_mtrx'
var_level = 'lvlM_mtrx'    # 1=below threshold, 0=above threshold
threshold_wnt = 'mn_75kmM'  # M=monthly threshold, S=seasonal threshold
xdis = 150
lat_rgn = [[43.5, 48], [40, 43.5], [35.5, 40], [30, 35.5]]
rgn_wnt = [1, 2, 3, 4]

# File2 (seasonal thresholds): x distance threshold
xdis_thrsh = 75

# contour levels
d_min = 7
d_max = 22
dx = 1
nlvl1 = np.arange(d_min-4*dx, d_max+8*dx, dx)

cbr_tcks = np.arange(d_min, d_max+3, 3)

# lvl_wnt dimension name of the xr.ds
var_dim = 'temp cutoff'

# time series name
var_vec = 'num_below_cutoff'

# map grid
xlm_list = [[-126.5, -122.5],
            [-126.5, -123.0],
            [-126.5, -120.5],
            [-123.5, -115.3]]
xg1_list = [np.arange(-127, -120, 1),
            np.arange(-127.0, -119.0, 1),
            np.arange(-127.0, -119.0, 1),
            np.arange(-126.0, -114.0, 1)]
yg1_list = [np.arange(43, 49, 1),
            np.arange(39.5, 44, 1),
            np.arange(35, 42, 1),
            np.arange(30, 37, 1)]

xg2_list = [xg1_list[0][0::2], xg1_list[1][1::2], xg1_list[2][2::3], xg1_list[3][1::3]]
yg2_list = [yg1_list[0][1::2], yg1_list[1][1::2], yg1_list[2][1::2], yg1_list[3][1::2]]

# colorbar label
clrbr_lbl = 'Temp ($\degree$C)'

# rows and columns of figure
num_rows = 6
num_clmn = 8

# months
mons = np.arange(1, 13)

# lat and lon position of text labels 
yr_lbl = [[-122.55, 47.3],[-123.1, 42.7],[-120.9, 39.2],[-115.4, 34.4]]
pnts_lbl = [[-123.15, 46.0], [-123.6, 41.7], [-121.4, 37.6], [-116.3, 32.6]]
lgnd_pad = [0.35, 0.35, 0.15, 0.25]

# figure size
fg_sz = [[8.0, 10], [9.75, 9.75], [9.25, 8.5], [11.25, 9.0]]

# annotate hci based on sd
lbl_sd = ['Low Compression', 'Medium Compression', 'High Compression']
clr_sd = ['springgreen', 'yellow', 'tomato']


# --IEA file names
file_pre = 'hci'

# --plot directory
dir_plot_out = './figures_gha/HCI_ROMS/'
# dir_plot_out = './figures_x13/HCI_ROMS/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# len of input variables
num_data = len(rgn)

# size of input variables
num_mons = len(mons)
num_rgn = len(lat_rgn)

# xr.ds variable name
data_var = '{}M_{}m'.format(var_wnt[0], abs(z_wnt))

for iii in np.arange(0, num_rgn):
    lat_rgn1 = lat_rgn[iii]
    lat_bgn = lat_rgn1[0]
    lat_end = lat_rgn1[1]

    # grid for the different regions
    xg1 = xg1_list[iii]
    yg1 = yg1_list[iii]

    # lat and lon labels
    xg2 = xg2_list[iii]
    yg2 = yg2_list[iii]

    # open file as xr.ds
    fn_in = '{}{}_{}_lat_{}_{}_xdis_{}km_thresh_xdis_{}km.nc'.format(dir_in, file_pre, var_roms, lat_bgn, lat_end, xdis, xdis_thrsh)
    ds1 = xr.open_dataset(fn_in)

    # get the 2M temperature data, put into xr.da
    da1M = ds1[var_roms]

    # get the theshold level matrix
    in_thrsh_wnt = np.where(ds1.threshold == threshold_wnt)[0]
    da1M_lvl = np.squeeze(ds1[var_level][:, :, :, in_thrsh_wnt])

    # get the seasonal theshold values as a monthly vector
    da1M_ssn_thrsh = np.squeeze(ds1['thresholdM_vec'][:, in_thrsh_wnt])

    # get the hci
    da1M_hci = np.squeeze(ds1['hci'][:, in_thrsh_wnt])

    # years
    yy = da1M.time.dt.year
    yrs = np.unique(yy)
    num_yrs = len(yrs)

    # mean and std for HCI
    hci_mn = da1M_hci.mean().data
    hci_sd = da1M_hci.std().data
    hci_p1sd = hci_mn + hci_sd
    hci_m1sd = hci_mn - hci_sd

    # lon and lat matrix
    lat = ds1[dim_wnt[0]].data
    lon = ds1[dim_wnt[1]].data
    ny = len(lat)
    nx = len(lon)

    # create new figure and close old one
    params = {'figure.figsize': fg_sz[iii]}
    rcParams.update(params)

    # subplots, rows = water mass + delt_p + error, columns = quarters
    gs1 = gridspec.GridSpec(num_rows, num_clmn)
    gs1.update(left=0.05, right=0.95, bottom=0.04, top=0.7,
               wspace=0.1, hspace=0.1)

    # get coastline
    land1 = cfeature.NaturalEarthFeature(
        'physical', 'land', '50m', edgecolor='face', facecolor='lightgray')

    for i in range(num_mons):
        # get the month wanted for the temp, level matrices, seasonal threshold, hci
        da1Mi = da1M.sel(time=da1M.time.dt.month.isin(mons[i]))
        nti = len(da1Mi.time)

        da1M_lvli = da1M_lvl.sel(time=da1M.time.dt.month.isin(mons[i]))
        da1M_ssn_thrshi = da1M_ssn_thrsh.sel(time=da1M.time.dt.month.isin(mons[i]))
        da1M_hcii = da1M_hci.sel(time=da1M.time.dt.month.isin(mons[i]))

        # create new figure and close old one
        plt.close()
        fig = plt.figure()

        # change the d_max but keep d_min
        # d_max = np.ceil(da1Mi.max())
        # nlvl1 = np.arange(d_min, d_max+dx, dx)

        for j in range(nti):
            print('i{}:{}, j{}:{}'.format(i, num_mons, j, nti))

            # get the temp, level, seasonal threshold for the given year
            data1 = np.squeeze(da1Mi.data[j, :, :])
            lvl1 = np.squeeze(da1M_lvli.data[j, :, :])
            lvl_wnt = da1M_ssn_thrshi.data[j]
            hci1 = da1M_hcii.data[j]

            # plus/minus 1 sd labels
            sd_lnd_lbl = lbl_sd[1]
            sd_lnd_clr = clr_sd[1]
            if hci1 < hci_mn:
                sd_lnd_lbl = lbl_sd[2]
                sd_lnd_clr = clr_sd[2]
            if hci1 > hci_p1sd:
                sd_lnd_lbl = lbl_sd[0]
                sd_lnd_clr = clr_sd[0]

            # contour
            ax = plt.subplot(gs1[j], projection=ccrs.PlateCarree(), frameon=False)
            if len(data1) == 0:
                data1 = np.zeros([ny, nx])
                plt.contourf(lon, lat, data1, cmap='RdYlBu_r')
            else:
                plt.contourf(lon, lat, data1, nlvl1, cmap='RdYlBu_r',
                             vmin=d_min, vmax=d_max)

            # contour lvl_wnt contour
            cs = plt.contour(lon, lat, data1,
                             np.arange(lvl_wnt, lvl_wnt+0.01, 0.01),
                             colors='k', linewidths=1, zorder=100)

            # lon and lat grid
            long, latg = np.meshgrid(lon, lat)

            # get locations below threshold
            in_lvl_thrsh = np.where(lvl1 == 1)
            xp = long[in_lvl_thrsh]
            yp = latg[in_lvl_thrsh]

            # plot points that are below lvl_wnt
            # plt.plot(xp, yp, ',', zorder=50)
            plt.plot(xp, yp, 'o', zorder=50, markersize=0.6, markerfacecolor=sd_lnd_clr,color='red',markeredgewidth=0)

            # Text on land: year
            plt.text(yr_lbl[iii][0], yr_lbl[iii][1], yrs[j], fontsize=12, ha='right', bbox=dict(boxstyle='square,pad=0.0', fc='lightgray', ec='none'))

            # Text on land: number of points
            num1 = hci1

            plt.text(pnts_lbl[iii][0], pnts_lbl[iii][1],
                     'HCI \n{:3.2f}'.format(num1), fontsize=7, ha='center',
                     bbox=dict(boxstyle='circle,pad=0.3', fc=sd_lnd_clr, ec='none'))

            # coastlines
            ax.add_feature(land1)

            # grid lines
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                              linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
            gl.xlocator = mticker.FixedLocator(xg1)
            gl.ylocator = mticker.FixedLocator(yg1)

            # tick labels
            ax.set_xticks(xg2, crs=ccrs.PlateCarree())
            ax.set_yticks(yg2, crs=ccrs.PlateCarree())

            ax.set_xticklabels(xg2, color='black')
            ax.set_yticklabels(yg2, color='black')
            # if j < (num_rows-1)*num_clmn:
            if j < (nti-num_clmn):
                ax.set_xticklabels(xg2, color='None')

            if np.mod(j, num_clmn) > 0:
                ax.set_yticklabels(yg2, color='None')

            # fancy grid labels
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)

            # --remove grid box
            # ax.outline_patch.set_edgecolor('white')

            # xlimit
            plt.xlim(xlm_list[iii])

            # ylimit (extend the northern limit a little to get a grid line)
            plt.ylim([np.min(lat), np.max(lat)+0.05])

            # one Month Name title
            num_clmn_half = int(num_clmn/2) - 1
            if j == num_clmn_half:
                xlm = plt.xlim()
                ylm = plt.ylim()
                ttl_thrsh = '{}, Area Below {}$\degree$C Threshold'.format(clndr.month_name[i+1], lvl_wnt)
                plt.text(xlm[1], ylm[1], ttl_thrsh, ha='center', va='bottom', fontsize=14)

        # boundaries for colorbar, this will make discreate color intervals
        boundaries1 = np.arange(d_min, d_max+dx, dx)

        # colorbar
        ax_pos = ax.get_position()
        x_cb = ax_pos.x0 + ax_pos.width + (ax_pos.width)/5
        y_cb = ax_pos.y0
        y_hght = ax_pos.height
        cbaxes = plt.gcf().add_axes([x_cb, y_cb-(y_hght*0.1), 0.005, y_hght*0.8])

        m = plt.cm.ScalarMappable(cmap='RdYlBu_r')
        m.set_array(data1)
        m.set_clim(d_min, d_max)
        cbar = plt.colorbar(m, cax=cbaxes, label=clrbr_lbl, boundaries=boundaries1)

        # set tick on colorbar
        cbar.set_ticks(cbr_tcks)

        # legend for HCI SD
        ax_pos1 = ax.get_position()
        x_cb = ax_pos1.x0 + ax_pos1.width*2
        y_cb = ax_pos1.y0
        y_hght = ax_pos1.height
        cbaxes = plt.gcf().add_axes([x_cb, y_cb-(y_hght*0.1), ax_pos1.width*4, y_hght*0.8])
        cbaxes.axis('off')

        sd_ttl = 'HCI colors based on standard deviation (SD) and\nmean (MN) of all values over {} to {}'.format(da1M_hci.time.data[0].astype('datetime64[M]'), da1M_hci.time.data[-1].astype('datetime64[M]'))
        plt.text(0, 0.9, sd_ttl, fontsize=11)

        plt.text(0.05, 0.6, '{} (HCI > 1SD)'.format(lbl_sd[0]),
                 bbox=dict(boxstyle='round,pad={}'.format(lgnd_pad[iii]), fc=clr_sd[0], ec='none'))
        plt.text(0.05, 0.3, '{} (HCI > MN and HCI < 1SD)'.format(lbl_sd[1]),
                 bbox=dict(boxstyle='round,pad={}'.format(lgnd_pad[iii]), fc=clr_sd[1], ec='none'))
        plt.text(0.05, 0.0, '{} (HCI < MN)'.format(lbl_sd[2]),
                 bbox=dict(boxstyle='round,pad={}'.format(lgnd_pad[iii]), fc=clr_sd[2], ec='none'))

        # save figure
        # create plot output directory
        dir_plots = '{}/all_years/{:02d}/'.format(dir_plot_out, mons[i])

        # check if directory exist, if it doesn't then create
        try:
            os.makedirs(dir_plots)
        except OSError:
            if not os.path.isdir(dir_plots):
                raise
        fn_fig = '{}roms_temp{}m_rgn{}_{:02d}.png'.format(dir_plots, np.abs(z_wnt), rgn_wnt[iii], mons[i])
        plt.savefig(fn_fig)

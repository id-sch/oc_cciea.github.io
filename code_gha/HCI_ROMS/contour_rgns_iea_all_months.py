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

# --IEA file names
file_pre = 'hci'

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
            [-126.5, -121],
            [-123.5, -115.3]]
xg1_list = [np.arange(-127, -120, 1),
            np.arange(-127.0, -119.0, 1),
            np.arange(-127.0, -119.0, 1),
            np.arange(-126.0, -114.0, 1)]
yg1_list = [np.arange(43, 49, 1),
            np.arange(39.5, 44, 1),
            np.arange(35, 43, 1),
            np.arange(30, 37, 1)]

xg2_list = [xg1_list[0][0::2], xg1_list[1][1::2], xg1_list[2][2::3], xg1_list[3][1::3]]
yg2_list = [yg1_list[0][1::2], yg1_list[1][1::2], yg1_list[2][1::2], yg1_list[3][1::2]]

# colorbar label
clrbr_lbl = 'Temp ($\degree$C)'

# rows and columns of figure
num_rows = 1
num_clmn = 1

# months
mons = np.arange(1, 13)

# lat and lon position of text labels 
yr_lbl = [[-122.55, 47.2],[-123.1, 42.7],[-121.05, 39.2],[-115.4, 34.3]] 
pnts_lbl = [[-123.0, 46.35], [-123.4, 42.0], [-121.5, 38.35], [-116.0, 33.1]]
lgnd_pad = [0.35, 0.35, 0.15, 0.25]

# figure size
fg_sz = [[4.25, 5.5], [4.25, 5.5], [4.25, 5.5], [4.25, 5.5]]

# annotate hci based on sd
lbl_sd = ['High Compression', 'Medium Compression', 'Low Compression']
clr_sd = ['tomato', 'yellow', 'springgreen']

lbl_sd = ['Low Compression', 'Medium Compression', 'High Compression']
clr_sd = ['springgreen', 'yellow', 'tomato']

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

# subplots, rows = water mass + delt_p + error, columns = quarters
gs1 = gridspec.GridSpec(num_rows, num_clmn)
gs1.update(left=0.05, right=0.95, bottom=0.04, top=0.7,
           wspace=0.1, hspace=0.1)

# get coastline
land1 = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                     edgecolor='face', facecolor='lightgray')


for iii in np.arange(0, num_rgn):
# for iii in np.arange(2, 3):
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

    # lon and lat matrix
    lat = ds1[dim_wnt[0]].data
    lon = ds1[dim_wnt[1]].data
    ny = len(lat)
    nx = len(lon)

    # lon and lat grid
    long, latg = np.meshgrid(lon, lat)

    for i in range(num_mons):
        da1Mi = da1M.sel(time=da1M.time.dt.month.isin(mons[i]))
        da1M_lvli = da1M_lvl.sel(time=da1M_lvl.time.dt.month.isin(mons[i]))
        da1M_ssn_thrshi = da1M_ssn_thrsh.sel(time=da1M_ssn_thrsh.time.dt.month.isin(mons[i]))
        da1M_hcii = da1M_hci.sel(time=da1M_hci.time.dt.month.isin(mons[i]))

        # mean and std for HCI
        hci_mn = da1M_hci.mean().data
        hci_sd = da1M_hci.std().data
        hci_p1sd = hci_mn + hci_sd
        hci_m1sd = hci_mn - hci_sd

        # years
        yrs = da1Mi.time.dt.year.data
        num_yrs = len(yrs)

        # change figure size
        params = {'figure.figsize': fg_sz[iii]}
        rcParams.update(params)

        for j in range(num_yrs):
            # create new figure and close old one
            plt.close()
            fig = plt.figure()

            # dataarray with year wanted
            da1Yi = da1Mi.sel(time=da1Mi.time.dt.year.isin(yrs[j]))
            da1Y_lvli = da1M_lvli.sel(time=da1M_lvli.time.dt.year.isin(yrs[j]))
            da1Y_ssn_thrshi = da1M_ssn_thrshi.sel(time=da1M_ssn_thrshi.time.dt.year.isin(yrs[j]))
            da1Y_hcii = da1M_hcii.sel(time=da1M_hcii.time.dt.year.isin(yrs[j]))

            # get the temp, level, seasonal threshold for the given year
            data1 = np.squeeze(da1Yi)
            yy1 = da1Yi.time.dt.year.data[0]
            mm1 = da1Yi.time.dt.month.data[0]

            lvl1 = np.squeeze(da1Y_lvli.data)
            lvl_wnt = da1Y_ssn_thrshi.data
            hci1 = da1M_hcii.data[j]

            # get locations below threshold
            in_lvl_thrsh = np.where(lvl1 > 0)
            xp = long[in_lvl_thrsh]
            yp = latg[in_lvl_thrsh]

            # plus/minus 1 sd labels
            sd_lnd_lbl = lbl_sd[1]
            sd_lnd_clr = clr_sd[1]
            sd_lgnd_txt = '{} (HCI > MN and HCI < 1SD)'.format(lbl_sd[1])
            if hci1 < hci_mn:
                sd_lnd_lbl = lbl_sd[2]
                sd_lnd_clr = clr_sd[2]
                sd_lgnd_txt = '{} (HCI < MN)'.format(lbl_sd[2])
            if hci1 > hci_p1sd:
                sd_lnd_lbl = lbl_sd[0]
                sd_lnd_clr = clr_sd[0]
                sd_lgnd_txt = '{} (HCI > 1SD)'.format(lbl_sd[0])

            # contour
            ax = plt.subplot(gs1[0], projection=ccrs.PlateCarree(), frameon=False)
            if len(data1)==0:
                data1 = np.zeros([ny, nx])
                plt.contourf(lon, lat, data1, cmap='RdYlBu_r')
            else:
                plt.contourf(lon, lat, data1, nlvl1, cmap='RdYlBu_r',
                             vmin=d_min, vmax=d_max)

            # plot points that are below lvl_wnt
            plt.plot(xp, yp, 'o', color=sd_lnd_clr, markersize=0.6, zorder=50)

            # contour lvl_wnt contour
            cs = plt.contour(lon, lat, data1,
                             np.arange(lvl_wnt[0], lvl_wnt[0]+0.01, 0.01),
                             colors='k', linewidths=1, zorder=100)

            # clabel
            labels1 = plt.clabel(cs, fmt='%4.1f')
            for l in labels1:
                l.set_rotation(0)
                l.set_color('black')

            # Text on land: month and year
            mon_yr_txt = '{}\n{}'.format(clndr.month_name[mons[i]], yrs[j])
            plt.text(yr_lbl[iii][0], yr_lbl[iii][1], mon_yr_txt, fontsize=18, ha='right', bbox=dict(boxstyle='square,pad=0.0', fc='lightgray', ec='none'))

            # Text on land: number of points
            num1 = hci1

            plt.text(pnts_lbl[iii][0], pnts_lbl[iii][1],
                     'HCI \n{:3.2f}'.format(num1), fontsize=11, ha='center',
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

            # fancy grid labels
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)

            # xlimit
            plt.xlim(xlm_list[iii])

            # ylimit (extend the northern limit a little to get a grid line)
            plt.ylim([np.min(lat), np.max(lat)+0.05])

            # boundaries for colorbar, this will make discreate color intervals
            boundaries1 = np.arange(d_min, d_max+dx, dx)

            # colorbar
            ax_pos = ax.get_position()
            x_cb = ax_pos.x0 + ax_pos.width + (ax_pos.width)/15
            y_cb = ax_pos.y0
            y_hght = ax_pos.height
            cbaxes = plt.gcf().add_axes([x_cb, y_cb+(y_hght*0.1)/2.0, 0.025, y_hght*0.9])

            m = plt.cm.ScalarMappable(cmap='RdYlBu_r')
            m.set_array(data1)
            m.set_clim(d_min, d_max)
            cbar = plt.colorbar(m, cax=cbaxes, label=clrbr_lbl, format='%4.2f', boundaries=boundaries1)
            # set tick on colorbar
            cbar.set_ticks(cbr_tcks)

            # SD text on bottom
            x_cb = ax_pos.x0
            y_cb = ax_pos.y0
            cbaxes1 = plt.gcf().add_axes([x_cb, y_cb-0.2, ax_pos.width, 0.2])
            cbaxes1.axis('off')
            sd_ttl = 'HCI color based on standard deviation (SD) and mean (MN)\nof all values over {} to {}'.format(da1M_hci.time.data[0].astype('datetime64[M]'), da1M_hci.time.data[-1].astype('datetime64[M]'))
            plt.text(0, 0.45, sd_ttl, fontsize=12)
            plt.text(0.05, 0.25, sd_lgnd_txt,
                     bbox=dict(boxstyle='round,pad={}'.format(lgnd_pad[iii]),
                               fc=sd_lnd_clr, ec='none'))

            # save figure
            # create plot output directory
            dir_plots = '{}/all_months/{}/{:02d}/'.format(dir_plot_out, yrs[j], mons[i])

            # check if directory exist, if it doesn't then create
            try:
                os.makedirs(dir_plots)
            except OSError:
                if not os.path.isdir(dir_plots):
                    raise
            fn_fig = '{}roms_temp{}m_rgn{}_{}_{:02d}.png'.format(dir_plots, np.abs(z_wnt), rgn_wnt[iii], yrs[j], mons[i])
            plt.savefig(fn_fig)

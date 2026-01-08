import os
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
from C_iea import C_iea


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2025
# iea_yr = 2025


# dir of spatial IEA stats
dir_out = './data_gha/SSTspatial/'
# dir_out = './oc_cciea.github.io/data_gha/SSTspatial/'
dir_in = dir_out

# variables for fun_xr_ds2IEA_contour3
nlvl1 = [np.arange(-3.0, 3.5, 0.5), np.arange(-3.0, 3.5, 0.5), np.arange(-3, 3.5, 0.5)]
ttl = ['SST anom (°C)', '5-yr Mean / SD', '5-yr Trend / SD']

# interpolate every 1 degree intervals
intrp_type = (1, 1)

# tick labels
ytck = np.arange(35, 65, 10)
xtck = np.arange(-160, -100, 20)

# seasons to contour
season_lbl = ['winter', 'spring', 'summer', 'fall']

# output file name
file_pre_out = 'oc_sst_spatial'

# ds variables
ds1_var = ['coord_mtrx', 'data_mtrx', 'mrkr_mtrx', 'ts_mtrx']

# example timeseries variables
freq_wnt = 'YS'
wndw = 5
yy_end = iea_yr
yr_clim_bgn = 1982
yr_clim_end = yy_end

# figure size
fig_wdth = 11
fig_hght = 8.5

# column label
clmn_lbl = ['min', 'max']

# row label
row_lbl = ['anom', 'mean5', 'trend5']
sup_ttl_lbl = ['Anom', '5-year mean', '5-year trend']

# plot 3 time series along the coast at these latitudes
lat_ts = [45.125, 39.125, 33.125]

x_ts3 = [0.72, 0.72, 0.72]
y_ts3 = [0.77, 0.55, 0.33]

lenx_ts3 = [0.2, 0.2, 0.2]
leny_ts3 = [0.15, 0.15, 0.15]

# plot 2 min/max time series
x_ts2 = [0.125, 0.425]
y_ts2 = [0.77, 0.77]

lenx_ts2 = [0.2, 0.2]
leny_ts2 = [0.15, 0.15]
lbl_ts2 = ['Min', 'Max']

# --plot directory
dir_plot_out = './figures_gha/SSTspatial/'
# dir_plot_out = './figures_x13/SSTspatial/'
fig_type = '.png'

# figure size
fig_wdth = 11
fig_hght = 8.5

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

# create plot output directory
dir_plots = '{}{}/'.format(dir_plot_out, iea_yr)

# length of input variables
num_season = len(season_lbl)
num_ds_var = len(ds1_var)
num_map = len(row_lbl)

# define fireice colormap
cmap1 = mpl.colors.LinearSegmentedColormap.from_list(
    "", ["mediumblue", "dodgerblue", "cyan", "white", "yellow", "tomato", "red"])

# subplots using gridspec
num_row = 1
num_clmn = 1
gs1 = gridspec.GridSpec(num_row, num_clmn)
gs1.update(left=0.1, right=0.9, bottom=0.1,
           top=0.9, wspace=0.1, hspace=0.1)

# coastlines
land1 = cfeature.NaturalEarthFeature(
    'physical', 'land', '50m', edgecolor='lightgray',
    facecolor=[0.9, 0.9, 0.9], zorder=100, linewidth=0.5)

# open the data and contour
fn_list = list()
for i in range(num_season):
    fn_in_ssn = '{}anom_mn5_trnd5_{}_clim_{}_{}.nc'.format(
          dir_in, season_lbl[i], yr_clim_bgn, yr_clim_end)

    ds1 = xr.open_dataset(fn_in_ssn)

    # lon and lat
    lon = ds1[ds1_var[0]].data[:, 0] - 360
    lat = ds1[ds1_var[0]].data[:, 1]

    # if intrp_type is a tuple then create grid for bilinear interp
    if isinstance(intrp_type, tuple):
        dx = intrp_type[0]
        dy = intrp_type[1]
        lon_vec = np.arange(np.nanmin(lon), np.nanmax(lon)+dx, dx)
        lat_vec = np.arange(np.nanmin(lat), np.nanmax(lat)+dy, dy)
        xv, yv = np.meshgrid(lon_vec, lat_vec)

    for j in range(num_map):
        # open new figure of set size
        plt.close()
        fig = plt.figure(figsize=(fig_wdth, fig_hght))

        # data
        data1 = ds1[ds1_var[1]].data[:, j]

        # if intrp_type is a tuple then bilinear interp
        if isinstance(intrp_type, tuple):
            ind = np.isfinite(data1)
            data1g = griddata((lon[ind], lat[ind]),
                              data1[ind], (xv, yv), method='linear')

        # markers
        marker1 = ds1[ds1_var[2]].data[:, j]

        # set contour limits and extend nlvl1 so that there are no white contour areas where values exceed the nlvl1 limits
        dmin1 = nlvl1[j][0]
        dmax1 = nlvl1[j][-1]
        dnlvl1 = np.nanmean(np.diff(nlvl1))
        nlvl1_ext = np.arange(dmin1*2, dmax1*2 + dnlvl1, dnlvl1)

        # plot stations and contour data
        ax1 = plt.subplot(gs1[0],  projection=ccrs.PlateCarree())
        cf1 = plt.contourf(xv, yv, data1g, nlvl1_ext, vmin=dmin1, vmax=dmax1, cmap=cmap1, linestyles='-')

        plt.xlim([-164, -105])

        # remove box around figure
        #ax1.background_patch.set_visible(False)   # Background
        ax1.spines['geo'].set_edgecolor('white')
        # land
        ax1.add_feature(land1)

        # index of less/greater than 1 sd
        in1 = np.where(marker1 == '.')
        plt.plot(lon[in1], lat[in1], marker='.',
                 linestyle='none', color='grey', markersize=2)


        # index of smallest/greastest of time series
        in1 = np.where(marker1 == '+')
        plt.plot(lon[in1], lat[in1], marker='x',
                 linestyle='none', color='black', markersize=4,
                 markeredgewidth=0.3)

        # set x and y limits
        xlm = plt.xlim()
        ylm = plt.ylim()

        # remove puget sound
        in_lon_chck = np.where(lon > -127.5)[0]
        in_lat_chck = np.where(lat[in_lon_chck] > 46.5)[0]
        in_ps = in_lon_chck[in_lat_chck]

        data_chck = np.copy(data1)
        data_chck[in_ps] = np.nan

        in_lat52 = np.where(lat > 52)[0]
        data_chck[in_lat52] = np.nan

        in_keep = np.isfinite(data_chck).nonzero()[0]

        # 2 time min/max time series at top -- grid location and text
        # in_keep = np.where(lat < 52)
        data_keep = data1[in_keep]
        lat_keep = lat[in_keep]
        lon_keep = lon[in_keep]

        for iii in range(len(lbl_ts2)):
            if iii == 0:
                in_min_max = np.nanargmin(data_keep)
            else:
                in_min_max = np.nanargmax(data_keep)
            data_min_max = data_keep[in_min_max]
            lat_min_max = lat_keep[in_min_max]
            lon_min_max = lon_keep[in_min_max]

            plt.plot(lon_min_max, lat_min_max, '*', color='green', zorder=500)
            plt.text(lon_min_max+0.2, lat_min_max-0.2, lbl_ts2[iii], fontsize=18, zorder=1000, color='green', va='center')

        # 3 time series along the coast -- grid location and text
        for iii in range(len(lat_ts)):
            in_lat_ts = np.where(lat_ts[iii] == lat)[0]
            data1i = data1[in_lat_ts]
            ind1i = np.isfinite(data1i).nonzero()[0]
            in_ts = in_lat_ts[ind1i[-1]]
            lat_edge = lat[in_ts]
            lon_edge = lon[in_ts]
            ts = ds1['ts_mtrx'][in_ts, :]
            yrs = ds1['time'].dt.year.data

            plt.plot(lon_edge, lat_edge, '*', color='green', zorder=500)
            plt.text(lon_edge+0.2, lat_edge-0.2, str(iii+1), fontsize=18, zorder=1000, color='green', va='center')

        plt.xlim(xlm)
        plt.ylim(ylm)

        # 2 time min/max time series at top -- IEA time series
        # in_keep = np.where(lat < 50)[0]
        data_keep = data1[in_keep]
        lat_keep = lat[in_keep]
        lon_keep = lon[in_keep]
        ts_keep = ds1['ts_mtrx'][in_keep, :]
        coords_keep = ds1[ds1_var[0]][in_keep, :]

        for iii in range(len(lbl_ts2)):
            if iii == 0:
                in_min_max = np.nanargmin(data_keep)
            else:
                in_min_max = np.nanargmax(data_keep)
            data_min_max = data_keep[in_min_max]
            lat_min_max = lat_keep[in_min_max]
            lon_min_max = lon_keep[in_min_max]
            ts = ts_keep[in_min_max, :]
            coords1 = coords_keep[in_min_max, :]

            ax2 = fig.add_axes([x_ts2[iii], y_ts2[iii], lenx_ts2[iii], leny_ts2[iii]])
            ps1 = pd.Series(data=ts.data, index=ts.time.data).asfreq(freq_wnt)

            z1 = C_iea(ps1, yr_clim_bgn=yr_clim_bgn, yr_clim_end=yr_clim_end,
                       wndw=wndw, yy_end=yy_end, lon=coords1.data[0],
                       lat=coords1.data[1], fs_mrkr=14, fs_tcks=12)
            z1.plot()
            plt.title('')
            dy_ts = (np.max(z1.data) - np.min(z1.data))/10
            plt.ylim([np.min(z1.data) - dy_ts, np.max(z1.data) + dy_ts])
            num_yrs = len(z1.data.values)
            plt.title('Grid Location {}'.format(lbl_ts2[iii]), pad=0)

        # 3 time series along the coast -- IEA time series
        for iii in range(len(lat_ts)):
            in_lat_ts = np.where(lat_ts[iii] == lat)[0]
            data1i = data1[in_lat_ts]
            ind1i = np.isfinite(data1i).nonzero()[0]
            in_ts = in_lat_ts[ind1i[-1]]
            lat_edge = lat[in_ts]
            lon_edge = lon[in_ts]
            ts = ds1['ts_mtrx'][in_ts, :]
            yrs = ds1['time'].dt.year.data

            ax2 = fig.add_axes([x_ts3[iii], y_ts3[iii], lenx_ts3[iii], leny_ts3[iii]])
            ps1 = pd.Series(data=ts.data, index=ts.time.data).asfreq(freq_wnt)
            coords1 = ds1[ds1_var[0]][in_ts, :]
            z1 = C_iea(ps1, yr_clim_bgn=yr_clim_bgn, yr_clim_end=yr_clim_end,
                       wndw=wndw, yy_end=yy_end, lon=coords1.data[0],
                       lat=coords1.data[1], fs_mrkr=14, fs_tcks=12)
            z1.plot()
            plt.title('')
            dy_ts = (np.max(z1.data) - np.min(z1.data))/8
            plt.ylim([np.min(z1.data) - dy_ts, np.max(z1.data) + dy_ts])
            num_yrs = len(z1.data.values)
            plt.title('Grid Location {}'.format(iii+1), pad=0)

        # colorbar
        ax_pos = ax1.get_position()
        cbaxes = plt.gcf().add_axes([0.4, 0.19, 0.3, 0.015])

        m = plt.cm.ScalarMappable(cmap=cmap1)
        m.set_array(data1g)
        m.set_clim(dmin1, dmax1)

        ttl_lbl = '{} {}'.format(season_lbl[i].capitalize(), ttl[j])
        cbar = plt.colorbar(m, cax=cbaxes,
                            orientation="horizontal", shrink=0.95,
                            pad=0.04, label=ttl_lbl, aspect=30,
                            boundaries=np.arange(dmin1, dmax1+dnlvl1, dnlvl1))

        # cbaxes.xaxis.tick_bottom()
        cbar.ax.tick_params(axis='x', direction='in', labelsize=14)
        cbar.ax.xaxis.set_ticks_position('bottom')
        cbar.ax.xaxis.set_label_position('bottom')
        cbar.set_label(ttl_lbl, labelpad=4)

        # save figure
        dir_plots = '{}{}/example_ts/{}/'.format(dir_plot_out, iea_yr, season_lbl[i])

        # check if directory exist, if it doesn't then create
        try:
            os.makedirs(dir_plots)
        except OSError:
            if not os.path.isdir(dir_plots):
                raise
        fn_fig = '{}{}_{}_example_ts.png'.format(dir_plots, season_lbl[i], row_lbl[j])
        plt.savefig(fn_fig)

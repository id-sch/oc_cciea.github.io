import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import gridspec
from matplotlib import rcParams
from matplotlib import interactive
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
interactive(True)


# -----------------------------------------------------------------------------
# Input variables, change these
# -----------------------------------------------------------------------------
# .1) sst time series
iea_yr = 2026
fn1 = './data_gha/SSTspatial/TS_monthly.nc'


# .2) SST
# x and y limits
xlm = [-128, -116.8]
ylm = [32, 48]

# variable
var_wnt = ['sst']

# dimension
dim_lbl = ['lat_vec', 'lon_vec', 'time']

# .3) months
mon_list = [[1, 3], [4, 6], [7, 9], [10, 12]]

# season label
ssn_lbl = ['Winter', 'Spring', 'Summer', 'Fall']
xssn = -117.5
yssn = 46

# .4)
# figure usually has 1 columns, 2 rows
num_clmn = 1
num_row = 1

# figure size
fig_wdth = 8.5
fig_hght = 6

# nlevels
dmin = -2
dmax = 2
dl = 0.25
lvls = np.arange(dmin*10, 10*dmax+10*dl, 10*dl)/10

lvls2 = [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2]


# x and y limits
xlm = [-128, -116.8]
ylm = [32, 48]

# lat and lon labels
xg2 = np.arange(-130, -110, 4)
yg2 = np.arange(30, 54, 4)

# colorbar label
clrbr_lbl = 'SST Anom (°C)'

# .5) plot directory
dir_plot_out = './figures_gha/miniESR/'
# dir_plot_out = './figures_x13/miniESR'

# -----------------------------------------------------------------------------
# END: Input variables, change these
# -----------------------------------------------------------------------------
# input variable size
num_mon_list = len(mon_list)

ds1 = xr.open_dataset(fn1).sel(lon_vec=slice(xlm[0]+360, xlm[1]+360)).sel(lat_vec=slice(ylm[0], ylm[1]))

# get the years up to iea_yr
yrs = ds1.time.dt.year
yr_bgn = yrs[0]
yr_end = iea_yr

da1 = ds1[var_wnt[0]].sel(time=ds1.time.dt.year.isin(np.arange(yr_bgn, yr_end+1)))

# lat and lon
lat = ds1[dim_lbl[0]].data
lon = ds1[dim_lbl[1]].data


# setup subplots, spacing and figure size
gs1 = gridspec.GridSpec(num_row, num_clmn)
gs1.update(left=0.05, right=0.9, bottom=0.05,
           top=0.65, wspace=0.1, hspace=0.1)
fig = plt.figure(figsize=(fig_wdth, fig_hght))

mpl.rcParams.update({'font.size': 7})

# get coastline
land1 = cfeature.NaturalEarthFeature(
    'physical', 'land', '50m', edgecolor='slategray', facecolor='lightgray')

# loop over the seasons
for iii in range(num_mon_list):
    # create new figure and close old one
    plt.close()
    fig = plt.figure()

    # create season, anom
    mon_bgn = mon_list[iii][0]
    mon_end = mon_list[iii][1]
    if mon_bgn > mon_end:
        num_mon_avg = (12 - mon_bgn + 1) + mon_end
    else:
        num_mon_avg = mon_end - mon_bgn + 1

    da1_mn = da1.rolling(time=num_mon_avg, center=False).mean().sel(time=da1.time.dt.month.isin(mon_end)).groupby('time.year').mean('time')

    # check last year of the seasonal mean, sometimes there is not a full year of sst data
    yr_mn_end = da1_mn.year[-1].data
    if yr_mn_end >= yr_end:

        # calculate the clim
        da1_clim = da1_mn.mean('year')

        da1_anom = da1_mn.sel(year=yr_end) - da1_clim
        anom1 = da1_anom.data

        min_anom = np.nanmin(anom1)
        max_anom = np.nanmax(anom1)

        if min_anom < dmin and max_anom > dmax:
            extnd = 'both'
        if min_anom < dmin and max_anom < dmax:
            extnd = 'min'
        if min_anom > dmin and max_anom > dmax:
            extnd = 'max'
        if min_anom > dmin and max_anom < dmax:
            extnd = 'neither'

        ax = fig.add_subplot(
            gs1[0], projection=ccrs.PlateCarree())

        cs2 = plt.contour(lon, lat, anom1, levels=lvls2, transform=ccrs.PlateCarree(), colors=[0.2,0.2,0.2], linewidths=0.5)
        plt.clabel(cs2, fontsize=3)
        cs = plt.contourf(lon, lat, anom1, cmap='bwr', levels=lvls,
                          extend='both', vmax=dmax, vmin=dmin,
                          transform=ccrs.PlateCarree())

        # coastlines
        ax.add_feature(land1, linewidth=0.5, zorder=40)

        # tick labels
        ax.set_xticks(xg2, crs=ccrs.PlateCarree())
        ax.set_yticks(yg2, crs=ccrs.PlateCarree())

        # fancy grid labels
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)

        # x and y limits
        plt.xlim(xlm)
        plt.ylim(ylm)

    #    plt.text(xssn, yssn, '{} {}'.format(
    #        ssn_lbl, yrs_wnt[i]), fontsize=10, horizontalalignment='right')


        # add season and year
        plt.text(xssn, yssn, '{} {}'.format(
            ssn_lbl[iii], iea_yr), fontsize=10, horizontalalignment='right', zorder=50)

        # colorbar
        ax_pos = ax.get_position()
        x_cb = ax_pos.x0 + ax_pos.width + ax_pos.width/20.0
        y_cb = ax_pos.y0
        y_hght = ax_pos.height
        cbaxes = plt.gcf().add_axes(
            [x_cb, y_cb+(0.2*y_hght)/2.0, 0.005, y_hght*0.8])

        m = plt.cm.ScalarMappable(cmap='bwr')
        m.set_array(anom1)
        m.set_clim(dmin, dmax)
        plt.colorbar(m, cax=cbaxes, label=clrbr_lbl,
                     format='%2.1f', extend=extnd, boundaries=lvls)

        # create plot output directory
        dir_year = '{}/{}/'.format(dir_plot_out, iea_yr)
        dir_plots = dir_year + 'sst_anom/'
        # check if directory exist, if it doesn't then create
        try:
            os.makedirs(dir_plots)
        except OSError:
            if not os.path.isdir(dir_plots):
                raise

        # --Save the figure
        fn_fig = '{}oc_miniESR_{}_mons_{}_{}_contour.png'.format(dir_plots, var_wnt[0], mon_bgn, mon_end)
        plt.savefig(fn_fig, dpi=300, bbox_inches='tight')

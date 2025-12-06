import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.gridspec as gridspec
from matplotlib import interactive
from matplotlib import ticker
interactive(True)


# turn of toolbar on fig, set font
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 11})
mpl.rcParams.update({'savefig.bbox': 'tight'})


def fun_xr_ds2IEA_contour3(ds1, nlvl1, nlvl2, ttl, intrp_type, xtck, ytck,
                           dir_plots, fn_plots):
    """contour the anom, window mean, window trend with IEA symbols
    input:
    1) ds1 = xarray dataset, created in 'create_oc_sst_S_spatial.py'
    2) nlvl1 = contourf intervals for the anom can be an array or list
        example:
            nlvl1=np.arange(-2,2.1,0.1) or
            nlvl1=[np.arange(-1,1.1,0.1), np.arange(-1,1.1,0.1), np.arange(-2,2.1,0.1)]
    3) nlvl2 = contour intervals for the anom can by an array or a list
        example:
            nlvl2 = np.arange(-4, 4.5, 0.5)
            nlvl2=[np.arange(-1,1.1,0.5), np.arange(-1,1.1,0.5), np.arange(-2,2.1,0.5)]
    4) ttl = colorbar title for the anom can be a list of size 1 or size 3
        example:
            ttl = ['Winter SST anom ($\degree$C)']
            ttl = ['Winter SST anom ($\degree$C)', 'Winter 5-yr Mean / SD', 'Winter 5-yr Mean /SD']
    5) intrp_dim = grid dimension of interpolation
       a) to do bi-linear interpolation give x and y dimensions in a tuple (dx,dy)
       b) if oi interpolation was done in octave then give file name 'oi_sample.nc'
    6) xtck = map xticks
    7) ytck = map yticks
    8) dir_plots = output directory for the figure
    9) fn_plots = output filename for the figure (should contain type 'test.png'
    """

    # ---------------------------------------------------------------------------------
    # BEGIN: CHANGE THESE
    # ---------------------------------------------------------------------------------
    # variables in xarray dataset
    ds1_var = ['coord_mtrx', 'data_mtrx', 'mrkr_mtrx']

    # anom, mean5, trend5
    num_map = 3

    # figure size
    fig_wdth = 11
    fig_hght = 3.0
    # ---------------------------------------------------------------------------------
    # END: CHANGE THESE
    # ---------------------------------------------------------------------------------

    # define fireice colormap
    cmap1 = mpl.colors.LinearSegmentedColormap.from_list(
        "", ["mediumblue", "dodgerblue", "cyan", "white", "yellow", "tomato", "red"])

    # subplots
    num_clmn = num_map
    num_row = 1
    gs1 = gridspec.GridSpec(num_row, num_clmn)

    # change spacing
    gs1.update(left=0.1, right=0.9, bottom=0.05,
               top=0.9, wspace=0.1, hspace=0.1)

    # coastlines
    land1 = cfeature.NaturalEarthFeature(
        'physical', 'land', '50m', edgecolor='gray',
        facecolor='lightgray', zorder=100, linewidth=0.5)

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

    # open new figure of set size
    plt.close()
    plt.figure(figsize=(fig_wdth, fig_hght))

    # if ttl is only length one then need to append mean5 and trend5 labels
    if len(ttl) == 1:
        ttl.append("5-yr Mean / SD")
        ttl.append("5-yr Trend / SD")

    # loop over the 3 map types
    for i in range(num_map):
        # data
        data1 = ds1[ds1_var[1]].data[:, i]

        # if intrp_type is a tuple then bilinear interp
        if isinstance(intrp_type, tuple):
            ind = np.isfinite(data1)
            data1g = griddata((lon[ind], lat[ind]),
                              data1[ind], (xv, yv), method='linear')

        # markers
        marker1 = ds1[ds1_var[2]].data[:, i]

        # set contour limits and extend nlvl1 so that there are no white
        # contour areas where values exceed the nlvl1 limits
        if isinstance(nlvl1, list):
            dmin1 = nlvl1[i][0]
            dmax1 = nlvl1[i][-1]
            dnlvl1 = np.nanmean(np.diff(nlvl1[i]))
            nlvl1_ext = np.arange(dmin1*2, dmax1*2 + dnlvl1, dnlvl1)
        else:
            # contour intervals, only have user set levels for the first plot
            if i > 0:
                max_lvl = np.nanmax(np.abs(data1g))
                rng1 = np.round(max_lvl*10)/10
                nlvl1 = np.arange(-1*rng1, rng1+0.05, 0.05)
                nlvl2 = np.arange(-2.5, 3, 0.5)
            dmin1 = nlvl1[0]
            dmax1 = nlvl1[-1]
            dnlvl1 = np.nanmean(np.diff(nlvl1))
            nlvl1_ext = np.arange(dmin1*2, dmax1*2 + dnlvl1, dnlvl1)

        # plot stations and contour data
        ax1 = plt.subplot(gs1[i],  projection=ccrs.PlateCarree())
        cf1 = plt.contourf(xv, yv, data1g, nlvl1_ext, vmin=dmin1,
                           vmax=dmax1, cmap=cmap1, linestyles='-')

        # contour and labels
        if isinstance(nlvl2, list):
            mpl.rcParams['contour.negative_linestyle'] = 'solid'
            cn1 = plt.contour(xv, yv, data1g, nlvl2[i],
                              colors='dimgray', linewidths=0.5)
        else:
            mpl.rcParams['contour.negative_linestyle'] = 'solid'
            cn1 = plt.contour(xv, yv, data1g, nlvl2,
                              colors='dimgray', linewidths=0.5)

        # label contours and rotate to 0
        labels1 = plt.clabel(cn1, inline=1, inline_spacing=-5,
                             fontsize=6, fmt='%3.1f', colors='black')

        # rotate and change integers (ie 1.0) to whole number (ie 1)
        for l in cn1.labelTexts:
            txt1 = l.get_text()
            if float(txt1).is_integer():
                l.set_text(txt1.split('.')[0])
            l.set_rotation(0)

        # index of less/greater than 1 sd
        in1 = np.where(marker1 == '.')
        plt.plot(lon[in1], lat[in1], marker='.',
                 linestyle='none', color='black', markersize=1)

        # index of smallest/greastest of time series
        in1 = np.where(marker1 == '+')
        plt.plot(lon[in1], lat[in1], marker='x',
                 linestyle='none', color='black', markersize=2,
                 markeredgewidth=0.3)

        # put ticks on all axis
        plt.tick_params(bottom=True, top=True, left=True, right=True,
                        labelbottom=False, labeltop=False, labelleft=False,
                        labelright=False)

        # land
        ax1.add_feature(land1)

        # gridlines and lat and lon labels
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
        # gl.xlabels_bottom = False
        # gl.ylabels_right = False
        gl.xlines = False
        gl.ylines = False
        gl.xlocator = mpl.ticker.FixedLocator(xtck)
        gl.ylocator = mpl.ticker.FixedLocator(ytck)
        gl.bottom_labels = False
        gl.left_labels = False
        gl.right_labels = False

        # no left lat labels for trend5 and mean5
        if i == 0:
            gl.left_labels = True
        if i == 2:
            gl.right_labels = True

        # format the lat and lon labels (ie put degN and degW)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # xtick and yticks
        plt.gca().set_xticks(xtck, crs=ccrs.PlateCarree())
        plt.gca().set_yticks(ytck, crs=ccrs.PlateCarree())

        # colorbar
        m = plt.cm.ScalarMappable(cmap=cmap1)
        m.set_array(data1g)
        m.set_clim(dmin1, dmax1)

        cb = plt.colorbar(m, ax=plt.gca(),
                          orientation="horizontal", shrink=0.95,
                          pad=0.04, label=ttl[i], aspect=30,
                          boundaries=np.arange(dmin1, dmax1+dnlvl1, dnlvl1))

        # only have maximum 5 intervals of colorbar ticks
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()

    # --create output directory for plots and corr netcdf files
    dir_fn = '{}{}'.format(dir_plots, fn_plots)

    # filename out for plots and corr netcdf files
    plt.savefig(dir_fn, dpi=300)

    return dir_fn

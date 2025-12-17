import os
import shutil
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import rcParams
from scipy.interpolate import griddata
from matplotlib import ticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER



# plot paramaters
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'RdYlBu_r',
    'axes.grid': False,
    'savefig.dpi': 300,  # to adjust notebook inline plot size
    'xtick.top':        False,  # shold the top and bottom have tick marks
    'xtick.bottom':     True,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'ytick.direction': 'out',
    'xtick.direction': 'out',
    'axes.labelsize': 10,  # fontsize for x and y labels 
    'axes.titlesize': 10,
    'font.size': 10,  # was 10
    'legend.fontsize': 10,  # was 10
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.figsize': [8.5, 11],
    'font.family': 'STIXGeneral',
    'mathtext.fontset': 'stix',
    'toolbar': 'None',
    'savefig.bbox': 'tight'
}
rcParams.update(params)

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2025

# prelim ctd OxAve_StaCorr data has a few problems
# eg winter 2023 90-110, Spring 2023  86.7-60 and 86.7-70
# other stations don't look great, filter by unusual large values at the three depths
chck_lrg = [10, 4.3, 10]

# chck_top = 15
# chck_bttm = 5


# dir of spatial IEA stats
dir_in = './data_x13/DissolvedOxygen/'

# the file names have z1 and z2 (depth range to find minimum values)
z1 = 0
z2 = 500

# --variable name
var_name = ['R_O2']
var_ttl = 'DO'

# --names of variables to create monthly and seasonal means
dpth_wnt = np.array([50, 150])

# variables for fun_xr_ds2IEA_contour3
intrp_type = (0.05, 0.05)
intrp_type = (0.015, 0.015)
ytck = np.arange(30.5, 36, 2)
xtck = np.arange(-124, -116, 2)

# seasons to contour
season_name = ['Winter', 'Spring', 'Summer']
season_lbl = ['1', '2', '3']

# output file name
file_pre_out = 'oc_do_spatial'

# ds variables
ds_var = ['coord_mtrx', 'data_mtrx', 'mrkr_mtrx', 'ts_mtrx']

# figure size
fig_wdth = 11
fig_hght = 8.5

# clmn label (data, min, depth of min)
nc_file_wnt = ['.nc', '_min.nc', '_z_min.nc']
nc_file_wnt = ['.nc', '.nc', '_bottom.nc']

num_map = 3

# have multiple colormaps
cmap_list = ['RdYlBu', 'RdYlBu', 'RdYlBu']

dx_shallow = np.array([0.2, -0.2, 0.1, 0.1, -0.0, 0.2, 0.1, 0.15, -0.1, 0.1, 0.2])
dy_shallow = np.array([0.0, 0.05, 0.1, 0.1, -0.15, -0.1, 0.1, 0.0, 0.1, 0.1, -0.1])

# --plot directory
dir_plot_out = './figures_gha/DissolvedOxygen/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_nc_file_wnt = len(nc_file_wnt)

# create plot output directory
dir_plots = '{}/{}/'.format(dir_plot_out, iea_yr)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_plots)
except OSError:
    if not os.path.isdir(dir_plots):
        raise

# CalCOFI DO can be at different depths, setup variable want based on depth
num_dpth = len(dpth_wnt)
var_wnt = [var_name[0]+'{}m'.format(x) for x in dpth_wnt]
num_var = len(var_wnt)

# length of input variables
num_season = len(season_lbl)
num_ds_var = len(ds_var)

# open the data and contour
fn_list = list()

sttn_lrg_list = list()
for i in range(num_season):
    # 50 m
    in_dpth_wnt1 = 0
    fn_dpth = '{}data_mn5_trnd5_qrtr{}_{}{}m_between_{}_{}m.nc'.format(
        dir_in, season_lbl[i], var_name[0], dpth_wnt[in_dpth_wnt1], z1, z2)
    ds1_dpth1 = xr.open_dataset(fn_dpth)

    # lon and lat
    lon = ds1_dpth1.coord_mtrx[:, 0].data - 360
    lat = ds1_dpth1.coord_mtrx[:, 1].data

    # 150 m
    in_dpth_wnt2 = 1
    fn_dpth = '{}data_mn5_trnd5_qrtr{}_{}{}m_between_{}_{}m.nc'.format(
        dir_in, season_lbl[i], var_name[0], dpth_wnt[in_dpth_wnt2], z1, z2)
    ds1_dpth2 = xr.open_dataset(fn_dpth)

    # bottom
    fn_bottom = '{}data_mn5_trnd5_qrtr{}_{}_bottom_between_{}_{}m.nc'.format(
        dir_in, season_lbl[i], var_name[0], z1, z2)
    ds1_bottom = xr.open_dataset(fn_bottom)


    # station depth
    sttn_depth = ds1_bottom['sttn_depth'].data


    # variable labels are kind of complicated, combination of season,
    # dpth_wnt, bottom
    var_lbl_wnt = list()
    var_lbl_wnt.append('{} {} at {} m (ml/L)'.format(
        season_name[i], var_ttl, dpth_wnt[in_dpth_wnt1]))
    var_lbl_wnt.append('{} {} at {} m (ml/L)'.format(
        season_name[i], var_ttl, dpth_wnt[in_dpth_wnt2]))
    var_lbl_wnt.append(
        '{} bottom {} (ml/L)'.format(season_name[i], var_ttl))

    # only interested in the 'data' and not mean5 and trend5, so construct new
    # data and marker matrix
    data_mtrx = np.zeros([ds1_dpth1.index.shape[0], num_map])
    mrkr_mtrx = np.zeros([ds1_dpth1.index.shape[0], num_map], dtype='str')

    data_mtrx[:, 0] = ds1_dpth1.data_mtrx[:, 0].data
    data_mtrx[:, 1] = ds1_dpth2.data_mtrx[:, 0].data
    data_mtrx[:, 2] = ds1_bottom.data_mtrx[:, 0].data

    dm_chck = np.copy(data_mtrx)
    dm_chck_mn = np.nanmean(dm_chck, axis=1)
    
    mrkr_mtrx[:, 0] = ds1_dpth1.mrkr_mtrx[:, 0].data
    mrkr_mtrx[:, 1] = ds1_dpth2.mrkr_mtrx[:, 0].data
    mrkr_mtrx[:, 2] = ds1_bottom.mrkr_mtrx[:, 0].data

    # subplots
    num_clmn = num_map
    num_row = 1
    gs1 = gridspec.GridSpec(num_row, num_clmn)

    # change spacing
    gs1.update(left=0.1, right=0.9, bottom=0.05,
               top=0.9, wspace=0.1, hspace=0.1)

    # coastlines
    land1 = cfeature.NaturalEarthFeature(
        'physical', 'land', '50m', edgecolor='lightgray',
        facecolor='lightgray', zorder=100, linewidth=0.5)


    # shallow stations
    in_shallow = np.where(sttn_depth < 500)[0]
    lon_shallow = lon[in_shallow]
    lat_shallow = lat[in_shallow]
    dpth_shallow = sttn_depth[in_shallow]

    # if intrp_type is a tuple then create grid for bilinear interp
    if isinstance(intrp_type, tuple):
        dx = intrp_type[0]
        dy = intrp_type[1]
        lon_vec = np.arange(np.nanmin(lon), np.nanmax(lon)+dx, dx)
        lat_vec = np.arange(np.nanmin(lat), np.nanmax(lat)+dy, dy)
        xv, yv = np.meshgrid(lon_vec, lat_vec)

    # open new figure of set size
    plt.close()
    fig = plt.figure(figsize=(fig_wdth, fig_hght))

    # n levels
    min_data = np.round(np.nanmin(data_mtrx))
    max_data = np.ceil(np.nanmax(data_mtrx))

    # check unusual surface or bottom depths
    
    chck_data = np.isfinite(data_mtrx).nonzero()[0]
    if len(chck_data) > 0:
        # if max(data_mtrx[:, 0]) > chck_top:
        #     in_lrg = np.where(data_mtrx[:, 0] > chck_top)
        #     sttn_lrg = ds1_dpth1['sttn_mtrx'][in_lrg[0][0],:].data
        #     print("Remove from map: {} {} --> large surface value".format(sttn_lrg[0], sttn_lrg[1]))
        #     data_mtrx[in_lrg[0][0], :] = np.nan*np.zeros(3)

        # if max(data_mtrx[:, 2]) > chck_bttm:
        #     in_lrg = np.where(data_mtrx[:, 2] > chck_bttm)
        #     sttn_lrg = ds1_dpth1['sttn_mtrx'][in_lrg[0][0],:].data
        #     print("Remove from map: {} {} --> large bottom value".format(sttn_lrg[0], sttn_lrg[1]))
        #     data_mtrx[in_lrg[0], :] = np.nan*np.zeros(3)

        in_lrg = np.where(dm_chck_mn > chck_lrg[i])[0]
        if len(in_lrg) > 0:
            data_mtrx[in_lrg, :] = np.nan
            sttn_lrg = ds1_dpth1['sttn_mtrx'][in_lrg,:].data
            sttn_lrg_list.append(sttn_lrg)

    # need to recalculate min/max as data_mtrx might have changed
    # if a cast was removed
    min_data = np.round(np.nanmin(data_mtrx))
    max_data = np.ceil(np.nanmax(data_mtrx))
    nlvl1 = np.arange(min_data, max_data+0.2, 0.2)

    for iii in range(num_map):
        # lon and lat
        lon = ds1_dpth1.coord_mtrx[:, 0].data - 360
        lat = ds1_dpth1.coord_mtrx[:, 1].data

        # data
        data1 = data_mtrx[:, iii]

        # if intrp_type is a tuple then bilinear interp
        if isinstance(intrp_type, tuple):
            ind = np.isfinite(data1)
            numd = len(ind.nonzero()[0])
            if numd > 0:
                data1g = griddata((lon[ind], lat[ind]),
                                  data1[ind], (xv, yv), method='linear')
            else:
                data1g = xv*0.0

        # markers
        marker1 = mrkr_mtrx[:, iii]

        # plot stations and contour data
        ax1 = plt.subplot(gs1[iii],  projection=ccrs.PlateCarree())

        cf1 = plt.contourf(xv, yv, data1g, nlvl1,
                           cmap=cmap_list[iii], linestyles='-')

        # contour and labels
        mpl.rcParams['contour.negative_linestyle'] = 'solid'

        cn1 = plt.contour(xv, yv, data1g, [1.4, 1.401],
                          colors='dimgray', linewidths=0.5)

        # label contours and rotate to 0
        labels1 = plt.clabel(cn1, inline=1, inline_spacing=-5,
                             fontsize=5, fmt='%3.1f', colors='black')

        # rotate and change integers (ie 1.0) to whole number (ie 1)
        for l in cn1.labelTexts:
            txt1 = l.get_text()
            if float(txt1).is_integer():
                l.set_text(txt1.split('.')[0])
            l.set_rotation(0)

        # land
        ax1.add_feature(land1)

        marker1 = marker1[ind]
        lon = lon[ind]
        lat = lat[ind]
        # index of less than 1 sd (-)
        in1 = np.where(marker1 == '\u2212')[0]
        for jjj in range(len(in1)):
            plt.text(lon[in1][jjj], lat[in1][jjj], '\u2212',
                     color='black', ha='center', va='center', fontsize=5)
        # index of greater than 1 sd (+)
        in1 = np.where(marker1 == '\u002B')[0]
        for jjj in range(len(in1)):
            plt.text(lon[in1][jjj], lat[in1][jjj], '\u002B',
                     color='black', ha='center', va='center', fontsize=5)
        # index of smallest of time series (circle minus)
        in1 = np.where(marker1 == '\u2296')[0]
        for jjj in range(len(in1)):
            plt.text(lon[in1][jjj], lat[in1][jjj], '\u2296',
                     color='black', ha='center', va='center', fontsize=5)
        # index of largest of time series (circle plus)
        in1 = np.where(marker1 == '\u2295')[0]
        for jjj in range(len(in1)):
            plt.text(lon[in1][jjj], lat[in1][jjj], '\u2295',
                     color='black', ha='center', va='center', fontsize=5)

        # otherwise just plot a 'bullet' to mark station location
        in1 = np.where(marker1 == '')[0]
        for jjj in range(len(in1)):
            plt.text(lon[in1][jjj], lat[in1][jjj], '\u2022',
                     color='black', ha='center', va='center', fontsize=5)

        # otherwise just plot a 'bullet' to mark station location
        in1 = np.where(marker1 == ' ')[0]
        for jjj in range(len(in1)):
            plt.text(lon[in1][jjj], lat[in1][jjj], '\u2022',
                     color='black', ha='center', va='center', fontsize=5)


        # add depth labels for bottom depth
        if iii == 2:
            for l in range(len(in_shallow)):
                plt.text(lon_shallow[l]+dx_shallow[l],
                         lat_shallow[l]+dy_shallow[l],
                         int(dpth_shallow[l]),
                         color='green', zorder=1000, rotation=0,
                         ha='center', va='center', fontsize=5)

        # gridlines and lat and lon labels
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
        gl.bottom_labels = False
        gl.ylabels_right = False
        gl.xlines = False
        gl.ylines = False
        gl.xlocator = mpl.ticker.FixedLocator(xtck)
        gl.ylocator = mpl.ticker.FixedLocator(ytck)
        gl.left_labels = False
        gl.right_labels = False

        # no left lat labels for trend5 and mean5
        if iii == 0:
            gl.left_labels = True
        if iii == 2:
            gl.right_labels = True

        # put right lat labesl for mean5
        if iii == 2:
            gl.right_labels = True

        # format the lat and lon labels (ie put degN and degW)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # xtick and yticks
        plt.gca().set_xticks(xtck, crs=ccrs.PlateCarree())
        plt.gca().set_yticks(ytck, crs=ccrs.PlateCarree())

        # put ticks on all axis
        plt.tick_params(bottom=True, top=True, left=True, right=True,
                        labelbottom=False, labeltop=False, labelleft=False,
                        labelright=False)

        # colorbar
        ttl_clbr = var_lbl_wnt[iii]
        cb = plt.colorbar(cf1, ax=plt.gca(),
                          orientation="horizontal", shrink=0.95,
                          pad=0.04, label=ttl_clbr, aspect=30)

        # only have maximum 5 intervals of colorbar ticks
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()

        # x and y limits
        xlm = [np.nanmin(lon)-0.1, np.nanmax(lon)+0.1]
        ylm = [np.nanmin(lat)-0.1, np.nanmax(lat)+0.1]
        plt.xlim(xlm)
        plt.ylim(ylm)

    # --create output directory for plots and corr netcdf files
    dir_fn = '{}data_{}_{}_{}m_{}m_bottom.png'.format(
        dir_plots, var_name[0], season_name[i], dpth_wnt[in_dpth_wnt1],
        dpth_wnt[in_dpth_wnt2])

    # filename out for plots and corr netcdf files
    plt.savefig(dir_fn)


# remove the directory and files that has the downloaded
shutil.rmtree(dir_in)

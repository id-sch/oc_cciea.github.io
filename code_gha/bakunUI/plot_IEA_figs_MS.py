import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from fun_pd_df2IEA_fig_blue import fun_pd_df2IEA_fig_blue


# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Times New Roman']})

# turn of toolbar on fig, set font
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
rcParams.update({'font.size': 12})



# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2025

# --lat
# lat_wnt = [33, 39, 45]
# lat_str = ['33', '39', '45']
lat_wnt = [45, 39, 33]
lat_str = ['45', '39', '33']


# --input directory
dir_M = './csv_for_erddap/'
dir_S = './csv_for_erddap/'

# --IEA file names
file_pre = 'oc_ui'
num_pre = len(file_pre)

# --IEA year clim
yr_clim_bgn = 1967
yr_clim_end = iea_yr

# --IEA window size
wndw = 5
yy_wnt = yr_clim_end

# --subplot info
nr = 4
nc = 1

# --plot directory
dir_plot_out = './figures_gha/bakunUI/'
fig_type = '.png'


# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_lat = len(lat_wnt)

# create plot output directory
dir_plots = '{}{}/'.format(dir_plot_out, iea_yr)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_plots)
except OSError:
    if not os.path.isdir(dir_plots):
        raise


# Monthly
file_M = '{}{}_M.csv'.format(dir_M, file_pre)
dfM = pd.read_csv(file_M)
ordr_list = [[1], [2], [3]]
for i in range(len(ordr_list)):
    plt.clf()
    fun_pd_df2IEA_fig_blue(
        dfM, nr, nc, ordr_list[i], yr_clim_bgn, yr_clim_end, wndw, yy_wnt, marker_flag=0)
    fn_fig_M = '{}{}_{}_Monthly{}'.format(dir_plots, file_pre, lat_wnt[i], fig_type)
    plt.savefig(fn_fig_M, dpi=300, bbox_inches='tight')

# Seasons
file_S = '{}{}_S.csv'.format(dir_S, file_pre)
dfS = pd.read_csv(file_S)
ordr_list = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]
for i in range(len(ordr_list)):
    plt.clf()
    fun_pd_df2IEA_fig_blue(
        dfS, nr, nc, ordr_list[i], yr_clim_bgn, yr_clim_end, wndw, yy_wnt)
    fn_fig_S = '{}{}_{}_Season{}'.format(dir_plots, file_pre, lat_wnt[i], fig_type)
    plt.savefig(fn_fig_S, dpi=300, bbox_inches='tight')

# Plot monthly on one figure
fun_pd_df2IEA_fig_blue(dfM, nr, nc, [1, 2, 3], yr_clim_bgn,
                       yr_clim_end, wndw, yy_wnt, marker_flag=0)

lat_all = '_'.join(lat_str)
fn_fig_M_all = '{}{}_{}_Monthly{}'.format(
    dir_plots, file_pre, lat_all, fig_type)

plt.savefig(fn_fig_M_all, dpi=300, bbox_inches='tight')

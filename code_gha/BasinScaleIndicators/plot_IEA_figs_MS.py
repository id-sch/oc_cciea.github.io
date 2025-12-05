import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from fun_pd_df2IEA_fig_blue import fun_pd_df2IEA_fig_blue
# pylint: disable=C0103


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

# --basin index wanted
basin_wnt = ['ONI', 'PDO', 'NPGO']
num_basin = len(basin_wnt)

# --input directory
dir_M = './csv_for_erddap/'
dir_S = './csv_for_erddap/'

# --IEA file names
file_pre = 'oc'
num_pre = len(file_pre)

# --IEA year clim
yr_clim_bgn = 1950
# yr_clim_bgn = 1990
yr_clim_end = iea_yr

# --IEA window size
wndw = 5
yy_wnt = yr_clim_end

# --subplot info
nr = 4
nc = 1

# --plot directory


# dir_plots = '/home/isaac/data_files/Work/IEA/2018/fig_samples_for_SoCC/'
fig_type = '.png'


# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

# create plot output directory
dir_plots = '{}{}/'.format(dir_plot_out, iea_yr)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_plots)
except OSError:
    if not os.path.isdir(dir_plots):
        raise

for i in range(0, num_basin):

    # A) Monthly
    # --open Monthly
    file_M = '{}{}_{}_M.csv'.format(dir_M, file_pre, basin_wnt[i])
    dfM = pd.read_csv(file_M)
    fun_pd_df2IEA_fig_blue(
        dfM, nr, nc, [1], yr_clim_bgn, yr_clim_end, wndw, yy_wnt)

    fn_fig_M = '{}{}_{}_Monthly{}'.format(
        dir_plots, file_pre, basin_wnt[i], fig_type)
    plt.savefig(fn_fig_M, dpi=300, bbox_inches='tight')

    # B) seasons
    # --open Seasons
    file_S = '{}{}_{}_S.csv'.format(dir_S, file_pre, basin_wnt[i])

    dfS = pd.read_csv(file_S)
    fun_pd_df2IEA_fig_blue(
        dfS, nr, nc, [1, 2, 3, 4], yr_clim_bgn, yr_clim_end, wndw, yy_wnt)

    fn_fig_S = '{}{}_{}_Season{}'.format(
        dir_plots, file_pre, basin_wnt[i], fig_type)

    plt.savefig(fn_fig_S, dpi=300, bbox_inches='tight')


# C) all on one fig

# plot monthly on one figure
for i in range(0, num_basin):
    file_M = '{}{}_{}_M.csv'.format(dir_M, file_pre, basin_wnt[i])
    dfMi = pd.read_csv(file_M)
    nt = dfMi.shape[0]
    dfMi['order'].values[0:nt] = i+1

    if i == 0:
        dfM_all = dfMi
    else:
        dfM_all = pd.concat([dfM_all, dfMi], sort=False)

fun_pd_df2IEA_fig_blue(dfM_all, nr, nc, [1, 2, 3], yr_clim_bgn,
                       yr_clim_end, wndw, yy_wnt, marker_flag=0)

basin_all = '_'.join(basin_wnt)
fn_fig_M_all = '{}{}_{}_Monthly{}'.format(
    dir_plots, file_pre, basin_all, fig_type)

plt.savefig(fn_fig_M_all, dpi=300, bbox_inches='tight')

import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from fun_pd_df2IEA_fig_threshold_blue import fun_pd_df2IEA_fig_threshold_blue
from fun_mkdir import fun_mkdir


mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 12})


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2025

# --station wanted, two stations newport hydrographic at 5 and 25 km
sttn_wnt = ['NH05', 'NH25']
num_sttn = len(sttn_wnt)

# --input directory
dir_M = './csv_for_erddap/'
dir_S = './csv_for_erddap/'

# --IEA file names
file_pre = 'oc_arg_Newport'
num_pre = len(file_pre)

# --IEA year clim
yr_clim_bgn = 1998
yr_clim_end = iea_yr

# --IEA window size
wndw = 5
yy_wnt = yr_clim_end

# --subplot info
nr = 4
nc = 1

# aragonite threshold
threshold = 1

# --plot directory
dir_plot_out = './figures_gha/aragoniteNH/'
fig_type = '.png'

# figure size
fig_wdth = 11
fig_hght = 8.5
# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

# open new figure
plt.close()
fig = plt.figure(figsize=(fig_wdth, fig_hght))

# create plot output directory
dir_plots = '{}{}/'.format(dir_plot_out, iea_yr)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_plots)
except OSError:
    if not os.path.isdir(dir_plots):
        raise

# --open Monthly
file_M = '{}{}_M.csv'.format(dir_M, file_pre)
dfM_index = pd.read_csv(file_M)
dfM = dfM_index.rename(index=str, columns={"index": "data"})

for i in range(0, num_sttn):
    plt.clf()
    fun_pd_df2IEA_fig_threshold_blue(
        dfM, nr, nc, [i+1], yr_clim_bgn, yr_clim_end, wndw, yy_wnt, threshold)

    fn_fig_M = '{}{}_{}_Monthly_threshold{}'.format(
        dir_plots, file_pre, sttn_wnt[i], fig_type)
    plt.savefig(fn_fig_M, dpi=300, bbox_inches='tight')

    plt.clf()
    # --open Seasons
    file_S = '{}{}_{}_S.csv'.format(dir_S, file_pre, sttn_wnt[i])
    dfS_index = pd.read_csv(file_S)
    dfS = dfS_index.rename(index=str, columns={"index": "data"})
    fun_pd_df2IEA_fig_threshold_blue(dfS, nr, nc, [1, 2, 3, 4],
                                yr_clim_bgn, yr_clim_end, wndw, yy_wnt, threshold)

    fn_fig_S = '{}{}_{}_Season_threshold{}'.format(
        dir_plots, file_pre, sttn_wnt[i], fig_type)
    plt.savefig(fn_fig_S, dpi=100, bbox_inches='tight')

# remove the directory and files that has the downloaded
os.remove("./NH05_AragSat.nc")
os.remove("./NH25_AragSat.nc")

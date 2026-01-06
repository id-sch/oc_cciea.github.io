import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from fun_pd_df2IEA_fig_blue import fun_pd_df2IEA_fig_blue


# turn of toolbar on fig, set font
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 12})


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# end year
iea_yr = 2025

# regions
rgn = np.arange(1, 5)

# iea directory of the current yearly update
dir_iea = 'HCI_coastwide'

# --input directory
dir_M = './csv_for_erddap/'
dir_S = './csv_for_erddap/'

# --IEA file names
file_pre = 'oc_hci'
num_pre = len(file_pre)

# --IEA year clim
yr_clim_bgn = 1982
yr_clim_end = iea_yr

# --IEA window size
wndw = 5
yy_wnt = yr_clim_end

# --subplot info
nr = 4
nc = 1

# --plot directory
dir_plot_out = './figures_gha/HCI/'
fig_type = '.png'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# len of input variables
num_data = len(rgn)

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

for i in range(num_data):
    plt.clf()
    fun_pd_df2IEA_fig_blue(
        dfM, nr, nc, [i+1], yr_clim_bgn, yr_clim_end, wndw, yy_wnt, marker_flag=0)
    fn_fig_M = '{}{}_rgn{}_Month{}'.format(dir_plots, file_pre, rgn[i], fig_type)
    print('i={}, filename={}'.format(i, fn_fig_M))
    plt.savefig(fn_fig_M)


# Seasons
file_S = '{}{}_S.csv'.format(dir_S, file_pre)
dfS = pd.read_csv(file_S)
ordr_list = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]]
for i in range(len(ordr_list)):
    plt.clf()
    fun_pd_df2IEA_fig_blue(
        dfS, nr, nc, ordr_list[i], yr_clim_bgn, yr_clim_end, wndw, yy_wnt)
    fn_fig_S = '{}{}_rgn{}_Season{}'.format(dir_plots, file_pre, rgn[i], fig_type)
    plt.savefig(fn_fig_S, dpi=300, bbox_inches='tight')


# All on one fig
# plot monthly on one figure
fun_pd_df2IEA_fig_blue(dfM, nr, nc, [1, 2, 3, 4], yr_clim_bgn,
                       yr_clim_end, wndw, yy_wnt, marker_flag=0)

rgn_all = '_'.join(rgn.astype('str'))
fn_fig_M_all = '{}{}_rgns_{}_Monthly{}'.format(
    dir_plots, file_pre, rgn_all, fig_type)

plt.savefig(fn_fig_M_all)

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from fun_pd_df2IEA_fig_threshold_blue import fun_pd_df2IEA_fig_threshold_blue


# turn of toolbar on fig, set font
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 12})
# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2025

# --station wanted
sttn_wnt = ['933_300', '900_900', '800_800']

# --input directory
dir_Q = './data_x13/DissolvedOxygen/'

# --IEA file names
file_pre = 'oc_do_calcofi'
num_pre = len(file_pre)

# --IEA year clim
yr_clim_bgn = 1984
yr_clim_end = iea_yr

# --IEA window size
wndw = 5
yy_wnt = yr_clim_end

# --subplot info
nr = 4
nc = 1

# do threshold
threshold = 1.4

# --plot directory
dir_plot_out = './figures_gha/DissolvedOxygen/'
fig_type = '.png'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_sttn = len(sttn_wnt)

# figure size
fig_wdth = 11
fig_hght = 8.5

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

# Quarterly
file_Q = '{}{}_Q.csv'.format(dir_Q, file_pre)
dfQ = pd.read_csv(file_Q)
ordr_list = [[1], [2], [3]]
for i in range(len(ordr_list)):
    plt.clf()
    fun_pd_df2IEA_fig_threshold_blue(
        dfQ, nr, nc, ordr_list[i], yr_clim_bgn, yr_clim_end, wndw, yy_wnt, threshold, marker_flag=0)
    fn_fig_Q = '{}{}_{}_Quarter{}'.format(dir_plots, file_pre, sttn_wnt[i], fig_type)
    plt.savefig(fn_fig_Q, dpi=300, bbox_inches='tight')

# Plot monthly on one figure
fun_pd_df2IEA_fig_threshold_blue(
    dfQ, nr, nc, [1, 2, 3], yr_clim_bgn, yr_clim_end, wndw, yy_wnt, threshold, marker_flag=0)

sttn_wnt_all = '_'.join(sttn_wnt)
fn_fig_Q_all = '{}{}_{}_Quarter{}'.format(
    dir_plots, file_pre, sttn_wnt_all, fig_type)

plt.savefig(fn_fig_Q_all, dpi=300, bbox_inches='tight')

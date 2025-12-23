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

# iea directory of the current yearly update
dir_iea = 'nph'

# --input directory
dir_M = ['csv_for_erddap/']
dir_S = ['csv_for_erddap/']

# --IEA file names
file_pre = 'oc_nph'
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
# dir_plots = 'plots/'
dir_plot_out = './figures_gha/NPH/'
fig_type = '.png'


# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# len of input variables
num_data = len(dir_M)

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


for i in range(0, num_data):
    plt.clf()
    # --open Monthly
    file_M = '{}{}_M.csv'.format(dir_M[i], file_pre)
    dfM = pd.read_csv(file_M)
    fun_pd_df2IEA_fig_blue(dfM, nr, nc, [1], yr_clim_bgn, yr_clim_end, wndw, yy_wnt, marker_flag=0)

    fn_fig_M = '{}{}_Month{}'.format(dir_plots, file_pre, fig_type)

    plt.savefig(fn_fig_M, dpi=300, bbox_inches='tight')

    # --open Seasons
    file_S = '{}{}_S.csv'.format(dir_S[i], file_pre)

    dfS = pd.read_csv(file_S)
    fun_pd_df2IEA_fig_blue(dfS, nr, nc, [1, 2, 3, 4],
                      yr_clim_bgn, yr_clim_end, wndw, yy_wnt)
    fn_fig_S = '{}{}_Season{}'.format(dir_plots, file_pre, fig_type)
    plt.savefig(fn_fig_S, dpi=300, bbox_inches='tight')

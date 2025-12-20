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
# sttn_wnt = ['do_933_300', 'do_900_900']
sttn_wnt = ['933_300', '900_900', '800_800']
# sttn_wnt = ['867_450', '800_900']

num_sttn = len(sttn_wnt)

# --input directory
dir_S = './data_x13/csv_files/'

# --IEA file names
file_pre = 'oc_do'
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

for i in range(0, num_sttn):
    # --open Seasons
    file_S = dir_S + file_pre + '_' + sttn_wnt[i] + '_S.csv'

    dfS_index = pd.read_csv(file_S)
    dfS = dfS_index.rename(columns={'index': 'data'})
    fun_pd_df2IEA_fig_threshold_blue(dfS, nr, nc, [1, 2, 3, 4],
                                     yr_clim_bgn, yr_clim_end, wndw, yy_wnt,
                                     threshold)

    depth = str(int(dfS['depth'][0]))
    fn_fig_S = dir_plots + file_pre + '_' + sttn_wnt[i] + '_' + depth + '_Season' + fig_type
    plt.savefig(fn_fig_S, dpi=300, bbox_inches='tight')

# remove the directory and files that has the downloaded
shutil.rmtree(dir_in)

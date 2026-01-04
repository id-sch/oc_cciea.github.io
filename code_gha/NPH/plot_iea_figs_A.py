import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import rcParams
from fun_pd_df2IEA_fig_blue import fun_pd_df2IEA_fig_blue



rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'
rcParams.update({'font.size': 12})

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# end year
iea_yr = 2024

# iea directory of the current yearly update
dir_iea = 'nph'

# --input directory
# dir_A = ['/home/isaac/data_files/Work/IEA/{}/{}/'.format(iea_yr, dir_iea)]
dir_A = ['./csv_for_erddap/']

# --IEA file names
file_pre = 'oc_nph_jf'
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
dir_plot_out = './figures_gha/NPH/'
fig_type = '.png'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
# len of input variables
num_data = len(dir_A)

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

    # --open Seasons
    file_A = '{}{}_A.csv'.format(dir_A[i], file_pre)

    dfA = pd.read_csv(file_A)
    fun_pd_df2IEA_fig_blue(dfA, nr, nc, [1],
                           yr_clim_bgn, yr_clim_end, wndw, yy_wnt)
    fn_fig_A = '{}{}_Annual{}'.format(dir_plots, file_pre, fig_type)
    plt.savefig(fn_fig_A, dpi=300, bbox_inches='tight')

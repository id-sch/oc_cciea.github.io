import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from fun_pd_df2IEA_fig_blue import fun_pd_df2IEA_fig_blue


# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Times New Roman']})

# turn of toolbar on fig, set font
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 12})

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2024

# --basin index wanted
pheno_wnt = ['sti', 'lusi', 'tumi', 'fti']
num_basin = len(pheno_wnt)

# --input directory
dir_in = './csv_for_erddap/'

# --IEA file names
file_pre = 'oc'
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
dir_plot_out = './figures_gha/phenologyUI/'
fig_type = '.png'


# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_pheno = len(pheno_wnt)

# create plot output directory
dir_plots = '{}{}/'.format(dir_plot_out, iea_yr)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_plots)
except OSError:
    if not os.path.isdir(dir_plots):
        raise

# A) annual
for i in range(0, num_pheno):
    file_A = '{}{}_{}_A.csv'.format(dir_in, file_pre, pheno_wnt[i])

    dfA = pd.read_csv(file_A)
    fun_pd_df2IEA_fig_blue(dfA, nr, nc, [3, 2, 1], yr_clim_bgn, yr_clim_end, wndw, yy_wnt)

    fn_fig_A = '{}{}_{}_Annual{}'.format(
        dir_plots, file_pre, pheno_wnt[i], fig_type)

    plt.savefig(fn_fig_A, dpi=300, bbox_inches='tight')

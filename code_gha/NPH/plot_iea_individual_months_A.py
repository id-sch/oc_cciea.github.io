import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import calendar as clndr
from fun_pd_df2IEA_fig_blue import fun_pd_df2IEA_fig_blue


mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams.update({'font.size': 12})

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# end year
iea_yr = 2024

# iea directory of the current yearly update
dir_iea = 'nph'

# --input directory
dir_M = ['./csv_for_erddap/']

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
dir_plot_out = './figures_gha/NPH/months/'
fig_type = '.png'

# months
mon = np.arange(1, 13)
# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_mon = len(mon)

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


# --open Monthly
file_M = '{}{}_M.csv'.format(dir_M[0], file_pre)
dfM = pd.read_csv(file_M)

for i in range(num_mon):
    plt.clf()
    dfM_mon_indx_lbl = dfM.loc[dfM['month'] == mon[i]]
    dfM_mon = dfM_mon_indx_lbl.rename(columns={'index': 'data'})
    new_ttl = '{} North Pacific High Area'.format(clndr.month_name[mon[i]])
    dfM_mon['timeseries'].replace({dfM_mon.timeseries.values[0]: new_ttl}, inplace=True)


    fun_pd_df2IEA_fig_blue(dfM_mon, nr, nc, [1],
                      yr_clim_bgn, yr_clim_end, wndw, yy_wnt)
    fn_fig_mon = '{}{}_Mon{}{}'.format(dir_plots, file_pre, mon[i], fig_type)
    plt.savefig(fn_fig_mon, dpi=300, bbox_inches='tight')

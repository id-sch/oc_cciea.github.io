import os
import numpy as np
import xarray as xr
from fun_xr_ds2IEA_contour3 import fun_xr_ds2IEA_contour3


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# iea year
iea_yr = 2025

# dir of spatial IEA stats
dir_in = '/home/isaac/data_files/Work/IEA/{}/spatial/oi_sst/'.format(iea_yr)

# variables for fun_xr_ds2IEA_contour3
# nlvl for old contour figures
nlvl1_old = [np.arange(-3, 3.1, 0.1),
             np.arange(-3, 3.1, 0.1),
             np.arange(-3, 3.1, 0.1)] 
nlvl2_old = [np.arange(-3, 3.5, 0.5),
             np.arange(-3, 3.5, 0.5),
             np.arange(-3, 3.5, 0.5)]

ttl = ['SST anom ($\degree$C)', '5-yr Mean / SD', '5-yr Trend /SD']
intrp_type = (1, 1)
ytck = np.arange(35, 65, 10)
xtck = np.arange(-160, -100, 20)

# seasons to contour
season_lbl = ['winter', 'spring', 'summer', 'fall']

# output file name
file_pre_out = 'oc_sst_spatial'

# ds variables
ds_var = ['coord_mtrx', 'data_mtrx', 'mrkr_mtrx', 'ts_mtrx']

# example timeseries variables
freq_wnt = 'AS'
wndw = 5
yy_end = iea_yr
yr_clim_bgn = 1982
yr_clim_end = yy_end

# figure size
fig_wdth = 11
fig_hght = 8.5

# column label
clmn_lbl = ['min', 'max']

# row label
row_lbl = ['anom', 'mean5', 'trend5']

# --plot directory
dir_plot_out = './figures_gha/SSTspatial/'
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

# length of input variables
num_season = len(season_lbl)
num_ds_var = len(ds_var)

# open the data and contour
fn_list = []
for i in range(num_season):
    fn_in_ssn = '{}anom_mn5_trnd5_{}_clim_{}_{}.nc'.format(
        dir_in, season_lbl[i], yr_clim_bgn, yr_clim_end)
    ds1 = xr.open_dataset(fn_in_ssn)

    # contour the 3 IEA maps old
    fn_out_ssn_old = '{}_{}_clim_{}_{}_old.png'.format(
        file_pre_out, season_lbl[i], yr_clim_bgn, yr_clim_end)

    fn_out_ssn_old = '{}_{}_clim_{}_current_year.png'.format(
        file_pre_out, season_lbl[i], yr_clim_bgn)

    ttl_clrbar = ['{} {}'.format(season_lbl[i].capitalize(), ttl[0]),
                  '{} {}'.format(season_lbl[i].capitalize(), ttl[1]),
                  '{} {}'.format(season_lbl[i].capitalize(), ttl[2])]
    fn1 = fun_xr_ds2IEA_contour3(
        ds1, nlvl1_old, nlvl2_old, ttl_clrbar, intrp_type, xtck, ytck,
        dir_plots, fn_out_ssn_old)
    fn_list.append(fn1)

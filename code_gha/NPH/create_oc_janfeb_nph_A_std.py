import xarray as xr
import pandas as pd
import numpy as np
from fun_pd_df2csvR_time import fun_pd_df2csvR_time



# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------

# --input directory
fn_in = ['./data_gha/NPH/nph_area_jan_feb.nc']

# --variable
var_df = ['Area Jan-Feb mean']
sd_df = ['Area Jan-Feb sd']

# --IEA file names
file_pre = 'oc_nph_jf'

# --IEA labels
y_lbl = ['Area (10^6 km^2)']
ts_lbl = ['North Pacific High January - February Area']


# --IEA file columns
lat = [np.nan]
lon = [np.nan]
depth = [np.nan]

# --names of time, depth, lat dimensions
var_time = 'time'

# --output CSV file
yr_csv_bgn = 1967
dir_out = './csv_for_erddap/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_data = len(fn_in)

# create list of dataframes
df_list = list()
for i in range(0, num_data):
    # --open the netcdf file
    ds1 = xr.open_dataset(fn_in[i])

    # year
    yr_csv_end = ds1.year.data[-1]
    
    # make the xr.dataset into a pd.dataframe
    df1 = ds1.to_dataframe()

    # get the var and sd
    var1 = [var_df[i], sd_df[i]]
    dff = df1[var1]

    # --combine mean and sd into one pd.df
    dfA = dff.rename(index=str, columns={var_df[i]: "data"})
    dfA = dfA.rename(index=str, columns={sd_df[i]: "sd"})
    dfA['Datetime'] = pd.to_datetime(dfA.index)
    dfA = dfA.set_index('Datetime')

    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfQ, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list.append(dfA)


# write R style CSV
num_order = len(df_list)
fn_out = '{}_A.csv'.format(file_pre)
metric_lbl = y_lbl[i]
clmns_iea = ['year', 'time', 'index', 'error', 'SElo', 'SEup', 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']
fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_list, lat[0], lon[0], depth[0], metric_lbl, ts_lbl, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

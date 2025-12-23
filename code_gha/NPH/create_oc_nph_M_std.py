import xarray as xr
import pandas as pd
import numpy as np
from fun_pd_df2csvR_time import fun_pd_df2csvR_time


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------

# --input directory
fn_in = ['./data_gha/NPH/nph_area_month.nc']

# --variable
var_df = ['Area Monthly']

# --IEA file names
file_pre = 'oc_nph'

# --IEA labels
y_lbl = ['Area (10^6 km^2)']
ts_lbl = ['North Pacific High Area']

# --IEA file columns
lat = [np.nan]
lon = [np.nan]
depth = [np.nan]

# --names of time, depth, lat dimensions
var_time = 'time'

# --output CSV file
dir_out = './csv_for_erddap/'
yr_csv_bgn = 1967

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_data = len(fn_in)

df_list = list()
for i in range(0, num_data):
    # --open the netcdf file
    ds1 = xr.open_dataset(fn_in[i])

    # years
    yr_csv_end = ds1.time.dt.year.data[-1]
    
    # make the xr.dataset into a pd.dataframe
    df1 = ds1.to_dataframe()

    # --resample by taking 'M' means, i.e. monthly means
    dfMmn = df1.resample('ME').mean()
    dfMsd = df1.resample('ME').std()

    # --combine mean and sd into one pd.df
    dfM = dfMmn.rename(index=str, columns={var_df[i]: "data"})
    dfM['Datetime'] = pd.to_datetime(dfM.index)
    dfM = dfM.set_index('Datetime')
    dfM['sd'] = dfMsd.values

    ts_lbl1 = 'Monthly {}'.format(ts_lbl[i])
    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfM, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list = [dfM]
    num_order = len(df_list)

    fn_out = '{}_M.csv'.format(file_pre)

    metric_lbl = y_lbl[i]
    ts_lbl = [ts_lbl1]
    clmns_iea = ['year', 'month', 'time', 'index', 'error', 'SElo', 'SEup',
                 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

    fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_list, lat, lon, depth,
                                metric_lbl, ts_lbl, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

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
yr_csv_bgn = 1967
dir_out = './csv_for_erddap/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_data = len(fn_in)

for i in range(0, num_data):
    # --open the netcdf file
    ds1 = xr.open_dataset(fn_in[i])

    # year
    yr_csv_end = ds1.time.dt.year.data[-1]

    # make the xr.dataset into a pd.dataframe
    df1 = ds1.to_dataframe()

    # --resample by taking 'Q' means, i.e. quarterly means
    dfQmn = df1.resample('QE-MAR').mean()
    dfQsd = df1.resample('QE-MAR').std(ddof=0)

    # --combine mean and sd into one pd.df
    dfQ = dfQmn.rename(index=str, columns={var_df[i]: "data"})
    dfQ['Datetime'] = pd.to_datetime(dfQ.index)
    dfQ = dfQ.set_index('Datetime')
    dfQ['sd'] = dfQsd.values

    # -- seperate by season
    in1 = dfQ.index.month == 3
    in2 = dfQ.index.month == 6
    in3 = dfQ.index.month == 9
    in4 = dfQ.index.month == 12
    dfQ1 = dfQ[in1]
    dfQ2 = dfQ[in2]
    dfQ3 = dfQ[in3]
    dfQ4 = dfQ[in4]
    ts_lbl1 = 'Winter {}'.format(ts_lbl[i])
    ts_lbl2 = 'Spring {}'.format(ts_lbl[i])
    ts_lbl3 = 'Summer {}'.format(ts_lbl[i])
    ts_lbl4 = 'Fall {}'.format(ts_lbl[i])

    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfQ, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list = [dfQ1, dfQ2, dfQ3, dfQ4]
    num_order = len(df_list)

    fn_out = '{}_S.csv'.format(file_pre)

    metric_lbl = y_lbl[i]
    ts_lbl = [ts_lbl1, ts_lbl2, ts_lbl3, ts_lbl4]
    clmns_iea = ['year', 'time', 'index', 'error', 'SElo', 'SEup',
                 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

    fn_out_csv = fun_pd_df2csvR_time(clmns_iea, df_list, lat[i], lon[i], depth[i],
                                metric_lbl, ts_lbl, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

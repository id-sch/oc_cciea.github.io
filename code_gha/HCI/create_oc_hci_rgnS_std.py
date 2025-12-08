import pandas as pd
import numpy as np
from fun_pd_df2csvR_time import fun_pd_df2csvR_time


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# parameters for HCI coastwide
thrsh_scenario = 'mn_75kmM'
xdis = 150
lat_rgn = [[43.5, 48], [40, 43.5], [35.5, 40], [30, 35.5]]
lat_rgn_lbl = [[43.5, 48], [40, 43.5], [35.5, 40], [30, 35.5]]
clim_bgn = 1980
clim_end = 2010

# dir out (need this to set input directory, will be changed later)
dir_out = './data_gha/HCI/'

# input directory
dir_in = dir_out

# --variable
var_df = ['fraction below threshold "{}"'.format(thrsh_scenario)]

# --IEA file names
file_pre = 'oc_hci'

# --IEA labels
y_lbl = 'Habitat Area (fraction below monthly threshold)'
ts_lbl = 'Habitat Compression Index'

# --IEA file columns
lat = [np.nan]
lon = [np.nan]

# --names of time, depth, lat dimensions
var_time = 'time'

# --output CSV file
dir_out = './csv_for_erddap/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_data = len(lat_rgn)

for i in range(num_data):
    # lat region
    lat_rgn1 = lat_rgn[i]
    lat_rgn_lbl1 = lat_rgn_lbl[i]

    # CSV filename
    fn_in = '{}fraction_below_Threshold_{}_Lat{}_{}_X{}_Clim{}_{}.csv'.format(
        dir_in, thrsh_scenario, lat_rgn1[0], lat_rgn1[1], xdis, clim_bgn, clim_end)

    # --open the CSV file as pd dataframe
    df1 = pd.read_csv(fn_in)

    # change the time column to datetime
    dfi = df1.set_index(pd.to_datetime(df1['time']))

    # get start and end year
    yr_csv_bgn = pd.to_datetime(df1['time']).dt.year.values[0]
    yr_csv_end = pd.to_datetime(df1['time']).dt.year.values[-1]    

    # --resample by taking 'Q' means, i.e. quarterly means
    dfQmn = dfi.resample('QE-MAR').mean(numeric_only=True)
    dfQsd = dfi.resample('QE-MAR').std(ddof=0, numeric_only=True)

    # --combine mean and sd into one pd.df
    dfQ = dfQmn.rename(index=str, columns={var_df[0]: "data"})
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

    ts_lbl1 = 'Winter {} over {}-{}N'.format(ts_lbl, lat_rgn_lbl1[0], lat_rgn_lbl1[1])
    ts_lbl2 = 'Spring {} over {}-{}N'.format(ts_lbl, lat_rgn_lbl1[0], lat_rgn_lbl1[1])
    ts_lbl3 = 'Summer {} over {}-{}N'.format(ts_lbl, lat_rgn_lbl1[0], lat_rgn_lbl1[1])
    ts_lbl4 = 'Fall {} over {}-{}N'.format(ts_lbl, lat_rgn_lbl1[0], lat_rgn_lbl1[1])

    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfQ, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list = [dfQ1, dfQ2, dfQ3, dfQ4]
    num_order = len(df_list)

    fn_out = '{}_rgn{}S.csv'.format(file_pre, i+1)

    depth = 0
    metric_lbl = y_lbl
    ts_lbl_rgn = [ts_lbl1, ts_lbl2, ts_lbl3, ts_lbl4]
    clmns_iea = ['year', 'time', 'index', 'error', 'SElo', 'SEup',
                 'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

    fn_out_csv = fun_pd_df2csvR_time(
        clmns_iea, df_list, lat[0], lon[0], depth, metric_lbl, ts_lbl_rgn,
        dir_out, fn_out, yr_csv_bgn, yr_csv_end)

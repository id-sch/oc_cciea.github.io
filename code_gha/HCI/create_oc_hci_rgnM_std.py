import pandas as pd
import numpy as np
from fun_pd_df2csvR_time import fun_pd_df2csvR_time


# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# parameters for HCI coastwide
thrsh_scenario = 'mn_75kmM'
xdis = 150
# lat_rgn = [[43.5, 47.9], [40, 43.5], [35.5, 40], [30, 35.5]]
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
# y_lbl = 'Habitat Area\n(fraction below monthly threshold)'
y_lbl = 'Habitat Area (fraction below monthly threshold)'
ts_lbl = 'Habitat Compression Index'

# --IEA file columns
lon = [np.nan]
lat = [np.nan]


# --names of time, depth, lat dimensions
var_time = 'time'

# --output CSV file
dir_out = './csv_for_erddap/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------
num_data = len(lat_rgn)

depth = 0
metric_lbl = y_lbl
clmns_iea = ['year', 'month', 'time', 'index', 'error', 'SElo', 'SEup',
             'metric', 'timeseries', 'lat', 'lon', 'depth', 'order']

df_list = []
ts_lbl_rgn = []
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

    # --resample by taking 'M' means, i.e. monthly means
    dfMmn = dfi.resample('ME').mean(numeric_only=True)
    dfMsd = dfi.resample('ME').std(numeric_only=True)

    # --combine mean and sd into one pd.df
    dfM = dfMmn.rename(index=str, columns={var_df[0]: "data"})
    dfM['Datetime'] = pd.to_datetime(dfM.index)
    dfM = dfM.set_index('Datetime')
    dfM['sd'] = dfMsd.values

    # --Create monthly CSV, 'fun_pd_df2csvR' expects a list of pd.df. For
    # --the case of monthly means this list only consists of dfM, but a
    # --CSV for seasonal means will need a list of 4 pd.df for each season.
    df_list.append(dfM)

    ts_lbl_rgn.append('Monthly {} over {}-{}N'.format(ts_lbl, lat_rgn_lbl1[0], lat_rgn_lbl1[1]))

# filename is based on Lynn's Uploader tool name
fn_out = '{}_M.csv'.format(file_pre)

fn_out_csv, df_all = fun_pd_df2csvR_time(
    clmns_iea, df_list, lat[0], lon[0], depth, metric_lbl, ts_lbl_rgn, dir_out, fn_out, yr_csv_bgn, yr_csv_end)

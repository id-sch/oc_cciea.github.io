import os
import numpy as np
import pandas as pd


# --------------------------------------------------------------
# --Write "R" style CSV files for IEA time series
#
# The function expects a list of pd df, each member of the list
# is a unique IEA time series. The time series are distinguished
# by a "order" number
#
# Input variables:
# 1) column names of the desired CSV file
# 2) list of pd df, each df is  (df have columns of year/mon, mean, sd)
# 3) lat, lon, depth (single value variables or list)
# 4) metric, timeseries (strings (or list) of y-label and list of titles)
# 5) dir_out, fn_out (directory and filename of the CSV file)
# 6) yr_bgn, yr_end (subset the df by the year range set by these limits)
# Output:
# 1) filename of CSV file (dir and fn)
# --------------------------------------------------------------
def fun_pd_df2csvR_time(clmns_iea, df_list, lat, lon, depth, metric_lbl, timeseries_lbl, dir_out, fn_out, yr_bgn, yr_end):
    '''
    Write "R" style CSV files for IEA time series

    The function expects a list of pd df, each member of the list
    is a unique IEA time series. The time series are distinguished
    by a "order" number

    Input variables:
    1) column names of the desired CSV file
    2) list of pd df, each df is  (df have columns of year/mon, mean, sd)
    3) lat, lon, depth (single value variables or list)
    4) metric, timeseries (strings (or list) of y-label and list of titles)
    5) dir_out, fn_out (directory and filename of the CSV file)
    6) yr_bgn, yr_end (subset the df by the year range set by these limits)
    Output:
    1) filename of CSV file (dir and fn)
    '''

    # number of time series
    num_order = len(df_list)

    # initialize a "empty" pd df, this df will have "order"
    # number of time series
    df_all = pd.DataFrame()
    for iii in range(0, num_order):
        df_iii = df_list[iii]

        # get only years desired
        in_yrs = ((df_iii.index.year >= yr_bgn) &
                  (df_iii.index.year <= yr_end))
        df_yrs = df_iii[in_yrs]

        # --Check to see if month column is needed
        in_mon = np.where(np.array(clmns_iea) == 'month')[0]

        # --Check to see if day column is needed
        in_day = np.where(np.array(clmns_iea) == 'day')[0]

        # --Each unique IEA timeseries is given an order number, so create a
        #   seperate "order pd.df"
        # --1) Year
        in_mrkr = 1
        df_order = pd.DataFrame(
            df_yrs.index.year.values, columns=[clmns_iea[0]])
        yrs = df_yrs.index.year.values

        # --2) Month, if requested (monthly has it, annual doesn't)
        if len(in_mon) == 1:
            df_order[clmns_iea[in_mon[0]]] = df_yrs.index.month
            mons = df_yrs.index.month.values
            in_mrkr = 2

        # --3) Day, if requested (if 'All' data is requested,
        #      ie don't do monthly or seasonal means)
        if len(in_day) == 1:
            df_order[clmns_iea[in_day[0]]] = df_yrs.index.day
            days = df_yrs.index.day.values
            in_mrkr = 3

        # --3.5) Time
        tt = []
        if in_mrkr == 1:
            for jjj in range(len(yrs)):
                tt.append('{}'.format(yrs[jjj]))
        if in_mrkr == 2:
            for jjj in range(len(yrs)):
                tt.append('{}-{:02d}'.format(yrs[jjj], mons[jjj]))
        if in_mrkr == 3:
            for jjj in range(len(yrs)):
                tt.append(
                    '{}-{:02d}-{:02d}'.format(yrs[jjj], mons[jjj], days[jjj]))

        in_mrkr = in_mrkr + 1

        df_order['time'] = tt

        # --4) Data and error
        df_order[clmns_iea[in_mrkr]] = pd.Series(df_yrs.data.values)
        df_order[clmns_iea[in_mrkr+1]] = pd.DataFrame(df_yrs.sd.values)

        # --5) error low and error high
        df_order[clmns_iea[in_mrkr+2]
                 ] = pd.Series(df_yrs.data.values-df_yrs.sd.values)
        df_order[clmns_iea[in_mrkr+3]
                 ] = pd.Series(df_yrs.data.values+df_yrs.sd.values)

        # --6) metric and timeseries labels
        if isinstance(metric_lbl, str):
            lbl_array = np.chararray(
                df_order.shape[0], itemsize=len(metric_lbl), unicode=True)
            lbl_array[:] = metric_lbl
        if isinstance(metric_lbl, list):
            lbl_array = np.chararray(
                df_order.shape[0], itemsize=len(metric_lbl[iii]), unicode=True)
            lbl_array[:] = metric_lbl[iii]

        df_order[clmns_iea[in_mrkr+4]] = pd.Series(lbl_array)

        lbl_array1 = np.chararray(
            df_order.shape[0], itemsize=len(timeseries_lbl[iii]), unicode=True)
        lbl_array1[:] = timeseries_lbl[iii]
        df_order[clmns_iea[in_mrkr+5]] = pd.Series(lbl_array1)

        # --7) lat, lon, depth
        var_array = np.ndarray(df_order.shape[0])
        if isinstance(lat, list):
            var_array[:] = lat[iii]
        else:
            var_array[:] = lat
        df_order[clmns_iea[in_mrkr+6]] = pd.Series(var_array)

        if isinstance(lon, list):
            var_array[:] = lon[iii]
        else:
            var_array[:] = lon
        df_order[clmns_iea[in_mrkr+7]] = pd.Series(var_array)

        if isinstance(depth, list):
            var_array[:] = depth[iii]
        else:
            var_array[:] = depth
        df_order[clmns_iea[in_mrkr+8]] = pd.Series(var_array)

        # --8) order
        var_array[:] = iii+1
        df_order[clmns_iea[in_mrkr+9]] = pd.Series(var_array)

        # --9) append all "order pd.df" together
        # df_all = df_all.append(df_order)

        # Note: As of pandas 2.0, append (previously deprecated) was removed.
        #       You need to use concat instead (for most applications):
        df_all = pd.concat([df_all, df_order], ignore_index=True)

    # --check if directory exist, if it doesn't then create
    try:
        os.makedirs(dir_out)
    except OSError:
        if not os.path.isdir(dir_out):
            raise
    dir_fn = dir_out + fn_out
    df_all.to_csv(dir_fn, sep=',', columns=clmns_iea,
                  na_rep='nan', index=False)

    return dir_fn, df_all

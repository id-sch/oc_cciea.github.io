import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator


def fun_pd_df2IEA_fig_threshold_blue(df, nr, nc, order_list, yr_clim_bgn, yr_clim_end, wndw, yy_end, threshold, marker_flag=1):
    '''
    Function to plot a pandas dataframe (pd.df) that has
    been created by reading an "IEA CSV" file.    # input
    Inputs:
       nr = num rows of the figure
       nc = num columns of the figure
       order_list = which "IEA order" to plot
       yr_clim_bgn = start year, will be set to be xlim[0]
       yr_clim_end = end year
       wndw = the "IEA window" size in years
       yy_end = the "IEA" year, ie the current year
       threshold = the y-value where a horizontal line will be drawn
    Example:
       fun_pd_df2IEA_fig(4,1,[1,2],1984, 2017, 5, 2017, 1)
      plots 2 time series in a figure of 1 row, 4 columns
      and the xlim will be from 1984 to 2017, status and
      trends will be calculated over a 5 year window, draws
      horizontal threshold line at 1.
    '''
    # figure size
    fig_wdth = 8.5
    fig_hght = 11

    # IEA symbol boundary box
    # bbox_props = dict(boxstyle="circle,pad=0", fc="cyan", ec="cyan", lw=0)
    # year to begin window
    yy_bgn = yy_end - wndw + 1

    # --number of plots to make
    num_order_list = len(order_list)

    # --check to see if "Annual", "Monthly", or "All" IEA time series file
    clmns = df.columns.values
    in_ann = np.where(clmns == 'year')[0]
    in_month = np.where(clmns == 'month')[0]
    in_day = np.where(clmns == 'day')[0]
    if len(in_ann):
        time_flag = 1
    if len(in_month):
        time_flag = 2
    if len(in_day):
        time_flag = 3

    # --set up plot grid
    gs1 = gridspec.GridSpec(nr, nc)
    gs1.update(left=0.2, right=0.8, bottom=0.05,
               top=0.9,  wspace=0.3, hspace=0.4)

    # open new figure of set size
    plt.close()
    fig = plt.figure(figsize=(fig_wdth, fig_hght))

    for iii in range(0, num_order_list):
        print([iii, num_order_list])
        ax = plt.subplot(gs1[iii])
        in_order = df.order == order_list[iii]
        ttl_array = df.timeseries[in_order]
        ttl = ttl_array.values[0]
        ylbl_array = df.metric[in_order]
        ylbl = ylbl_array.values[0]

        # --get pd.df unique to the order number and create datetime index
        dfo = df[in_order].copy()
        yrsd = dfo.year.values
        if time_flag == 1:
            mond = np.ones(np.size(yrsd))
            dayd = np.ones(np.size(yrsd))
        if time_flag == 2:
            mond = dfo.month.values
            dayd = np.ones(np.size(yrsd))
        if time_flag == 3:
            mond = dfo.month.values
            dayd = dfo.day.values

        tt_dict = {'year': yrsd, 'month': mond, 'day': dayd}
        # dt = pd.to_datetime(tt_dict)
        dt = pd.DatetimeIndex(pd.to_datetime(tt_dict))
        dfo['Datetime'] = dt
        dfo = dfo.set_index('Datetime')

        # --subset by clim yers
        in_clim = ((dfo.index.year >= yr_clim_bgn) &
                   (dfo.index.year <= yr_clim_end))
        dfo = dfo[in_clim]

        # --some IEA time series have sd, plot if they do
        std_flag = np.isnan(dfo.error.mean())

        # dataf_sd1 = dfo.data.values - dfo.error.values
        # dataf_sd2 = dfo.data.values + dfo.error.values
        sd_wts = dfo.data.std(skipna=True)
        mn_wts = dfo.data.mean(skipna=True)

        # horizontal threshold line
        plt.axhline(threshold, color='blue')

        # --plot the data
        plt.plot(dfo.index, dfo.data, '-k', linewidth=1)
        if marker_flag == 1:
            plt.plot(dfo.index, dfo.data, '.', markersize=3, color='orange')

        # fill values less than threshold
        tt1 = dfo.index
        ts1 = dfo.data
        ind = np.isfinite(ts1)
        tt1d = tt1[ind]
        ts1d = ts1[ind]
        plt.fill_between(tt1d, threshold, ts1d, where=ts1d <= threshold,
                         color='skyblue', linewidth=0.1,
                         interpolate=True, alpha=0.7, zorder=30)

        # x-limits set to clim interval
        xlm1 = np.datetime64(str(yr_clim_bgn)+'-01')
        if time_flag > 1:
            xlm2 = np.datetime64(str(yr_clim_end)+'-12')
        else:
            xlm2 = np.datetime64(str(yr_clim_end)+'-03')
        plt.xlim([xlm1, xlm2])

        # --fill in the 5-year window in green
        x5yr = [np.datetime64(str(yy_bgn)+'-01'),
                np.datetime64(str(yy_end)+'-12')]
        y5yr1 = np.ones(2)*[mn_wts-sd_wts]
        y5yr2 = np.ones(2)*[mn_wts+sd_wts]
        plt.fill_between(x5yr, y5yr1, y5yr2, facecolor='dodgerblue', alpha=0.2)

        in_wndw = ((dfo.index.year >= yy_bgn) & (dfo.index.year <= yy_end))
        # nt_wndw = in_wndw.nonzero()[0].size
        # y_sd1 = np.ones(nt_wndw)*mn_wts - sd_wts
        # y_sd2 = np.ones(nt_wndw)*mn_wts + sd_wts
        # plt.fill_between(dfo.index[in_wndw], y_sd1, y_sd2, facecolor='palegreen', alpha=0.4)

        # --fill gray std window
        if ~std_flag:
            dataf_sd1 = dfo.data.values - dfo.error.values
            dataf_sd2 = dfo.data.values + dfo.error.values
            plt.fill_between(dfo.index, dataf_sd1, dataf_sd2,
                             facecolor=np.ones(3)*0.8, alpha=0.4)

        # --data over the 5-year window, remove missing
        ind = np.isfinite(dfo.data[in_wndw].values)
        tt5 = dfo.index[in_wndw][ind]
        ts5 = dfo.data[in_wndw][ind]

        # --get y-limits to position the IEA symbols
        ylm1 = ax.get_ylim()
        dy1 = (ylm1[1] - ylm1[0]) / 20
        ypt = mn_wts + dy1
        ypb = mn_wts - dy1

        # --IEA symbols, trend over last 5 years
        A1 = np.array([tt5.to_julian_date(), np.ones(ind.nonzero()[0].size)])
        w1 = np.linalg.lstsq(A1.T, ts5, rcond=-1)[0]
        trnd5 = w1[0] * tt5.to_julian_date() + w1[1]
        plt.plot(tt5, trnd5, '-r', linewidth=1)

        xp_sym = plt.xlim()[1] + np.diff(plt.xlim())*0.01
        dtrnd = trnd5[-1] - trnd5[0]
        if dtrnd > sd_wts:
            ax.text(xp_sym, ypt, r'$\nearrow$', fontweight='bold',
                    fontsize='large', verticalalignment='center')
        if dtrnd < -1 * sd_wts:
            ax.text(xp_sym, ypt, r'$\searrow$', fontweight='bold',
                    fontsize='large', verticalalignment='center')
        if dtrnd >= -1 * sd_wts and dtrnd <= sd_wts:
            ax.text(xp_sym, ypt, r'$\leftrightarrow$', fontweight='bold',
                    fontsize='large', verticalalignment='center')

        # --IEA symbols, last 5 years
        # bbox_props = dict(boxstyle="circle,pad=0", fc="cyan", ec="cyan", lw=2)
        mn5 = dfo.data[in_wndw][ind].mean() * np.ones(ind.nonzero()[0].size)
        plt.plot(tt5, mn5, '-m', linewidth=1)

        if mn5[0] > mn_wts + sd_wts:
            ax.text(xp_sym, ypb, r'$\plus$', fontweight='bold',
                    fontsize='large', verticalalignment='center')
        if mn5[0] < mn_wts - sd_wts:
            ax.text(xp_sym, ypb, r'$\minus$', fontweight='bold',
                    fontsize='large', verticalalignment='center')
        if mn5[0] >= mn_wts - sd_wts and mn5[0] <= mn_wts + sd_wts:
            ax.text(xp_sym, ypb, r'$\bullet$', fontweight='bold',
                    fontsize='large', verticalalignment='center')

        # --remove top and right axes lines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # make sure that only the year gets written as xticklabels
        monthyearFmt = mdates.DateFormatter('%Y')
        ax.xaxis.set_major_formatter(monthyearFmt)

        # set xticks, xticklabels at start and end of the window
        yr_all = yr_clim_end - yr_clim_bgn
        if yr_all < 38:
            dt_intrvl = -1*(wndw-1)
        else:
            dt_intrvl = -1*(2*wndw-1)
        xt1 = np.flip(np.arange(yr_clim_end, yr_clim_bgn, dt_intrvl), 0)
        xt1_dt = [np.datetime64(xt1[j].astype('str')+'-01')
                  for j in range(0, np.size(xt1))]

        if (iii == num_order_list-1):
            plt.xticks(xt1_dt)
            plt.xlabel('Year')
        else:
            plt.xticks(xt1_dt, '')
        minor_locator = AutoMinorLocator(4)
        ax.xaxis.set_minor_locator(minor_locator)

        # y-axis label and title
        ylbl_split = ylbl.split(' (')
        if len(ylbl_split) == 2:
            ylblf = '{}\n({}'.format(ylbl_split[0], ylbl_split[1])
        else:
            ylblf = ylbl
        plt.ylabel(ylblf)

        # title
        plt.title(ttl)

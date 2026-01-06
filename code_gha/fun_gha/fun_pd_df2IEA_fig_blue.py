import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib import interactive
from matplotlib.ticker import AutoMinorLocator
interactive(False)
# pylint: disable=C0103


def fun_pd_df2IEA_fig_blue(df, nr, nc, order_list, yr_clim_bgn, yr_clim_end,
                           wndw, yy_end, marker_flag=1,
                           fig_wdth=8.5, fig_hght=11):
    '''
    Function to plot a pandas dataframe (pd.df) that has
    been created by reading an "IEA CSV" file.    # input

    This creates figures with the new blue color theme.

    Inputs:
        nr = num rows of the figure
        nc = num columns of the figure
        order_list = which "IEA order" to plot
        yr_clim_bgn = start year, will be set to be xlim[0]
        yr_clim_end = end year
        wndw = the "IEA window" size in years
        yy_end = the "IEA" year, ie the current year
    Example:
        fun_pd_df2IEA_fig(4,1,[1,2],1984, 2017, 5, 2017)

        plots 2 time series in a figure of 1 row, 4 columns
        and the xlim will be from 1984 to 2017, status and
        trends will be calculated over a 5 year window.
    '''

    # year to begin window
    yy_bgn = yy_end - wndw + 1

    # --number of plots to make
    num_order_list = len(order_list)

    # --check to see if "Annual", "Monthly", or "All" IEA time series file
    clmns = df.columns.values
    in_ann = np.where(clmns == 'year')[0]
    in_month = np.where(clmns == 'month')[0]
    in_day = np.where(clmns == 'day')[0]

    in_data_or_index = np.logical_or(
        clmns == 'index', clmns == 'data').nonzero()[0]

    data_clmn_lbl = clmns[in_data_or_index]

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
    plt.figure(figsize=(fig_wdth, fig_hght))

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
        dt = pd.DatetimeIndex(pd.to_datetime(tt_dict))
        dfo['Datetime'] = dt
        dfo = dfo.set_index('Datetime')

        # --subset by clim years
        in_clim = ((dfo.index.year >= yr_clim_bgn) &
                   (dfo.index.year <= yr_clim_end))
        dfo = dfo[in_clim]

        # --some IEA time series have sd, plot if they do
        std_flag = np.isnan(dfo.error.mean())
        sd_wts = dfo[data_clmn_lbl].std(skipna=True).values[0]
        mn_wts = dfo[data_clmn_lbl].mean(skipna=True).values[0]

        # set xlim begin and end
        xlm1 = np.datetime64(str(yr_clim_bgn)+'-01')

        if time_flag > 1:
            xlm2 = np.datetime64(str(yr_clim_end)+'-12')
        else:
            xlm2 = np.datetime64(str(yr_clim_end)+'-03')

        # --horizontal lines at 0, plus/minus 1 std
        x_hl = np.arange(xlm1, xlm2)
        ones_hl = np.ones(len(x_hl))

        plt.plot(x_hl, ones_hl * mn_wts,
                 linestyle='--', color='k', linewidth=0.5)
        plt.plot(x_hl, ones_hl * (mn_wts - sd_wts),
                 linestyle='-', color='dodgerblue', linewidth=0.5)
        plt.plot(x_hl, ones_hl * (mn_wts + sd_wts),
                 linestyle='-', color='dodgerblue', linewidth=0.5)

        # --plot the data
        plt.plot(dfo.index, dfo[data_clmn_lbl], '-k', linewidth=1)
        if marker_flag == 1:
            plt.plot(dfo.index, dfo[data_clmn_lbl],
                     '.', markersize=3, color='orange')

        # set xticks, xticklabels at start and end of the window
        yr_all = yr_clim_end - yr_clim_bgn
        if yr_all < 38:
            dt_intrvl = -1*(wndw-1)
            num_minor = 4
        else:
            dt_intrvl = -1*(2*wndw)
            num_minor = 10
        xt1 = np.flip(
            np.arange(yr_clim_end, yr_clim_bgn+dt_intrvl, dt_intrvl), 0)
        xt1_dt = [np.datetime64(xt1[j].astype('str')+'-01')
                  for j in range(0, np.size(xt1))]
        xt1_dt = [np.datetime64(xt1[j].astype('str'))
                  for j in range(0, np.size(xt1))]

        if iii == num_order_list-1:
            plt.xticks(xt1_dt)
            plt.xlabel('Year')
            myFmt = DateFormatter("%Y")
            ax.xaxis.set_major_formatter(myFmt)
        else:
            plt.xticks(xt1_dt, '')
        minor_locator = AutoMinorLocator(num_minor)
        ax.xaxis.set_minor_locator(minor_locator)

        # x-limits set to clim interval
        plt.xlim([xlm1, xlm2])

        # --fill in the 5-year window in green
        x5yr = [np.datetime64(str(yy_bgn)+'-01'),
                np.datetime64(str(yy_end)+'-12')]
        y5yr1 = np.ones(2)*[mn_wts-sd_wts]
        y5yr2 = np.ones(2)*[mn_wts+sd_wts]
        plt.fill_between(x5yr, y5yr1, y5yr2, facecolor='dodgerblue', alpha=0.2)

        # # --fill gray std window
        # if ~std_flag:
        #     dataf_sd1 = dfo[data_clmn_lbl].values.squeeze() - \
        #         dfo['error'].values
        #     dataf_sd2 = dfo[data_clmn_lbl].values.squeeze() + \
        #         dfo['error'].values
        #     plt.fill_between(dfo.index.values, dataf_sd1, dataf_sd2,
        #                      facecolor=np.ones(3)*0.8, alpha=0.6)

        # # --data over the 5-year window, remove missing
        # in_wndw = ((dfo.year.values >= yy_bgn) & (dfo.year.values <= yy_end))
        # tt5_wndw = dfo.index[in_wndw]
        # ts5_wndw = dfo[data_clmn_lbl][in_wndw]

        # ind = np.isfinite(ts5_wndw).values.nonzero()[0]
        # tt5 = tt5_wndw[ind]
        # ts5 = ts5_wndw.values[ind]

        # # --get y-limits to position the IEA symbols
        # ylm1 = ax.get_ylim()
        # dy1 = (ylm1[1] - ylm1[0]) / 20
        # ypt = mn_wts + dy1
        # ypb = mn_wts - dy1

        # # --IEA symbols, trend over last 5 years
        # A1 = np.array([tt5.to_julian_date().values,
        #                np.ones(len(tt5))])
        # w1, w2 = np.linalg.lstsq(A1.T, ts5, rcond=None)[0]
        # trnd5 = w1 * tt5.to_julian_date().values + w2
        # plt.plot(tt5, trnd5, '-r', linewidth=1)

        # xp_sym = plt.xlim()[1] + np.diff(plt.xlim())*0.01
        # dtrnd = trnd5[-1] - trnd5[0]

        # if dtrnd > sd_wts:
        #     ax.text(xp_sym, ypt, r'$\nearrow$', fontweight='bold',
        #             fontsize='large', verticalalignment='center')
        # if dtrnd < -1 * sd_wts:
        #     ax.text(xp_sym, ypt, r'$\searrow$', fontweight='bold',
        #             fontsize='large', verticalalignment='center')
        # if dtrnd >= -1 * sd_wts and dtrnd <= sd_wts:
        #     ax.text(xp_sym, ypt, r'$\leftrightarrow$', fontweight='bold',
        #             fontsize='large', verticalalignment='center')

        # # --IEA symbols, last 5 years
        # mn5 = ts5.mean() * np.ones(len(tt5))
        # plt.plot(tt5, mn5, '-m', linewidth=1)

        # if mn5[0] > mn_wts + sd_wts:
        #     ax.text(xp_sym, ypb, r'$\plus$', fontweight='bold',
        #             fontsize='large', verticalalignment='center')
        # if mn5[0] < mn_wts - sd_wts:
        #     ax.text(xp_sym, ypb, r'$\minus$', fontweight='bold',
        #             fontsize='large', verticalalignment='center')
        # if mn5[0] >= mn_wts - sd_wts and mn5[0] <= mn_wts + sd_wts:
        #     ax.text(xp_sym, ypb, r'$\bullet$', fontweight='bold',
        #             fontsize='large', verticalalignment='center')

        # # --remove top and right axes lines
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)

        # # y-axis label and title
        # ylbl_split = ylbl.split(' (')
        # if len(ylbl_split) == 2:
        #     ylblf = '{}\n({}'.format(ylbl_split[0], ylbl_split[1])
        # else:
        #     ylblf = ylbl
        # plt.ylabel(ylblf)

        # # title
        # plt.title(ttl)

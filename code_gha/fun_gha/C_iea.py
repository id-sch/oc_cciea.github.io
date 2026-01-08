import numpy as np
import matplotlib.pyplot as plt
from matplotlib import interactive
from matplotlib.ticker import AutoMinorLocator
interactive(True)


# NOTE: M-x pyvenv-workon py_cart


class C_iea(object):
    """An IEA timeseries:

    Attributes (defined in __init__):
    yr_clim_bgn: single value (clim start year for the whole time series)
    yr_clim_end: single value (clim end year for the whole time series - wts)
    wndw: single value (number of years for the IEA trend and mean, default 5 years)
    yy_bgn: single value (window start year)
    yy_end: single value (window end year)
    tt5: vector (dates of all values in the wndw)
    ts5: vector (data in the wndw)
    mn5: single value (mean of all data in the wndw)
    trnd5: single value (trend in the wndw)
    anom_end: single value (the last non-NaN value in ts5)
    marker_anom: single value (' ' if less than 1 sd, '.' greater than 1 sd, '+' max/min value of wts)
    marker_mn5: single value (' ' if less than 1sd, '.' greater than 1 sd)
    marker_trnd5: single value (' ' if less than 1sd, '.' greater than 1 sd)

    fs_tcks: sometimes fontsize for tick labels are too small, adjust with this
    fs_mrkr: sometimes symbols are too small, use this to adjust the size

    """

    def __init__(self, ps1, yr_clim_bgn=np.nan, yr_clim_end=np.nan,
                 wndw=np.nan, yy_end=np.nan,
                 lon=0, lat=0, line=0, sttn=0, fs_mrkr=0, fs_tcks=0, name='index_name', units='units'):
        """Return a C_sw object:
           must have pandas Series, with a datetime index
        """

        # name
        self.name = name

        # units
        self.units = units

        # set lon and lat
        self.lat = lat
        self.lon = lon

        # set line and station
        self.line = line
        self.sttn = sttn

        # fontsize
        self.fs_mrkr = fs_mrkr
        self.fs_tcks = fs_tcks

        # if values of inputs are np.nan then infer the values
        # from "standard" IEA figures (eg window size of 5 years)
        if np.isnan(yr_clim_bgn):
            self.yr_clim_bgn = ps1.index[0].year
        else:
            self.yr_clim_bgn = yr_clim_bgn
        if np.isnan(yr_clim_end):
            self.yr_clim_end = ps1.index[-1].year
        else:
            self.yr_clim_end = yr_clim_end
        if np.isnan(wndw):
            self.wndw = 5
        else:
            self.wndw = wndw
        if np.isnan(yy_end):
            self.yy_end = self.yr_clim_end
        else:
            self.yy_end = yy_end

        # year to begin window
        self.yy_bgn = self.yy_end - self.wndw + 1
        yrs5 = np.arange(self.yy_bgn, yy_end + 1)

        # check the frequency of ps1, should be either M (monthly) or Annual (A)
        freq_str = ps1.index.freqstr[0]
        if freq_str == 'M':
            dt = 12
        elif freq_str == 'Y':
            dt = 1

        nt5 = len(yrs5)*dt

        # --subset by clim years
        in_clim = ((ps1.index.year >= self.yr_clim_bgn) &
                   (ps1.index.year <= self.yr_clim_end))
        self.ps2 = ps1[in_clim]

        # mean and standard deviation of the whole time series (wts)
        # self.mn_wts = ps1.mean(skipna=True)
        # self.sd_wts = ps1.std(skipna=True)
        self.mn_wts = self.ps2.mean(skipna=True)
        self.sd_wts = self.ps2.std(skipna=True)

        # data
        # self.data = self.ps2
        self.data = ps1

        # anom
        # self.anom = self.ps2 - self.mn_wts
        self.anom = self.data - self.mn_wts

        # data last year
        in_ly = np.where(ps1.index.year.values == yr_clim_end)[0]
        data_ly = ps1.values[in_ly]
        anom_ly = data_ly - self.mn_wts

        if len(data_ly) > 0:
            data_ly_mn = np.nanmean(data_ly)
            anom_ly_mn = np.nanmean(anom_ly)
        else:
            data_ly_mn = np.nan
            anom_ly_mn = np.nan

        # index of data in the window
        # in_wndw = ((self.ps2.index.year >= self.yy_bgn) &
        #            (self.ps2.index.year <= self.yy_end))
        in_wndw = ((self.data.index.year >= self.yy_bgn) &
                   (self.data.index.year <= self.yy_end))

        # --data over the 5-year window, remove missing
        # ind = np.isfinite(self.ps2.values[in_wndw])
        # self.tt5 = self.ps2.index[in_wndw][ind]
        # self.ts5 = self.ps2.values[in_wndw][ind]
        ind = np.isfinite(self.data.values[in_wndw])
        self.tt5 = self.data.index[in_wndw][ind]
        self.ts5 = self.data.values[in_wndw][ind]

        # mean and trend over the window
        if len(ind.nonzero()[0]) >= 0.6*nt5:
            # 1) mean
            # self.mn5 = self.ps2.values[in_wndw][ind].mean()
            self.mn5 = self.data.values[in_wndw][ind].mean()

            # 2) trend
            A1 = np.array([self.tt5.to_julian_date(),
                           np.ones(ind.nonzero()[0].size)])
            w1 = np.linalg.lstsq(A1.T, self.ts5, rcond=-1)[0]
            self.fit5 = w1[0] * self.tt5.to_julian_date() + w1[1]
            self.trnd5 = self.fit5[-1] - self.fit5[0]
        else:
            self.mn5 = np.nan
            self.trnd5 = np.nan

        # IEA MAP SYMBOLS
        # 1) Data or Anom symbols
        # self.data_end = self.data[-1]
        # self.anom_end = self.anom[-1]
        self.data_end = self.data.iloc[-1]
        self.anom_end = self.anom.iloc[-1]

        # less/greater than sd of the whole time series
        if abs(self.anom_end) >= self.sd_wts:
            self.map_marker_anom = '.'
        else:
            self.map_marker_anom = ' '

        # min/max value of the whole time series
        if self.anom_end == np.min(self.anom) or self.anom_end == np.max(self.anom):
            self.map_marker_anom = '+'

        # 2) 5-year mean symbols
        # anom over the 5-year window normalized by long-term std
        if self.sd_wts != 0:
            self.mn5_sd = (self.mn5 - self.mn_wts)/self.sd_wts
        else:
            self.mn5_sd = 0

        if abs(self.mn5_sd) > 1:
            self.map_marker_mn5 = '.'
        else:
            self.map_marker_mn5 = ' '

        # 3) 5-year trend symbols
        # trend over the 5-year window normalized by long-term std
        if self.sd_wts != 0:
            self.trnd5_sd = self.trnd5/self.sd_wts
        else:
            self.trnd5_sd = 0

        if abs(self.trnd5_sd) > 1:
            self.map_marker_trnd5 = '.'
        else:
            self.map_marker_trnd5 = ' '

        # IEA TIME SERIES SYMBOLS
        # 1) 5-year mean symbols
        # anom over the 5-year window normalized by long-term std
        if self.sd_wts != 0:
            self.mn5_sd = (self.mn5 - self.mn_wts)/self.sd_wts
        else:
            self.mn5_sd = 0

        if self.mn5_sd > 1:
            self.ts_marker_mn5 = r'$\plus$'
        if self.mn5_sd < -1:
            self.ts_marker_mn5 = r'$\minus$'
        if self.mn5_sd >= -1 and self.mn5_sd <= 1:
            self.ts_marker_mn5 = r'$\bullet$'

        # 2) 5-year trend symbols
        # trend over the 5-year window normalized by long-term std
        if self.sd_wts != 0:
            self.trnd5_sd = self.trnd5/self.sd_wts
        else:
            self.trnd5_sd = 0

        if self.trnd5_sd > 1:
            self.ts_marker_trnd5 = r'$\nearrow$'
        if self.trnd5_sd < -1:
            self.ts_marker_trnd5 = r'$\searrow$'
        if self.trnd5_sd >= -1 and self.trnd5_sd <= 1:
            self.ts_marker_trnd5 = r'$\leftrightarrow$'

        # create dictionary of status and trends
        data_dict = {'info': {
            'name': self.name,
            'units': self.units,
            'window_length': wndw,
            'window_year': [self.yy_bgn, self.yy_end],
            'clim_length': yr_clim_end - yr_clim_bgn + 1,
            'clim_year': [yr_clim_bgn, yr_clim_end],
            'time_range': ['{}-{:02d}'.format(ps1.index.year.values[0], ps1.index.month.values[0]),
                           '{}-{:02d}'.format(ps1.index.year.values[-1], ps1.index.month.values[-1])]
        }}

        data_dict['monthly'] = {
            'data': self.data.iloc[-1],
            'anom': self.anom.iloc[-1],
            'last_yr_mean': data_ly_mn,
            'last_yr_anom': anom_ly_mn
        }

        self.data_dict = data_dict

    def plot(self):
        # check the frequency of ps1, should be either M (monthly) or Annual (A)
        freq_str = self.ps2.index.freqstr[0]
        if freq_str == 'M':
            time_flag = 2
        elif freq_str == 'Y':
            time_flag = 1

        # horizontal lines at 0, plus/minus 1 std
        plt.plot(self.ps2.index,
                 np.ones(self.ps2.shape[0]) * self.mn_wts,
                 linestyle='--', color='k', linewidth=0.5)
        plt.plot(self.ps2.index,
                 np.ones(self.ps2.shape[0]) * (self.mn_wts - self.sd_wts),
                 linestyle='-', color='g', linewidth=0.5)
        plt.plot(self.ps2.index,
                 np.ones(self.ps2.shape[0]) * (self.mn_wts + self.sd_wts),
                 linestyle='-', color='g', linewidth=0.5)

        # plot the data
        plt.plot(self.ps2.index, self.ps2.values, '-r', linewidth=1)
        plt.plot(self.ps2.index, self.ps2.values,
                 '.k', markersize=3, clip_on=False)

        # x-limits set to clim interval
        xlm1 = np.datetime64(str(self.yr_clim_bgn)+'-01')
        if time_flag > 1:
            xlm2 = np.datetime64(str(self.yr_clim_end)+'-12')
        else:
            xlm2 = np.datetime64(str(self.yr_clim_end)+'-03')
        plt.xlim([xlm1, xlm2])

        # fill in the 5-year window in green
        x5yr = [np.datetime64(str(self.yy_bgn)+'-01'),
                np.datetime64(str(self.yy_end)+'-12')]
        y5yr1 = np.ones(2)*[self.mn_wts-self.sd_wts]
        y5yr2 = np.ones(2)*[self.mn_wts+self.sd_wts]
        plt.fill_between(x5yr, y5yr1, y5yr2, facecolor='palegreen', alpha=0.4)

        # get x position of the IEA symbols
        xp_sym = plt.xlim()[1] + np.diff(plt.xlim())*0.04
        xpt = xp_sym[0]

        # get y positions of the IEA symbols
        ylm1 = plt.gca().get_ylim()
        dy1 = (ylm1[1] - ylm1[0]) / 12
        ypt = self.mn_wts + dy1
        ypb = self.mn_wts - dy1

        # mean line and symbol
        mn5 = np.ones(len(self.tt5))*self.mn5
        plt.plot(self.tt5, mn5, '-m', linewidth=1)

        if self.fs_mrkr > 0:
            plt.gca().text(xpt, ypb, self.ts_marker_mn5, fontweight='bold',
                           fontsize=self.fs_mrkr, verticalalignment='center')
        else:
            plt.gca().text(xpt, ypb, self.ts_marker_mn5, fontweight='bold',
                           fontsize='smaller', verticalalignment='center')

        # trend line and symbol
        plt.plot(self.tt5, self.fit5, '-g', linewidth=1)

        if self.fs_mrkr > 0:
            plt.gca().text(xpt, ypt, self.ts_marker_trnd5, fontweight='bold',
                           fontsize=self.fs_mrkr, verticalalignment='center')
        else:
            plt.gca().text(xpt, ypt, self.ts_marker_trnd5, fontweight='bold',
                           fontsize='smaller', verticalalignment='center')

        # --remove top and right axes lines
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)

        # set the number of yticks
        plt.locator_params(axis='y', nbins=6)

        # set xticks, xticklabels at start and end of the window
        yr_all = self.yr_clim_end - self.yr_clim_bgn
        if yr_all < 28:
            dt_intrvl = -1*(self.wndw-1)
            xtck_mnr = 4
        else:
            dt_intrvl = -1*(2*self.wndw)
            xtck_mnr = 10

        xt1 = np.flip(np.arange(self.yr_clim_end,
                                self.yr_clim_bgn, dt_intrvl), 0)
        xt1_dt = [np.datetime64(xt1[j].astype('str')+'-01')
                  for j in range(0, np.size(xt1))]

        if self.fs_tcks > 0:
            plt.xticks(xt1_dt, fontsize=self.fs_tcks, labels=xt1)
            plt.yticks(fontsize=self.fs_tcks)
        else:
            plt.xticks(xt1_dt)

        minor_locator = AutoMinorLocator(xtck_mnr)
        plt.gca().xaxis.set_minor_locator(minor_locator)

        # title
        if self.lon > 180:
            lon1 = self.lon - 360
        else:
            lon1 = self.lon

        ttl = 'Lat {:.2f} Lon {:.2f}, Line {:.1f} Sttn {:.1f}'.format(
            lon1, self.lat, self.line, self.sttn)
        plt.title(ttl)

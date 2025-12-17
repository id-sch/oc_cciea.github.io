# ---------------------------------------------------------------------------------------
# --Hampel is a simple filter to remove outliers, based on window size around
# --a sample point and how it deviates from the sd of the points in the window
# --k = window size
# --t0 = number of standard deviations to exceed inorder to test for outlier
# --iflag = check to see if only one value, if there is then set all values
#           to NaN. This is useful for interpolation which needs 2 or more data
# --see:
# --https://www.mathworks.com/help/signal/ref/hampel.html?requestedDomain=www.mathworks.com
# ---------------------------------------------------------------------------------------
import numpy as np


def hampel(vals_orig, k=7, t0=3, iflag=0):
    '''
    vals: pandas series of values from which to remove outliers
    k: size of window (including the sample; 7 is equal to 3
    on either side of value)
    '''
    # Make copy so original not edited
    vals = vals_orig.copy()
    # Hampel Filter
    L = 1.4826
    rolling_median = vals.rolling(k).median()
    difference = np.abs(rolling_median-vals)
    median_abs_deviation = difference.rolling(k).median()
    threshold = t0 * L * median_abs_deviation
    outlier_idx = difference > threshold
    vals[outlier_idx] = np.nan

    if iflag == 1:
        if vals.count() == 1:
            vals[:] = np.nan

    return(vals)

import os
import requests
import numpy as np
from datetime import date



def fun_get_ts_erddap(yr_wnt1, yr_wnt2, mon1, mon2, dataset_ID, name1, fn1_out, file1_type):

    time_bgn = '{}-{:02d}-01'.format(yr_wnt1, mon1)
    time_end = '{}-{:02d}-01'.format(yr_wnt2, mon2)

    dataset_ID_lat = dataset_ID.format(name1)

    uri = (
        'https://oceanview.pfeg.noaa.gov/erddap/tabledap/{}.{}'
        '?time,{},station_id,latitude,longitude'
        '&time>={}'
        '&time<={}').format
    url = uri(dataset_ID_lat, file1_type, name1, time_bgn, time_end)

    # urllib only excepts ansi characters, has problems with commas, greater/less than symbols
    # use HTML URL encode sequences
    url = url.replace(',', '%2C')
    url = url.replace('>', '%3E')
    url = url.replace('<', '%3C')

    ddd = requests.get(url, timeout=500)

    open(fn1_out, 'wb').write(ddd.content)

    return url


# -------------------------------------------------------
# -- Start: Input variables, change these
lat_wnt = np.arange(24, 60, 3)

dataset_ID = 'erdUI{}6hr'
var_wnt = 'upwelling_index'

today = date.today()

mon_wnt1 = 1
mon_wnt2 = today.month

yr_bgn = 1967
yr_end = today.year

file_type = 'nc'
dir_out = './data_gha/bakunUI/'

# -- End: Input variables, change these
# -------------------------------------------------------

# size of input variables
num_lat_wnt = len(lat_wnt)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_out)
except OSError:
    if not os.path.isdir(dir_out):
        raise

# loop over all erddap data wanted
for i in range(0, num_lat_wnt):
    fn_out = dir_out + 'ui_{}_{}_{:02d}.{}'.format(lat_wnt[i], yr_end, mon_wnt2, file_type)

    # -------------------------------------------------------
    # --now, call function to download the data from ERDDAP
    # File path
    a = fn_out

    # Check if the file exists and is a file
    if os.path.isfile(fn_out):
        print("File exists, if you want to download then delete fn_out and run code again.")
    else:
        # remove old files
        files_old = os.listdir(dir_out)

        for j in range(len(files_old)):
            file_pre = 'ui_{}'.format(lat_wnt[i])
            if files_old[j].startswith(file_pre):
                os.remove('./{}/{}'.format(dir_out, files_old[j]))
        
        # download the data
        dataset_ID1 = dataset_ID.format(lat_wnt[i])
        url_erddap = fun_get_ts_erddap(yr_bgn, yr_end, mon_wnt1, mon_wnt2, dataset_ID1, var_wnt, fn_out, file_type)

        print(url_erddap)

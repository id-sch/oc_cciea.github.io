import os
import re
import urllib.request  # py36
import requests  # Aug 9, 2024 -- changed to request as urlib gave 404 error
from shutil import copyfile
from datetime import date
# pylint: disable=C0103



def fun_get_ts_erddap(yr_wnt1, yr_wnt2, mon1, mon2, name1, fn1_out, file1_type):

    time_bgn = '{}-{:02d}-01'.format(yr_wnt1, mon1)
    time_end = '{}-{:02d}-01'.format(yr_wnt2, mon2)

    uri = (
        'https://oceanview.pfeg.noaa.gov/erddap/tabledap/cciea_OC_{}.{}'
        '?time,{}'
        '&time>={}'
        '&time<={}').format
    url = uri(name1, file1_type, name1, time_bgn, time_end)

    # urllib only excepts ansi characters, has problems with commas, greater/less than symbols
    # use HTML URL encode sequences
    url = url.replace(',', '%2C')
    url = url.replace('>', '%3E')
    url = url.replace('<', '%3C')

    ddd = requests.get(url, timeout=20)

    open(fn1_out, 'wb').write(ddd.content)

    return url


# -------------------------------------------------------
# -- Start: Input variables, change these
name_wnt = ['ONI', 'PDO', 'NPGO']

today = date.today()

mon_wnt1 = 1
mon_wnt2 = today.month

yr_bgn = 1900
yr_end = today.year

file_type = 'nc'


# -- End: Input variables, change these
# -------------------------------------------------------

# size of input variables
num_wnt = len(name_wnt)

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_out)
except OSError:
    if not os.path.isdir(dir_out):
        raise

# loop over all erddap data wanted
for i in range(0, num_wnt):
    fn_out = dir_out + 'ts_{}_{}_{:02d}.{}'.format(name_wnt[i], yr_end, mon_wnt2, file_type)

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
            file_pre = 'ts_{}'.format(name_wnt[i])
            if files_old[j].startswith(file_pre):
                os.remove('./{}/{}'.format(dir_out, files_old[j]))
        
        # download the data
        url_erddap = fun_get_ts_erddap(yr_bgn, yr_end, mon_wnt1, mon_wnt2, name_wnt[i], fn_out, file_type)

        print(url_erddap)

import os
import subprocess
import numpy as np


# -------------------------------------------------------
# -- Start: Input variables, change these
url1 = 'https://mjacox.com/wp-content/uploads/CUTI_daily.nc'
url2 = 'https://mjacox.com/wp-content/uploads/BEUTI_daily.nc'

var1_wnt = 'CUTI'
var2_wnt = 'BEUTI'

file_type = 'nc'
dir_out = './data_gha/newUI/'
# dir_out = './data_x13/newUI'


# -- End: Input variables, change these
# -------------------------------------------------------

# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_out)
except OSError:
    if not os.path.isdir(dir_out):
        raise


url_wnt = [url1, url2]
var_wnt = [var1_wnt, var2_wnt]

num_var = len(var_wnt)

for i in range(num_var):
    file_out = '{}/{}_daily_mgj.nc'.format(dir_out, var_wnt[i])
    subprocess.run(["wget", "-O", file_out, url_wnt[i]])

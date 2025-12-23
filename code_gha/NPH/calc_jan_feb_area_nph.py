import xarray as xr
import pandas as pd
import numpy as np



# function to find an interval of months, use this to get
# intervals such as March-May
def is_mon_int(month, mon1, mon2):
    return (month >= mon1) & (month <= mon2)

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# monthly means
mon_wnt = [1, 2]

# NPH has very large values, multiply by scale factor to make them smaller
sf = 1.0/np.power(10, 6)

dir_out = './data_gha/NPH/'
dir_in = dir_out
# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

# read CSV file
fn_csv = '{}year_mon_area_max_x_y_lon_lat.csv'.format(dir_in)
df1 = pd.read_csv(fn_csv, header=None)

# add columns to pd.df, call the area by the name 'area_large' becuase it is large number
# and the final 'area' name will be scaled by the variable sf to make it a much
# smaller value
df2 = pd.DataFrame(df1.values, columns=[
                   'year', 'month', 'area_large', 'max', 'x', 'y', 'lon', 'lat'])

# create time and set as index
dt1 = pd.to_datetime(
    {'year': df2['year'].values, 'month': df2['month'].values, 'day': np.ones(df2.shape[0])})
df2['time'] = dt1
df3 = df2.set_index('time')


# multiply area by scale factor
nc = len(df3.columns)
df3.insert(loc=nc, column='area', value=sf*df3['area_large'].values)


# save to xarray
da1 = df3.to_xarray()

# save to file
fn_nph = '{}NPH.nc'.format(dir_out)
da1.to_netcdf(fn_nph)

# get month intervals wanted and take mean
ds1_mons = da1.sel(time=is_mon_int(da1['time.month'], mon_wnt[0], mon_wnt[1]))
ds1_mn = ds1_mons.groupby('time.year').mean('time')
ds1_sd = ds1_mons.groupby('time.year').std('time')

# only save area, put these into data array
da1_mons = da1['area']
da1_mn = ds1_mn['area']
da1_sd = ds1_sd['area']

# now create final xr dataset
ds1_out = da1_mons.to_dataset(name='Area Monthly')
ds2_out = da1_mn.to_dataset(name='Area Jan-Feb mean')
ds2_out['Area Jan-Feb sd'] = da1_sd


# save to netcdf
fn_mons = '{}nph_area_month.nc'.format(dir_out)
fn_mn = '{}nph_area_jan_feb.nc'.format(dir_out)

ds1_out.to_netcdf(fn_mons)
ds2_out.to_netcdf(fn_mn)

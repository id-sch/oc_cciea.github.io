import os
import xarray as xr
import numpy as np
import seawater as sw

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------

# ROMS mask file (I uploaded this to GitHUB on 7/1/26)
fn_mask = './data_gha/HCI_ROMS/mask_rho.nc'
# fn_mask = './oc_cciea.github.io/data_gha/HCI_ROMS/mask_rho.nc'

land_mask = 0
water_mask = 1

# x limit
xlm1 = -134
xlm2 = -115

# y limit
ylm1 = 30
ylm2 = 48

var_roms = ['lon_rho', 'lat_rho', 'mask_rho']

# bath file name
fn_bath = './data_gha/HCI_ROMS/bath.nc'
# fn_bath = './oc_cciea.github.io/data_gha/HCI/bath.nc'

# distance
dis_wnt = [75, 100, 150, 175, 200, 225, 250, 275, 300, 325, 350]

# dir_out
# dir_out = './data_x13/HCI_ROMS/'
dir_out = './data_gha/HCI_ROMS/'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

dir_list = os.listdir()
print("START -------------------------------")
print("Files and directories in  :")
print(dir_list)

# directory to save netcdf output
# check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_out)
except OSError:
    if not os.path.isdir(dir_out):
        raise

dis_str = '_'.join(list(map(str,  dis_wnt)))
fn_out = '{}roms_distance_to_shore_dis_{}km.nc'.format(dir_out, dis_str)

# only calculate the mask grid if fn_out doesn't exist
if os.path.isfile(fn_out):
#    print("Land Masks for {} km have already been calculated".format(dis_str))
#else:
    print("Calculating Land Masks for {} km".format(dis_str))


    # open Bath netcdf file as an xarray
    ds1_bath = xr.open_dataset(fn_bath)
    var_ds1_bath = list(ds1_bath.keys())
    coord_ds1_bath = list(ds1_bath.coords.keys())

    lon_bath = ds1_bath[coord_ds1_bath[1]].data
    lat_bath = ds1_bath[coord_ds1_bath[0]].data
    topo = ds1_bath[var_ds1_bath[0]].data

    da1_bath = ds1_bath[var_ds1_bath[0]]

    # ROMS mask
    ds_roms = xr.open_dataset(fn_mask)

    lon_vec = ds_roms[var_roms[0]].data
    lat_vec = ds_roms[var_roms[1]].data
    roms = np.squeeze(ds_roms[var_roms[2]])

    # create mask, land=0, water=1
    mask_roms = np.copy(roms)

    # meshgrid of lon_vec, lat_vec
    lon_roms, lat_roms = np.meshgrid(lon_vec, lat_vec)

    # subset by lat
    in_lat = np.logical_and(lat_roms[:, 0] >= ylm1, lat_roms[:, 0] <= ylm2)
    lon_roms1 = lon_roms[in_lat, :]
    lat_roms1 = lat_roms[in_lat, :]
    mask_roms1 = mask_roms[in_lat, :]

    # subset by lon
    in_lon = np.logical_and(lon_roms1[0, :] >= xlm1, lon_roms1[0, :] <= xlm2)
    lon_roms2 = lon_roms1[:, in_lon]
    lat_roms2 = lat_roms1[:, in_lon]
    mask_roms2 = mask_roms1[:, in_lon]

    da1 = xr.DataArray(mask_roms2, coords=[lat_roms2[:, 0], lon_roms2[0, :]], dims=['latitude', 'longitude'])
    ny, nx = da1.shape

    # lat and lon
    lat1 = da1['latitude'].data
    lon1 = da1['longitude'].data

    # size of grid cell in degress
    dcell = (lat1[1]*10 - lat1[0]*10)/10

    # find land edge
    xp_land = np.zeros([ny, 1])
    yp_land = np.zeros([ny, 1])
    lon_land = np.zeros([ny, 1])
    lat_land = np.zeros([ny, 1])
    for i in range(ny):
        # find the index of all water values
        in_wtr = np.where(da1[i,:]==water_mask)[0]

        # the land index will be the last water index + 1
        xp_land[i, 0] = in_wtr[-1] + 1
        yp_land[i, 0] = i

        # lon and lat of land
        lon_land[i, 0] = da1['longitude'].data[in_wtr[-1] + 1]
        lat_land[i, 0] = da1['latitude'].data[i]

    jpw, ipw = np.where(da1.data == water_mask)
    jpl, ipl = np.where(da1.data == land_mask)

    # for each grid point find the distance to all the shore locations (there are ny of them)
    dis_mtrx = np.zeros([ny, nx, ny])*np.nan
    for i in range(ny):
        print('i: {}, ny: {}'.format(i, ny))
        lati = da1['latitude'].data[i]
        for j in range(nx):
            loni = da1['longitude'].data[j]

            # calc distance to all land values if the grid point is water
            if da1.data[i, j] == water_mask:
                for k in range(ny):
                    dis1 = sw.dist((lati, lat_land[k][0]), (loni, lon_land[k][0]), units='km')
                    dis1_km = dis1[0][0]
                    dis_mtrx[i, j, k] = dis1_km

    # now that each grid location has the distance to all shore locations, find the minimum
    # distance that is within the desired distance (set with dis_wnt array)
    mask_dis_shore = np.zeros([ny, nx, len(dis_wnt)])
    dis_shore = np.zeros([ny, nx, len(dis_wnt)])*np.nan
    depth_shore = np.zeros([ny, nx, len(dis_wnt)])*np.nan
    mask_area = np.zeros([ny, nx, len(dis_wnt)])*np.nan
    for i in range(ny):
        lati = da1['latitude'].data[i]

        for j in range(nx):
            loni = da1['longitude'].data[j]
            dis_ij = dis_mtrx[i, j, :]

            for k in range(len(dis_wnt)):
                # find shore locations within the distance wanted
                in_lt = np.where(dis_ij <= dis_wnt[k])[0]

                if len(in_lt) > 0:
                    # set mask to 1 if it is within the distance wanted
                    mask_dis_shore[i, j, k] = 1

                    # distance to shore is the minimim value of all the shore locations
                    dis_shore[i, j, k] = np.min(dis_ij[in_lt])

                    # find nearest lat and lon for the bath
                    da1_bath_ij = da1_bath.sel(
                        latitude=lati, method='nearest').sel(longitude=loni, method='nearest')

                    depth_shore[i, j, k] = da1_bath_ij.data

                    # calculate are of grid cell
                    lati_box = [lati-dcell/2, lati+dcell/2]
                    loni_box = [loni-dcell/2, loni+dcell/2]

                    dis_left = sw.dist(lati_box, [loni_box[0], loni_box[0]], units='km')
                    dis_right = sw.dist(lati_box, [loni_box[1], loni_box[1]], units='km')

                    dis_top = sw.dist([lati_box[1], lati_box[1]], loni_box, units='km')
                    dis_bot = sw.dist([lati_box[0], lati_box[0]], loni_box, units='km')
                    area_box = 0.5*(dis_top[0] + dis_bot[0])*dis_left[0]
                    mask_area[i, j, k] = area_box[0]

    # create dataarray
    da1_out = xr.DataArray(mask_dis_shore,
                           coords=[lat1, lon1, dis_wnt],
                           dims=['latitude', 'longitude', 'distance'])

    da2_out = xr.DataArray(dis_shore,
                           coords=[lat1, lon1, dis_wnt],
                           dims=['latitude', 'longitude', 'distance'])

    da3_out = xr.DataArray(depth_shore,
                           coords=[lat1, lon1, dis_wnt],
                           dims=['latitude', 'longitude', 'distance'])

    da4_out = xr.DataArray(mask_area,
                           coords=[lat1, lon1, dis_wnt],
                           dims=['latitude', 'longitude', 'distance'])

    # dataarray
    ds1_out = da1_out.to_dataset(name='mask_mtrx')
    ds1_out['distance_mtrx'] = da2_out
    ds1_out['depth_mtrx'] = da3_out
    ds1_out['area_mtrx'] = da4_out

    # save to netcdf
    # dis_str = '_'.join(list(map(str,  dis_wnt)))
    # fn_out = '{}roms_distance_to_shore_dis_{}km.nc'.format(dir_out, dis_str)
    ds1_out.to_netcdf(fn_out)


    dir_list_end = os.listdir()
    print("END -------------------------------")
    print("Files and directories in  :")
    print(dir_list_end)

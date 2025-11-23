import pandas as pd
import numpy as np


def fun_script_IDS(df1, cell2_wnt, ds_id, ts_id_list):
    '''
    Filters data from metadata dataframe.

    Input variables:
    1) df1 = dataframe of the metadata
    2) cell2_wnt = list of columns wanted in df1
    3) ds_id = dataset id
    4) ts_id_list = 
    Output:
    1) vec2 = dataframe with general information
    2) ds_id_lbl = erddap id
    3) rgn_lbl = lat/lon or region
    '''

    df2 = df1[cell2_wnt]

    data2 = []
    for j in range(len(cell2_wnt)):
        data2.append(df2[cell2_wnt[j]][0])

    vec2 = pd.DataFrame(np.expand_dims(np.array(data2),axis=0), columns=cell2_wnt)    

    df = pd.DataFrame(vec2, cell2_wnt)    

    num_ts_id_list = len(ts_id_list)

    rgn_lbl_vec = []

    for i in range(len(ts_id_list)):
        df3 = df1.loc[df1['CCIEA_timeseries_ID'] == ts_id_list[i]]
        rgn = df3['region'].values[0]
        lat = df3['latitude'].values[0]
        lon = df3['longitude'].values[0]
        lat2 = df3['latitude2'].values[0]
        lon2 = df3['longitude2'].values[0]
        if not np.isnan(lat):
            if np.isnan(lon):
                latlon = '{} ({:4.1f}N, {:4.1f}W)'.format(rgn, lat, lon)
            else:
                latlon = '{} ({:4.1f}N)'.format(rgn, lat)
            if not np.isnan(lat2):
                if not np.isnan(lon2):
                    latlon = '{} ({:4.1f}-{:4.1f}N, {:4.1f}-{:4.1f}W)'.format(rgn, lat, lat2, lon, lon2)
                else:
                    latlon = '{} ({:4.1f}-{:4.1f}N)'.format(rgn, lat, lat2)	    
        else:
            latlon = 'NA'
        rgn_lbl_vec.append(latlon)
    rgn_lbl = ', '.join(rgn_lbl_vec)
    ds_id_lbl = ', '.join(ts_id_list)

    return vec2, ds_id_lbl, rgn_lbl

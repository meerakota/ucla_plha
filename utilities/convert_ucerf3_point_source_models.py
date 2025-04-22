import numpy as np
import pandas as pd

def get_points(points_input_filename, weights_input_filename, points_output_filename, node_index_output_filename, ruptures_output_filename):
    '''
    Convert numpy array of points on surface of earth to Cartesian coordinates, returned as 
    Numpy array with same dimensions as inputs. We use only point source events with magnitude
    higher than 4.0 since smaller earthquakes are unlikely to cause liquefaction.
    '''
    columns = ['Node Index', 'Latitude', 'Longitude', '4.05', '4.15', '4.25', '4.35', '4.45', '4.55', '4.65', '4.75', '4.85',
       '4.95', '5.05', '5.15', '5.25', '5.35', '5.45', '5.55', '5.65', '5.75',
       '5.85', '5.95', '6.05', '6.15', '6.25', '6.35', '6.45', '6.55', '6.65',
       '6.75', '6.85', '6.95', '7.05', '7.15', '7.25', '7.35', '7.45', '7.55',
       '7.65', '7.75', '7.85', '7.95', '8.05', '8.15', '8.25', '8.35', '8.45',
       '8.55', '8.65', '8.75', '8.85', '8.95']
    df_points = pd.read_csv(points_input_filename, usecols=columns)
    node_index = df_points['Node Index'].values
    lat = df_points['Latitude'].values
    lon = df_points['Longitude'].values
    d = np.zeros(len(lat), dtype=float)
    points_xyz = np.empty((len(lat), 3), dtype=float)
    rad = np.pi/180.0
    a = 6378.1370 # Earth's equatorial radius in km
    b = 6356.7523 # Earth's polar radius in km
    r = np.sqrt(((a**2 * np.cos(lat * rad))**2 + (b**2 * np.sin(lat * rad))**2) / ((a * np.cos(lat * rad))**2 + (b * np.sin(lat * rad))**2))
    points_xyz[:, 0] = (
        (r - d)
        * np.cos(lat * rad)
        * np.cos(lon * rad)
    )
    points_xyz[:, 1] = (
        (r - d)
        * np.cos(lat * rad)
        * np.sin(lon * rad)
    )
    points_xyz[:, 2] = (r - d) * np.sin(lat * rad)
    np.save(points_output_filename, points_xyz)
    np.save(node_index_output_filename, node_index)

    df_weights = pd.read_csv(weights_input_filename)

    m = np.arange(4.05, 9.0, 0.1)
    node_index_all = np.repeat(node_index, len(m))
    rates = df_points.values.T[3:].T.reshape((len(node_index_all)))
    m_array = np.tile(m, len(df_points))
    unique_node_index = np.unique(node_index)
    weight_node_index = df_weights['Node Index'].values
    filter = []
    for wni in weight_node_index:
        if(wni in unique_node_index):
            filter.append(True)
        else:
            filter.append(False)
    rate_ss = np.repeat(df_weights['Fraction Strike-Slip'].values[filter], len(m)) * rates
    rate_rs = np.repeat(df_weights['Fraction Reverse'].values[filter], len(m)) * rates
    rate_ns = np.repeat(df_weights['Fraction Normal'].values[filter], len(m)) * rates
    rate = np.hstack((rate_ss, rate_rs, rate_ns))
    style_ss = np.full(len(rate_ss), 3)
    style_ns = np.full(len(rate_ns), 2)
    style_rs = np.full(len(rate_rs), 1)
    style = np.hstack((style_ss, style_ns, style_rs))
    df_out = pd.DataFrame()
    df_out['node_index'] = np.hstack((node_index_all, node_index_all, node_index_all))
    df_out['m'] = np.hstack((m_array, m_array, m_array))
    df_out['rate'] = rate
    df_out['style'] = style
    df_out[df_out['rate']> 0].to_pickle(ruptures_output_filename, compression='gzip')
    return

### Run function for UCERF3 input files
get_points(
    'FM3_1_branch_averaged/solution/grid_sub_seis_mfds.csv', 
    'FM3_1_branch_averaged/solution/grid_mech_weights.csv', 
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm31_grid_sub_seis/points.npy', 
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm31_grid_sub_seis/node_index.npy',
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm31_grid_sub_seis/ruptures.pkl'
)
get_points(
    'FM3_2_branch_averaged/solution/grid_sub_seis_mfds.csv', 
    'FM3_2_branch_averaged/solution/grid_mech_weights.csv', 
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm32_grid_sub_seis/points.npy', 
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm32_grid_sub_seis/node_index.npy',
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm32_grid_sub_seis/ruptures.pkl'
)
get_points(
    'FM3_1_branch_averaged/solution/grid_unassociated_mfds.csv', 
    'FM3_1_branch_averaged/solution/grid_mech_weights.csv', 
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm31_grid_unassociated/points.npy', 
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm31_grid_unassociated/node_index.npy',
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm31_grid_unassociated/ruptures.pkl'
)
get_points(
    'FM3_2_branch_averaged/solution/grid_unassociated_mfds.csv', 
    'FM3_2_branch_averaged/solution/grid_mech_weights.csv', 
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm32_grid_unassociated/points.npy', 
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm32_grid_unassociated/node_index.npy',
    '../src/ucla_plha/source_models/point_source_models/ucerf3_fm32_grid_unassociated/ruptures.pkl'
)
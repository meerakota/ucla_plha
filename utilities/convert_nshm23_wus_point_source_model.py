'''
Convert the 2023 USGS NSHM western U.S. branch-averaged gridded seismicity model into the
Numpy files used by ucla_plha. Output files have the same names and contents as the UCERF3
point source models produced by convert_ucerf3_point_source_models.py.

Input files (download into utilities/nshm23_wus_branch_avg/) from
https://code.usgs.gov/ghsc/nshmp/nshms/nshm-conus/-/tree/main/active-crust/grid
    grid-data/branch-avg-grid.csv: lon, lat, focal mechanism weights (ss_wt, r_wt, n_wt), and
        incremental annual rates for magnitude bins 5.05 to 7.85. This is the total gridded
        seismicity (sub-seismogenic plus unassociated), so there is a single NSHM23 point source model.
    features/grid-system-active.geojson: polygon defining the active crust region where the
        WUS gridded seismicity is applied. Nodes outside the polygon are in the stable crust and
        are modeled by the central and eastern U.S. sources, so they are excluded here.
'''
import json
import os
import numpy as np
import pandas as pd
import shapely

input_dir = 'nshm23_wus_branch_avg'
output_dir = '../src/ucla_plha/source_models/point_source_models/nshm23_wus_grid'

def get_points(points_input_filename, region_filename, points_output_filename, node_index_output_filename, ruptures_output_filename):
    '''
    Convert grid nodes on surface of earth to Cartesian coordinates, and save the magnitude, style of faulting,
    and rate of each point source event with non-zero rate.
    '''
    df_points = pd.read_csv(points_input_filename)
    region = shapely.geometry.shape(json.load(open(region_filename))['geometry'])
    df_points = df_points[shapely.contains_xy(region, df_points['lon'].values, df_points['lat'].values)]

    node_index = np.arange(len(df_points))
    lat = df_points['lat'].values
    lon = df_points['lon'].values
    d = np.zeros(len(lat), dtype=float)
    points_xyz = np.empty((len(lat), 3), dtype=float)
    rad = np.pi/180.0
    a = 6378.1370 # Earth's equatorial radius in km
    b = 6356.7523 # Earth's polar radius in km
    r = np.sqrt(((a**2 * np.cos(lat * rad))**2 + (b**2 * np.sin(lat * rad))**2) / ((a * np.cos(lat * rad))**2 + (b * np.sin(lat * rad))**2))
    points_xyz[:, 0] = (r - d) * np.cos(lat * rad) * np.cos(lon * rad)
    points_xyz[:, 1] = (r - d) * np.cos(lat * rad) * np.sin(lon * rad)
    points_xyz[:, 2] = (r - d) * np.sin(lat * rad)
    np.save(points_output_filename, points_xyz)
    np.save(node_index_output_filename, node_index)

    m_columns = [c for c in df_points.columns if c not in ['lon', 'lat', 'ss_wt', 'r_wt', 'n_wt']]
    m = np.asarray(m_columns, dtype=float)
    rates = df_points[m_columns].values.reshape(len(df_points) * len(m))
    node_index_all = np.repeat(node_index, len(m))
    m_array = np.tile(m, len(df_points))

    # fault_type: 1 = reverse, 2 = normal, 3 = strike slip
    rate_ss = np.repeat(df_points['ss_wt'].values, len(m)) * rates
    rate_rs = np.repeat(df_points['r_wt'].values, len(m)) * rates
    rate_ns = np.repeat(df_points['n_wt'].values, len(m)) * rates
    rate = np.hstack((rate_ss, rate_rs, rate_ns))
    style = np.hstack((np.full(len(rate_ss), 3), np.full(len(rate_rs), 1), np.full(len(rate_ns), 2)))
    node_index_all = np.hstack((node_index_all, node_index_all, node_index_all))
    m_all = np.hstack((m_array, m_array, m_array))
    keep = rate > 0
    np.savez_compressed(ruptures_output_filename, node_index=node_index_all[keep], m=m_all[keep], style=style[keep], rate=rate[keep])

    return

os.makedirs(output_dir, exist_ok=True)
get_points(
    os.path.join(input_dir, 'branch-avg-grid.csv'),
    os.path.join(input_dir, 'grid-system-active.geojson'),
    os.path.join(output_dir, 'points.npy'),
    os.path.join(output_dir, 'node_index.npy'),
    os.path.join(output_dir, 'ruptures.npz')
)

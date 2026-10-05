'''
Convert the 2023 USGS NSHM western U.S. fault system solution (branch averaged) into the
Numpy files used by ucla_plha. Output files have the same names and contents as the UCERF3
fault source models produced by convert_ucerf3_fault_source_models.py.

Input files (download into utilities/nshm23_wus_branch_avg/) from
https://code.usgs.gov/ghsc/nshmp/nshms/nshm-conus/-/tree/main/active-crust/fault/wus-system/branch-avg
    sections.geojson: fault subsection traces and properties
    ruptures.csv: one row per rupture with m, rate, depth (ztor, km), dip (deg), width (down-dip, km),
        rake (deg), and indices (subsection indices participating in the rupture)

Notation: NSHM23 uses "index" for fault subsections. As for UCERF3, we call these "segment_id",
and reserve "fault" for a named fault and "section" for a section of a fault.
'''
import json
import os
import numpy as np
import pandas as pd

input_dir = os.environ.get('NSHM23_WUS_INPUT_DIR', 'nshm23_wus_branch_avg')
output_dir = os.environ.get('NSHM23_WUS_OUTPUT_DIR', '../src/ucla_plha/source_models/fault_source_models/nshm23_wus')

# NSHM23 reduces the seismogenic area of each subsection by moving the upper depth down by
# aseismicity * (lower_depth - upper_depth). The rupture depths in ruptures.csv include this
# reduction, so apply it to the geometry too so that Rrup, Rjb, and ztor are consistent.
# As in nshmp-lib (DefaultGriddedSurface), the reduction removes the top of the dipping fault
# plane: the upper edge moves down-dip, by aseismicity * (lower_depth - upper_depth) / tan(dip)
# horizontally in the dip direction, and the lower edge does not move.
apply_aseismicity = True

def get_lat_lon(filename):
    '''
    Read NSHM23 sections geojson file and return Numpy arrays of segment_id and the latitude,
    longitude, and depth of the four corners of each planar piece of each subsection.
    Points 1 and 2 are the top of the piece along strike, and points 3 and 4 are directly
    down-dip of points 1 and 2, respectively.
    '''
    features = json.load(open(filename))['features']
    features = sorted(features, key=lambda f: f['properties']['index'])
    index = np.asarray([f['properties']['index'] for f in features])
    # ucla_plha indexes distance arrays by segment_id, so segment_id must be 0, 1, ..., N-1
    assert np.array_equal(index, np.arange(len(index)))

    segment_id = []
    lat1 = []
    lon1 = []
    lat2 = []
    lon2 = []
    dip = []
    dip_dir = []
    lower_depth = []
    upper_depth = []
    for f in features:
        p = f['properties']
        coords = np.asarray(f['geometry']['coordinates'])
        # A few traces repeat consecutive points. Remove them because zero-length pieces cause
        # division by zero when computing Rx, Rx1, and Ry0.
        coords = coords[np.r_[True, np.any(np.diff(coords, axis=0) != 0, axis=1)]]
        udepth = p['upper-depth']
        ldepth = p['lower-depth']
        shift = 0.0
        if apply_aseismicity:
            reduction = p.get('aseismicity', 0.0) * (ldepth - udepth)
            udepth = udepth + reduction
            shift = reduction / np.tan(np.radians(p['dip']))
        # top corners: the trace moved down-dip by the aseismic reduction
        tlat, tlon = destination_point(coords[:, 1], coords[:, 0], p['dip-direction'], shift)
        if shift == 0.0:
            tlat, tlon = coords[:, 1], coords[:, 0]
        for i in range(coords.shape[0] - 1):
            segment_id.append(p['index'])
            lon1.append(tlon[i])
            lat1.append(tlat[i])
            lon2.append(tlon[i + 1])
            lat2.append(tlat[i + 1])
            dip.append(p['dip'])
            dip_dir.append(p['dip-direction'])
            lower_depth.append(ldepth)
            upper_depth.append(udepth)

    segment_id = np.asarray(segment_id, dtype='intc')
    lat1 = np.asarray(lat1)
    lon1 = np.asarray(lon1)
    lat2 = np.asarray(lat2)
    lon2 = np.asarray(lon2)
    dip = np.asarray(dip)
    dip_dir = np.asarray(dip_dir)
    lower_depth = np.asarray(lower_depth)
    upper_depth = np.asarray(upper_depth)

    # Horizontal projection of the down-dip width, and the bottom corners found by moving
    # the top corners a distance w along the dip direction (spherical destination point
    # formula using the local radius of the oblate spheroid).
    w = (lower_depth - upper_depth) / np.tan(np.radians(dip))
    lat3, lon3 = destination_point(lat1, lon1, dip_dir, w)
    lat4, lon4 = destination_point(lat2, lon2, dip_dir, w)

    return(segment_id, lat1, lon1, upper_depth, lat2, lon2, upper_depth, lat3, lon3, lower_depth, lat4, lon4, lower_depth, dip)

def earth_radius(lat):
    '''
    Radius of oblate spheroid in km at latitude lat (degrees)
    '''
    rad = np.pi/180.0
    a = 6378.1370 # Earth's equatorial radius in km
    b = 6356.7523 # Earth's polar radius in km
    return np.sqrt(((a**2 * np.cos(lat * rad))**2 + (b**2 * np.sin(lat * rad))**2) / ((a * np.cos(lat * rad))**2 + (b * np.sin(lat * rad))**2))

def destination_point(lat, lon, azimuth, distance):
    '''
    Return latitude and longitude (degrees) of the point located distance (km) from (lat, lon)
    along azimuth (degrees clockwise from north).
    '''
    lat_r = np.radians(lat)
    lon_r = np.radians(lon)
    az = np.radians(azimuth)
    delta = distance / earth_radius(lat)
    lat_out = np.arcsin(np.sin(lat_r) * np.cos(delta) + np.cos(lat_r) * np.sin(delta) * np.cos(az))
    lon_out = lon_r + np.arctan2(np.sin(az) * np.sin(delta) * np.cos(lat_r), np.cos(delta) - np.sin(lat_r) * np.sin(lat_out))
    return(np.degrees(lat_out), np.degrees(lon_out))

def latlonel_to_xyz(geom):
    '''
    Accept N x M x 3 Numpy array of lat, lon, depth data, where N is the number of geometric objects,
    M is the number of points per geometric object, and 3 are the lat, lon, depth for the point.
    Convert to Cartesian coordinates assuming the earth is an oblate spheroid.
    Return N x M x 3 Numpy array of points in Cartesian coordinates
    '''
    xyz = np.empty(geom.shape)
    rad = np.pi/180.0
    d = geom[:, :, 2]
    lat = geom[:, :, 0]
    lon = geom[:, :, 1]
    r = earth_radius(lat)
    xyz[:, :, 0] = (r - d) * np.cos(lat * rad) * np.cos(lon * rad)
    xyz[:, :, 1] = (r - d) * np.cos(lat * rad) * np.sin(lon * rad)
    xyz[:, :, 2] = (r - d) * np.sin(lat * rad)
    return xyz

def get_triangles(lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4):
    '''
    Accept Numpy arrays of latitude, longitude, and depth for four points defining the segment corners
    Return a float array of triangles for computing Rrup, and a float array of triangles for computing Rjb
    '''
    tri1_rrup = np.asarray([[lat1, lat2, lat4], [lon1, lon2, lon4], [d1, d2, d4]]).T
    tri2_rrup = np.asarray([[lat1, lat3, lat4], [lon1, lon3, lon4], [d1, d3, d4]]).T
    tri_rrup = np.concatenate((tri1_rrup, tri2_rrup))
    tri_rrup_xyz = latlonel_to_xyz(tri_rrup)
    zero_depth = np.zeros(len(lat1))
    tri1_rjb = np.asarray([[lat1, lat2, lat4], [lon1, lon2, lon4], [zero_depth, zero_depth, zero_depth]]).T
    tri2_rjb = np.asarray([[lat1, lat3, lat4], [lon1, lon3, lon4], [zero_depth, zero_depth, zero_depth]]).T
    tri_rjb = np.concatenate((tri1_rjb, tri2_rjb))
    tri_rjb_xyz = latlonel_to_xyz(tri_rjb)
    return(tri_rrup_xyz, tri_rjb_xyz)

def get_rectangles(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4):
    '''
    Return a float array of rectangles for surface projection of fault for computing Rx, Rx1, and Ry0
    '''
    zero_depth = np.zeros(len(lat1))
    rect_rjb = np.asarray([[lat1, lat2, lat3, lat4], [lon1, lon2, lon3, lon4], [zero_depth, zero_depth, zero_depth, zero_depth]]).T
    rect_rjb_xyz = latlonel_to_xyz(rect_rjb)
    return(rect_rjb_xyz)

def parse_indices(indices):
    '''
    Parse an NSHM23 rupture indices string into a list of subsection indices. Groups are separated
    by "-", and each group is either a single index or an inclusive range "a:b" that may be
    ascending or descending (e.g., "2168:2169-1068:1056-17").
    '''
    out = []
    for group in indices.split('-'):
        if ':' in group:
            a, b = map(int, group.split(':'))
            step = 1 if b >= a else -1
            out.extend(range(a, b + step, step))
        else:
            out.append(int(group))
    return out

def get_ruptures_segments(rupture_df, output_file):
    '''
    Organize the subsections participating in each rupture into Numpy arrays of rupture_index
    and segment_index, and save in compressed npz format.
    '''
    segment_index_list = [parse_indices(s) for s in rupture_df['indices'].values]
    num_sections = np.asarray([len(s) for s in segment_index_list])
    rupture_index_all = np.repeat(np.arange(len(rupture_df), dtype=np.int32), num_sections)
    segment_index_all = np.fromiter((i for s in segment_index_list for i in s), dtype=np.int32, count=num_sections.sum())
    np.savez_compressed(output_file, rupture_index=rupture_index_all, segment_index=segment_index_all)
    return

def get_rupture_data(rupture_df, output_file):
    '''
    Save magnitude, rate, style of faulting, dip, ztor, and zbor for each rupture. NSHM23 provides
    rupture-level depth to top (km), dip (degrees), and down-dip width (km), so zbor = ztor + width * sin(dip).
    '''
    m = rupture_df['m'].values
    rate = rupture_df['rate'].values
    rake = rupture_df['rake'].values
    dip = rupture_df['dip'].values
    ztor = rupture_df['depth'].values
    zbor = ztor + rupture_df['width'].values * np.sin(np.radians(dip))
    # fault_type: 1 = reverse, 2 = normal, 3 = strike slip
    fault_type = np.full(len(rupture_df), 1)
    fault_type[(rake > -150) & (rake < -30)] = 2
    fault_type[(rake >= -180) & (rake <= -150)] = 3
    fault_type[(rake >= -30) & (rake <= 30)] = 3
    fault_type[(rake >= 150) & (rake <= 180)] = 3
    np.savez_compressed(output_file, m=m, rate=rate, fault_type=fault_type, dip=dip, ztor=ztor, zbor=zbor)
    return

os.makedirs(output_dir, exist_ok=True)

### Compute triangles representing fault segments, and array of segment_id values
segment_id, lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4, dip = get_lat_lon(os.path.join(input_dir, 'sections.geojson'))
tri = get_triangles(lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4)
np.save(os.path.join(output_dir, 'tri_segment_id.npy'), np.concatenate((segment_id, segment_id)))
np.save(os.path.join(output_dir, 'tri_rrup.npy'), tri[0])
np.save(os.path.join(output_dir, 'tri_rjb.npy'), tri[1])
rect = get_rectangles(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4)
np.save(os.path.join(output_dir, 'rect_segment_id.npy'), segment_id)
np.save(os.path.join(output_dir, 'rect_rjb.npy'), rect)

### Compute array mapping rupture and segment indices, and rupture properties
rupture_df = pd.read_csv(os.path.join(input_dir, 'ruptures.csv'))
get_ruptures_segments(rupture_df, os.path.join(output_dir, 'ruptures_segments.npz'))
get_rupture_data(rupture_df, os.path.join(output_dir, 'ruptures.npz'))

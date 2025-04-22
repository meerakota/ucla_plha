import numpy as np
import pandas as pd
import pandas_geojson as pdg

def get_lat_lon(filename):
    '''
    Read UCERF3 source geojson file, convert it to a Pandas dataframe, and gather together data
    and output a Python list containing Numpy arrays for segment_id, and the latitudes and longitudes.

    Notation: UCERF3 uses the variable name "FaultID" to refer to the subsections of the fault.
    We reserve the word "fault" for a particular named fault (e.g., Airport Lake), and "section" for
    a section of the fault. We therefore have renamed "FaultID" to "segment_id" in our dataframe.
    '''
    df = pdg.read_geojson(filename).to_dataframe()
    geometry_coordinates = df['geometry.coordinates'].values
    FaultID = df['properties.FaultID'].values
    DipDeg = df['properties.DipDeg'].values
    DipDir = df['properties.DipDir'].values
    LowDepth = df['properties.LowDepth'].values
    UpDepth = df['properties.UpDepth'].values
    
    lat1 = []
    lon1 = []
    lat2 = []
    lon2 = []
    segment_id = []
    dip = []
    dip_dir = []
    lower_depth = []
    upper_depth = []
    
    for gm, fid, ddeg, ddir, ldepth, udepth in zip(geometry_coordinates,FaultID,DipDeg,DipDir,LowDepth,UpDepth):
        gm_array = np.asarray(gm)
        for i in range(gm_array.shape[0]-1):
            segment_id.append(fid)
            lon1.append(gm_array[i,0])
            lat1.append(gm_array[i,1])
            lon2.append(gm_array[i+1,0])
            lat2.append(gm_array[i+1,1])
            dip.append(ddeg)
            dip_dir.append(ddir)
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

    rad = np.pi/180
    a = 6378.1370 # Earth's equatorial radius in km
    b = 6356.7523 # Earth's polar radius in km
    w = (lower_depth - upper_depth) / np.tan(dip * rad)
    wy = w * np.cos(dip_dir * rad)
    r1 = np.sqrt(((a**2 * np.cos(lat1 * rad))**2 + (b**2 * np.sin(lat1 * rad))**2) / ((a * np.cos(lat1*rad))**2 + (b * np.sin(lat1*rad))**2))
    r2 = np.sqrt(((a**2 * np.cos(lat2 * rad))**2 + (b**2 * np.sin(lat2 * rad))**2) / ((a * np.cos(lat2*rad))**2 + (b * np.sin(lat2*rad))**2))
    
    dlat1 = wy / r1 / rad
    dlat2 = wy / r2 / rad
    lat3 = lat1 + dlat1
    lat4 = lat2 + dlat2
    r3 = np.sqrt(((a**2 * np.cos(lat3 * rad))**2 + (b**2 * np.sin(lat3 * rad))**2) / ((a * np.cos(lat3*rad))**2 + (b * np.sin(lat3*rad))**2))
    r4 = np.sqrt(((a**2 * np.cos(lat4 * rad))**2 + (b**2 * np.sin(lat4 * rad))**2) / ((a * np.cos(lat4*rad))**2 + (b * np.sin(lat4*rad))**2))
    
    invangle13 = 1.0 - (2.0 * np.sin(w/(2.0*r3))**2.0 + np.cos(dlat1 * rad) - 1.0) / (np.cos(lat1 * rad) * np.cos(lat3 * rad))
    # floating point precision may render inverse angles higher than 1.0 (or less than -1.0, which doesn't happen here but 
    # is mathematically possible). So impose ranges on the inverse angles
    invangle13[invangle13 > 1.0] = 1.0
    invangle13[invangle13 < -1.0] = -1.0
    lon3 = lon1 + np.arccos(invangle13) / rad
    
    invangle24 = 1.0 - (2.0 * np.sin(w/(2.0*r4))**2.0 + np.cos(dlat2 * rad) - 1.0) / (np.cos(lat2 * rad) * np.cos(lat4 * rad))
    invangle24[invangle24 > 1.0] = 1.0
    invangle24[invangle24 < -1.0] = -1.0
    lon4 = lon2 + np.arccos(invangle24) / rad

    return(segment_id, lat1, lon1, upper_depth, lat2, lon2, upper_depth, lat3, lon3, lower_depth, lat4, lon4, lower_depth, dip)

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
    a = 6378.1370 # Earth's equatorial radius in km
    b = 6356.7523 # Earth's polar radius in km
    lat = geom[:, :, 0]
    lon = geom[:, :, 1]
    r = np.sqrt(((a**2 * np.cos(lat * rad))**2 + (b**2 * np.sin(lat * rad))**2) / ((a * np.cos(lat * rad))**2 + (b * np.sin(lat * rad))**2))
    xyz[:, :, 0] = (
        (r - d)
        * np.cos(lat * rad)
        * np.cos(lon * rad)
    )
    xyz[:, :, 1] = (
        (r - d)
        * np.cos(lat * rad)
        * np.sin(lon * rad)
    )
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
    Accept lat-lon array from the get_lat_lon() function, and return a Python list containing
    an integer array of segment_id values, a float array of rectangles for surface projection 
    of fault for computing Rx, Rx1, and Ry0
    '''
    zero_depth = np.zeros(len(lat1))
    rect_rjb = np.asarray([[lat1, lat2, lat3, lat4], [lon1, lon2, lon3, lon4], [zero_depth, zero_depth, zero_depth, zero_depth]]).T
    rect_rjb_xyz = latlonel_to_xyz(rect_rjb)
    return(rect_rjb_xyz)

def get_rupture_data(rupture_file, rate_file, ruptures_segments_file, ztor_segments, zbor_segments, dip_segments, output_file):
    '''
    Read UCERF3 rupture data file, and rate data file. Organize data into a Pandas dataframe containing
    rupture_index, magnitude, rate, and style of faulting. Save Pandas dataframe in pickle file format. 
    '''
    rupture_df =  pd.read_csv(rupture_file)
    rate_df = pd.read_csv(rate_file)
    ruptures_segments_df = pd.read_pickle(ruptures_segments_file, compression='gzip')
    segment_index = ruptures_segments_df['segment_index'].values
    ruptures_segments_df['ztor_all'] = ztor_segments[segment_index]
    ruptures_segments_df['zbor_all'] = zbor_segments[segment_index]
    ruptures_segments_df['dip_all'] = dip_segments[segment_index]
    ztor = ruptures_segments_df.groupby('rupture_index')['ztor_all'].min().values
    zbor = ruptures_segments_df.groupby('rupture_index')['zbor_all'].max().values
    dip = ruptures_segments_df.groupby('rupture_index')['dip_all'].mean().values
    rate = rate_df['Annual Rate'].values
    fault_type = np.full(len(rupture_df), 1)
    rake = rupture_df['Average Rake (degrees)'].values
    m = rupture_df['Magnitude'].values
    fault_type[(rake > -150) & (rake < -30)] = 2
    fault_type[(rake >= -180) & (rake <= -150)] = 3
    fault_type[(rake >= -30) & (rake <= 30)] = 3
    fault_type[(rake >= 150) & (rake <= 180)] = 3
    df_out = pd.DataFrame()
    df_out['rupture_index'] = rupture_df['Rupture Index'].values
    df_out['m'] = m
    df_out['rate'] = rate
    df_out['style'] = fault_type
    df_out['rake'] = rake
    df_out['dip'] = dip
    df_out['ztor'] = ztor
    df_out['zbor'] = zbor
    df_out.to_pickle(output_file, compression='gzip')
    return

def get_ruptures_segments(rupture_indices_file, output_file):
    '''
    Read UCERF3 file containing the list of all of the segments associated with each rupture.
    Organize the data into a single Numpy array containing rupture_index and segment_index.
    The array is very large, so use the smallest possible integer container, and gzip the pickle files.
    '''
    indices = pd.read_csv(rupture_indices_file, engine='python')
    rupture_index = indices['Rupture Index'].values
    segment_index = indices['Num Sections'].values
    indicesT = indices.drop(labels=['Rupture Index', 'Num Sections'], axis=1).transpose()
    segment_index_array = np.empty(len(rupture_index), dtype='object')
    rupture_index_array = np.empty(len(rupture_index), dtype='object')
    for i, it in indicesT.items():
        segment_index_array[i] = np.asarray(it.values[0:segment_index[i]], dtype='short')
        rupture_index_array[i] = np.full(segment_index[i], rupture_index[i])
    rupture_index_all = np.hstack(rupture_index_array)
    segment_index_all = np.hstack(segment_index_array)
    ruptures_segments_df = pd.DataFrame()
    ruptures_segments_df['rupture_index'] = rupture_index_all.astype('intc')
    ruptures_segments_df['segment_index'] = segment_index_all.astype('short')
    ruptures_segments_df.to_pickle(output_file, compression='gzip')


### Compute triangles representing fault segments, and array of segment_id values
segment_id, lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4, dip = get_lat_lon('FM3_1_branch_averaged/ruptures/fault_sections.geojson')
tri_fm31 = get_triangles(lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4)
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm31/tri_segment_id.npy', np.concatenate((segment_id, segment_id)))
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm31/tri_rrup.npy', tri_fm31[0])
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm31/tri_rjb.npy', tri_fm31[1])
rect_fm31 = get_rectangles(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4)
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm31/rect_segment_id.npy', segment_id)
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm31/rect_rjb.npy', rect_fm31)

### Compute array mapping rupture and segment indices
fm31_rupture_indices_file = 'FM3_1_branch_averaged/ruptures/indices.csv'
fm31_output_file = '../src/ucla_plha/source_models/fault_source_models/ucerf3_fm31/ruptures_segments.pkl'
get_ruptures_segments(fm31_rupture_indices_file, fm31_output_file)

### Compute ruptures.pkl file that contains rupture_index, magnitude, rate and style of faulting for each event
fm31_rupture_file = 'FM3_1_branch_averaged/ruptures/properties.csv'
fm31_rate_file = 'FM3_1_branch_averaged/solution/rates.csv'
fm31_ruptures_segments_file = '../src/ucla_plha/source_models/fault_source_models/ucerf3_fm31/ruptures_segments.pkl'
fm31_output_file = '../src/ucla_plha/source_models/fault_source_models/ucerf3_fm31/ruptures.pkl'
get_rupture_data(fm31_rupture_file, fm31_rate_file, fm31_ruptures_segments_file, d1, d3, dip, fm31_output_file)

# Now repeat for fm32
segment_id, lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4, dip = get_lat_lon('FM3_2_branch_averaged/ruptures/fault_sections.geojson')
tri_fm32 = get_triangles(lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4)
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm32/tri_segment_id.npy', np.concatenate((segment_id, segment_id)))
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm32/tri_rrup.npy', tri_fm32[0])
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm32/tri_rjb.npy', tri_fm32[1])
rect_fm32 = get_rectangles(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4)
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm32/rect_segment_id.npy', segment_id)
np.save('../src/ucla_plha/source_models/fault_source_models/ucerf3_fm32/rect_rjb.npy', rect_fm32)

fm32_rupture_indices_file = 'FM3_2_branch_averaged/ruptures/indices.csv'
fm32_output_file = '../src/ucla_plha/source_models/fault_source_models/ucerf3_fm32/ruptures_segments.pkl'
get_ruptures_segments(fm32_rupture_indices_file, fm32_output_file)

fm32_rupture_file = 'FM3_2_branch_averaged/ruptures/properties.csv'
fm32_rate_file = 'FM3_2_branch_averaged/solution/rates.csv'
fm32_ruptures_segments_file = '../src/ucla_plha/source_models/fault_source_models/ucerf3_fm32/ruptures_segments.pkl'
fm32_output_file = '../src/ucla_plha/source_models/fault_source_models/ucerf3_fm32/ruptures.pkl'
get_rupture_data(fm32_rupture_file, fm32_rate_file, fm32_ruptures_segments_file, d1, d3, dip, fm32_output_file)
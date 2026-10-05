'''
Convert the Cascadia subduction interface source model of the USGS NSHM for the conterminous U.S.
(nshm-conus release 6.2.0, subduction/interface) into the Numpy files used by ucla_plha. The
interface logic tree is collapsed to mean rates, so each output rupture rate is the NSHM rupture
rate multiplied by the weights of all logic tree branches leading to it. Geometry branches (top,
middle, and bottom depth models) and rupture branches (full, segmented, and unsegmented models)
are kept as separate ruptures.

Two source models are written:
    fault_source_models/nshm23_cascadia_interface: independent (Poissonian) interface ruptures.
        These are the ruptures that the USGS hazard web service reports as the "Interface" component.
    fault_source_models/nshm23_cascadia_interface_cluster: the full-rupture "cluster-in" branch
        (weight 0.1), in which Cascadia ruptures as a cluster (sequence) of 7 or 8 adjacent M~8
        ruptures. nshmp-lib builds these as ClusterRuptureSets (source type FAULT_CLUSTER), which the
        hazard web service reports in the "FaultCluster" component. They are not independent ruptures,
        see source_info.json in that directory for how nshmp-lib computes cluster hazard.

Input files are downloaded from https://code.usgs.gov/ghsc/nshmp/nshms/nshm-conus (tag 6.2.0) into
utilities/nshm23_cascadia_interface/ if they are not already there.

The conversion follows nshmp-lib (https://code.usgs.gov/ghsc/nshmp/nshmp-lib, commit 44728a7d):
    ModelLoader.Interface, Deserialize, SourceTree: logic tree traversal. source-tree.json branches
        have weights and source-group.json branches have scale factors (additive branches); the
        weight of a rupture set is the product of the branch weights and scales from the root to
        the leaf, rounded to 8 decimal places.
    InterfaceRuptureSet: the sections of a rupture set are joined into an upper and a lower trace
        (removing duplicate points), and an ApproxGriddedSurface with 5 km spacing (interface-config.json
        "surface-spacing") is built between them. MFD branch rates are scaled by the MFD branch weight.
        SINGLE MFDs rupture the full surface (RuptureFloating.OFF). GR MFDs (unsegmented branches) float
        along strike only (RuptureFloating.STRIKE_ONLY) over the full down-dip width, with rupture length
        from RuptureScaling.NSHM_SUB_GEOMAT_LENGTH, L = 10^((M - 4.94) / 1.39) km, and the rate of each
        magnitude divided equally among the floating ruptures. GR magnitudes are bin centers from
        mMin + dm / 2 to mMax - dm / 2 with incremental rate 10^(a - b M).
    ApproxGriddedSurface, GriddedSubsetSurface: zTor is the mean depth of the upper trace vertices
        (full ruptures) or of the top row of grid points (floating ruptures), dip is the average of the
        plunges from the first and last upper trace points to the first and last lower trace points,
        and width = (rows - 1) * spacing. The GMM inputs are dip, width, zTor, zHyp = zTor + width / 2 *
        sin(dip), and rake = 90.

ucla_plha represents each surface by "segments": the strips between adjacent columns of the
nshmp-lib grid of points, each with four corners (two on the top row and two on the bottom row of
the grid). A full rupture uses all strips of its surface and a floating rupture uses the strips
between its first and last grid columns. nshmp-lib computes rRup and rJB as the minimum distance
to the grid points (5 km spacing); ucla_plha computes the distance to the two triangles of each
strip, which is slightly smaller (by about 0.1 km at 50 km distance).

Identical ruptures (same surface, columns, and magnitude) on different branches are merged by
adding their rates.
'''
import json
import os
import urllib.parse
import urllib.request
import numpy as np

nshm_ref = '6.2.0'
nshm_project = 6063
nshm_path = 'subduction/interface'
input_dir = 'nshm23_cascadia_interface'
output_dir = '../src/ucla_plha/source_models/fault_source_models/nshm23_cascadia_interface'
output_dir_cluster = '../src/ucla_plha/source_models/fault_source_models/nshm23_cascadia_interface_cluster'

EARTH_RADIUS_MEAN = 6371.0072  # km, nshmp-lib Coordinates.EARTH_RADIUS_MEAN
RAKE = 90.0  # InterfaceRuptureSet passes rake = 90 to the rupture floating model


def gitlab_get(url):
    # code.usgs.gov rejects the default Python user agent
    request = urllib.request.Request(url, headers={'User-Agent': 'curl/8.0'})
    return urllib.request.urlopen(request).read()


def download_inputs():
    '''
    Download the subduction interface source model files from code.usgs.gov unless they exist.
    '''
    if os.path.exists(os.path.join(input_dir, 'interface-config.json')):
        return
    api = f'https://code.usgs.gov/api/v4/projects/{nshm_project}/repository'
    blobs = []
    page = 1
    while True:
        tree = json.loads(gitlab_get(f'{api}/tree?path={nshm_path}&ref={nshm_ref}&per_page=100&recursive=true&page={page}'))
        if not tree:
            break
        blobs += [t['path'] for t in tree if t['type'] == 'blob']
        page += 1
    for path in blobs:
        out = os.path.join(input_dir, os.path.relpath(path, nshm_path))
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, 'wb') as f:
            f.write(gitlab_get(f'{api}/files/{urllib.parse.quote(path, safe="")}/raw?ref={nshm_ref}'))


### Geodesy, ported from nshmp-lib geo.Locations and geo.LocationVector (spherical earth).
### Locations are (lon, lat, depth) tuples in degrees and km.

def horz_distance_fast(p1, p2):
    lat1, lat2 = np.radians(p1[1]), np.radians(p2[1])
    dlat = lat1 - lat2
    dlon = (np.radians(p1[0]) - np.radians(p2[0])) * np.cos((lat1 + lat2) * 0.5)
    return EARTH_RADIUS_MEAN * np.sqrt(dlat * dlat + dlon * dlon)


def horz_distance(p1, p2):
    lat1, lat2 = np.radians(p1[1]), np.radians(p2[1])
    s_dlat = np.sin((lat2 - lat1) / 2.0)
    s_dlon = np.sin((np.radians(p2[0]) - np.radians(p1[0])) / 2.0)
    c = s_dlat * s_dlat + np.cos(lat1) * np.cos(lat2) * s_dlon * s_dlon
    return EARTH_RADIUS_MEAN * 2.0 * np.arctan2(np.sqrt(c), np.sqrt(1.0 - c))


def linear_distance_fast(p1, p2):
    h = horz_distance_fast(p1, p2)
    v = p2[2] - p1[2]
    return np.sqrt(h * h + v * v)


def azimuth_rad(p1, p2):
    lat1, lat2 = np.radians(p1[1]), np.radians(p2[1])
    dlon = np.radians(p2[0]) - np.radians(p1[0])
    az = np.arctan2(np.sin(dlon) * np.cos(lat2), np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon))
    return (az + 2.0 * np.pi) % (2.0 * np.pi)


def vector(p1, p2):
    # LocationVector.create(p1, p2): azimuth (rad), horizontal distance, vertical distance
    return (azimuth_rad(p1, p2), horz_distance(p1, p2), p2[2] - p1[2])


def location(p, az, dh, dv):
    lat1 = np.radians(p[1])
    ad = dh / EARTH_RADIUS_MEAN
    lat2 = np.arcsin(np.sin(lat1) * np.cos(ad) + np.cos(lat1) * np.sin(ad) * np.cos(az))
    lon2 = np.radians(p[0]) + np.arctan2(np.sin(az) * np.sin(ad) * np.cos(lat1), np.cos(ad) - np.sin(lat1) * np.sin(lat2))
    return (np.degrees(lon2), np.degrees(lat2), p[2] + dv)


def trace_length(trace):
    return sum(horz_distance_fast(trace[i], trace[i + 1]) for i in range(len(trace) - 1))


def java_round(x):
    # Math.round
    return int(np.floor(x + 0.5))


def resample_trace(trace, num):
    '''
    Port of nshmp-lib Surfaces.resampleTrace
    '''
    resamp_int = trace_length(trace) / num
    locs = [trace[0]]
    remaining = resamp_int
    last = trace[0]
    i = 1
    while i < len(trace):
        nxt = trace[i]
        length = linear_distance_fast(last, nxt)
        if length > remaining:
            az, dh, dv = vector(last, nxt)
            loc = location(last, az, dh * remaining / length, dv * remaining / length)
            locs.append(loc)
            last = loc
            remaining = resamp_int
        else:
            last = nxt
            i += 1
            remaining -= length
    if linear_distance_fast(trace[-1], locs[-1]) > resamp_int / 2:
        locs.append(trace[-1])
    return locs


class Surface:
    '''
    Port of nshmp-lib ApproxGriddedSurface. grid[r, c] = (lon, lat, depth).
    '''

    def __init__(self, upper, lower, spacing):
        self.upper = upper
        self.lower = lower
        self.spacing = spacing
        num = java_round((trace_length(upper) + trace_length(lower)) / 2 / spacing)
        top = resample_trace(upper, num)
        bot = resample_trace(lower, num)
        assert len(top) == len(bot) == num + 1, (len(top), len(bot), num)
        ave_dist = np.mean([linear_distance_fast(t, b) for t, b in zip(top, bot)])
        n_rows = java_round(ave_dist / spacing) + 1
        grid = np.empty((n_rows, num + 1, 3))
        for c in range(num + 1):
            az, dh, dv = vector(top[c], bot[c])
            dh_incr = dh / (n_rows - 1)
            dv_incr = dv / (n_rows - 1)
            loc = top[c]
            grid[0, c] = loc
            for r in range(1, n_rows):
                loc = location(loc, az, dh_incr, dv_incr)
                grid[r, c] = loc
        self.grid = grid
        self.n_rows = n_rows
        self.n_cols = num + 1
        v_first = vector(upper[0], lower[0])
        v_last = vector(upper[-1], lower[-1])
        self.dip = np.degrees((np.arctan(v_first[2] / v_first[1]) + np.arctan(v_last[2] / v_last[1])) / 2)
        self.depth = np.mean([p[2] for p in upper])  # LocationList.depth() of the upper trace
        self.width = spacing * (n_rows - 1)


def float_strike_only(surface, m):
    '''
    Port of nshmp-lib RuptureFloating.STRIKE_ONLY with RuptureScaling.NSHM_SUB_GEOMAT_LENGTH.
    Returns a list of (start column, number of columns) of the floating ruptures.
    '''
    length = 10.0 ** ((m - 4.94) / 1.39)
    col_size = int(np.rint(length / surface.spacing + 1))  # Math.rint
    along_count = surface.n_cols - col_size + 1
    if along_count <= 1:
        along_count = 1
        col_size = surface.n_cols
    # Down-dip: floaters extend over the full width, so there is one row position
    return [(c, col_size) for c in range(along_count)]


def read_json(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def read_features(dirname):
    features = {}
    for name in os.listdir(dirname):
        f = read_json(os.path.join(dirname, name))
        traces = f['geometry']['coordinates']
        features[f['id']] = [[tuple(p) for p in t] for t in traces]
    return features


def join_traces(features, sections):
    '''
    InterfaceRuptureSet.joinInterfaceSections: concatenate the upper (index 0) and lower (index -1)
    traces of the sections in order and remove duplicate points.
    '''
    out = []
    for index in [0, -1]:
        trace = []
        for s in sections:
            for p in features[s][index]:
                if p not in trace:
                    trace.append(p)
        out.append(trace)
    return out


def read_tree(dirname):
    '''
    Return a list of (child directory name, weight) for a source-tree.json or source-group.json,
    or None for a leaf directory.
    '''
    if os.path.exists(os.path.join(dirname, 'source-tree.json')):
        return [(b['id'], b['weight']) for b in read_json(os.path.join(dirname, 'source-tree.json'))]
    if os.path.exists(os.path.join(dirname, 'source-group.json')):
        return [(b['id'], b.get('scale', 1.0)) for b in read_json(os.path.join(dirname, 'source-group.json'))]
    return None


def walk(dirname, weights, path, leaves):
    '''
    Collect logic tree leaves (directory, list of branch weights, path of branch ids).
    '''
    children = read_tree(dirname)
    if children is None:
        leaves.append((dirname, weights, path))
        return
    for name, w in children:
        walk(os.path.join(dirname, name), weights + [w], path + [name], leaves)


def mfd_branches(mfd_tree):
    '''
    Return a list of (magnitudes, rates) scaled by the MFD branch weights, following
    InterfaceRuptureSet.buildMfdTree with mfd-config.json epistemic-tree and aleatory-properties null.
    '''
    out = []
    for branch in mfd_tree:
        v = branch['value']
        if v['type'] == 'SINGLE':
            # A SINGLE MFD without a rate (cluster rupture sets) has rate 1, so the rate is the branch weight
            out.append((np.array([v['m']]), np.array([v.get('rate', 1.0)]) * branch['weight'], 'SINGLE'))
        elif v['type'] == 'GR':
            dm = v['Δm']
            n = int(round((v['mMax'] - v['mMin']) / dm))
            m = np.round(v['mMin'] + dm / 2.0 + dm * np.arange(n), 5)
            out.append((m, 10.0 ** (v['a'] - v['b'] * m) * branch['weight'], 'GR'))
        else:
            raise ValueError(v['type'])
    return out


class Model:
    '''
    Accumulate surfaces (as strips between grid columns) and ruptures.
    '''

    def __init__(self, features, spacing):
        self.features = features
        self.spacing = spacing
        self.surfaces = {}  # tuple(sections) -> (Surface, first segment id)
        self.n_segments = 0
        self.ruptures = {}  # key -> dict of rupture properties

    def surface(self, sections):
        key = tuple(sections)
        if key not in self.surfaces:
            upper, lower = join_traces(self.features, sections)
            s = Surface(upper, lower, self.spacing)
            self.surfaces[key] = (s, self.n_segments)
            self.n_segments += s.n_cols - 1
        return self.surfaces[key]

    def add(self, sections, m, rate, start_col=None, col_size=None, extra=()):
        s, seg0 = self.surface(sections)
        if start_col is None:
            # full rupture: zTor is the mean depth of the upper trace vertices
            start_col, col_size, ztor = 0, s.n_cols, s.depth
        else:
            ztor = np.mean(s.grid[0, start_col:start_col + col_size, 2])
        key = (tuple(sections), start_col, col_size, round(float(m), 5)) + tuple(extra)
        if key in self.ruptures:
            self.ruptures[key]['rate'] += rate
            return
        self.ruptures[key] = dict(
            m=float(m), rate=rate, dip=s.dip, ztor=ztor, width=s.width,
            segments=np.arange(seg0 + start_col, seg0 + start_col + col_size - 1),
            extra=extra)

    def geometry(self):
        '''
        Return the four corners (lat, lon, depth) of each strip, ordered by segment id. Points 1 and 2
        are on the top row of the grid, and points 3 and 4 are on the bottom row below points 1 and 2.
        '''
        corners = np.empty((self.n_segments, 4, 3))
        for s, seg0 in self.surfaces.values():
            g = s.grid[:, :, [1, 0, 2]]  # lat, lon, depth
            n = s.n_cols - 1
            corners[seg0:seg0 + n, 0] = g[0, :-1]
            corners[seg0:seg0 + n, 1] = g[0, 1:]
            corners[seg0:seg0 + n, 2] = g[-1, :-1]
            corners[seg0:seg0 + n, 3] = g[-1, 1:]
        return corners


def earth_radius(lat):
    '''
    Radius of oblate spheroid in km at latitude lat (degrees)
    '''
    rad = np.pi/180.0
    a = 6378.1370 # Earth's equatorial radius in km
    b = 6356.7523 # Earth's polar radius in km
    return np.sqrt(((a**2 * np.cos(lat * rad))**2 + (b**2 * np.sin(lat * rad))**2) / ((a * np.cos(lat * rad))**2 + (b * np.sin(lat * rad))**2))


def latlonel_to_xyz(geom):
    '''
    Convert N x M x 3 array of lat, lon, depth to Cartesian coordinates (same as the other converters)
    '''
    xyz = np.empty(geom.shape)
    rad = np.pi/180.0
    lat = geom[:, :, 0]
    lon = geom[:, :, 1]
    r = earth_radius(lat) - geom[:, :, 2]
    xyz[:, :, 0] = r * np.cos(lat * rad) * np.cos(lon * rad)
    xyz[:, :, 1] = r * np.cos(lat * rad) * np.sin(lon * rad)
    xyz[:, :, 2] = r * np.sin(lat * rad)
    return xyz


def write_model(model, out_dir, extra_names=()):
    '''
    Write the geometry and rupture files in the same format as the other fault source models.
    '''
    os.makedirs(out_dir, exist_ok=True)
    corners = model.geometry()
    segment_id = np.arange(model.n_segments, dtype='intc')
    c1, c2, c3, c4 = (corners[:, i] for i in range(4))
    tri_rrup = np.concatenate((np.stack((c1, c2, c4), axis=1), np.stack((c1, c3, c4), axis=1)))
    tri_rjb = tri_rrup.copy()
    tri_rjb[:, :, 2] = 0.0
    rect = corners.copy()
    rect[:, :, 2] = 0.0
    np.save(os.path.join(out_dir, 'tri_segment_id.npy'), np.concatenate((segment_id, segment_id)))
    np.save(os.path.join(out_dir, 'tri_rrup.npy'), latlonel_to_xyz(tri_rrup))
    np.save(os.path.join(out_dir, 'tri_rjb.npy'), latlonel_to_xyz(tri_rjb))
    np.save(os.path.join(out_dir, 'rect_segment_id.npy'), segment_id)
    np.save(os.path.join(out_dir, 'rect_rjb.npy'), latlonel_to_xyz(rect))

    ruptures = list(model.ruptures.values())
    segments = [r['segments'] for r in ruptures]
    np.savez_compressed(
        os.path.join(out_dir, 'ruptures_segments.npz'),
        rupture_index=np.repeat(np.arange(len(ruptures), dtype=np.int32), [len(s) for s in segments]),
        segment_index=np.concatenate(segments).astype(np.int32))
    m = np.array([r['m'] for r in ruptures])
    dip = np.array([r['dip'] for r in ruptures])
    ztor = np.array([r['ztor'] for r in ruptures])
    width = np.array([r['width'] for r in ruptures])
    data = dict(
        m=m,
        rate=np.array([r['rate'] for r in ruptures]),
        # fault_type: 1 = reverse (interface events, rake 90)
        fault_type=np.full(len(ruptures), 1),
        dip=dip,
        ztor=ztor,
        # zbor from the nshmp-lib (nominal) down-dip width, so width = (zbor - ztor) / sin(dip)
        zbor=ztor + width * np.sin(np.radians(dip)),
        width=width,
        zhyp=ztor + width / 2.0 * np.sin(np.radians(dip)),
        rake=np.full(len(ruptures), RAKE))
    for i, name in enumerate(extra_names):
        data[name] = np.array([r['extra'][i] for r in ruptures], dtype=np.int32)
    np.savez_compressed(os.path.join(out_dir, 'ruptures.npz'), **data)
    return data


def main():
    download_inputs()
    config = read_json(os.path.join(input_dir, 'interface-config.json'))
    mfd_config = read_json(os.path.join(input_dir, 'mfd-config.json'))
    assert config['rupture-floating'] == 'STRIKE_ONLY'
    assert config['rupture-scaling'] == 'NSHM_SUB_GEOMAT_LENGTH'
    assert mfd_config['epistemic-tree'] is None and mfd_config['aleatory-properties'] is None
    root = os.path.join(input_dir, 'Cascadia')
    mfd_map = read_json(os.path.join(root, 'mfd-map.json'))
    features = read_features(os.path.join(root, 'features'))

    leaves = []
    walk(root, [], [], leaves)

    interface = Model(features, config['surface-spacing'])
    cluster = Model(features, config['surface-spacing'])
    clusters = []  # (cluster_id, rate, weight, name)
    n_cluster_sections = 0
    branch_rates = {}
    for dirname, weights, path in leaves:
        # Leaf weights are rounded to 8 decimal places in nshmp-lib SourceTree
        weight = round(float(np.prod(weights)), 8)
        if os.path.exists(os.path.join(dirname, 'rupture-set.json')):
            rs = read_json(os.path.join(dirname, 'rupture-set.json'))
            total = 0.0
            for m, rates, mfd_type in mfd_branches(mfd_map[rs['mfd-tree']]):
                for mi, ri in zip(m, rates):
                    if mfd_type == 'SINGLE':
                        interface.add(rs['sections'], mi, weight * ri)
                    else:
                        s, _ = interface.surface(rs['sections'])
                        floaters = float_strike_only(s, mi)
                        for start_col, col_size in floaters:
                            interface.add(rs['sections'], mi, weight * ri / len(floaters), start_col, col_size)
                    total += weight * ri
            branch_rates['/'.join(path)] = total
        else:
            cs = read_json(os.path.join(dirname, 'cluster-set.json'))
            for rate_branch in read_json(os.path.join(dirname, 'rate-tree.json')):
                cluster_id = len(clusters)
                w = round(weight * rate_branch['weight'], 8)
                clusters.append((cluster_id, rate_branch['value'], w, '/'.join(path) + ':' + rate_branch['id']))
                for rs in cs['rupture-sets']:
                    for m, rates, mfd_type in mfd_branches(mfd_map[rs['mfd-tree']]):
                        assert mfd_type == 'SINGLE'
                        # rate holds the magnitude branch weight, as in nshmp-lib cluster calculations
                        cluster.add(rs['sections'], m[0], rates[0], extra=(cluster_id, n_cluster_sections))
                    n_cluster_sections += 1

    data = write_model(interface, output_dir)
    data_cluster = write_model(cluster, output_dir_cluster, extra_names=('cluster_id', 'cluster_section'))
    np.savez_compressed(
        os.path.join(output_dir_cluster, 'clusters.npz'),
        cluster_id=np.array([c[0] for c in clusters], dtype=np.int32),
        rate=np.array([c[1] for c in clusters]),
        weight=np.array([c[2] for c in clusters]))

    write_source_info(data, data_cluster, clusters, branch_rates)
    print(f'interface: {len(data["m"])} ruptures, {interface.n_segments} segments, total rate {data["rate"].sum():.6g}')
    print(f'cluster: {len(clusters)} clusters, {len(data_cluster["m"])} ruptures, {cluster.n_segments} segments')


def write_source_info(data, data_cluster, clusters, branch_rates):
    common_notes = (
        'Converted from nshm-conus release 6.2.0 (identical to main) by utilities/convert_nshm23_cascadia_interface.py '
        'following nshmp-lib commit 44728a7d (InterfaceRuptureSet, ApproxGriddedSurface, RuptureFloating, '
        'RuptureScaling, ModelLoader). The source logic tree is collapsed to mean rates: rupture rates include all '
        'source logic tree branch weights and source-group scale factors (geometry branches top 0.2, middle 0.5, '
        'bottom 0.3 are separate ruptures). '
        'Segments are the strips between adjacent columns of the nshmp-lib 5 km gridded surface of each rupture set; '
        'nshmp-lib computes rRup and rJB as minimum distances to the grid points, rX from the upper edge of the rupture. '
        'GMM inputs as used by nshmp-lib: ztor, dip, width (ruptures.npz "width"; zbor = ztor + width sin(dip)), '
        'zhyp = ztor + width/2 sin(dip), rake = 90. fault_type = 1 (reverse) because these are interface events. '
        'GMM tree (subduction/interface/gmm-tree.json): AM_09_INTERFACE_BASIN 0.125, ZHAO_06_INTERFACE_BASIN 0.125, '
        'AG_20_CASCADIA_INTERFACE_ADJUSTED_BASIN 0.25, KBCG_20_CASCADIA_INTERFACE_BASIN 0.25, '
        'PSBAH_20_CASCADIA_INTERFACE_BASIN 0.25; gmm-config.json max-distance 1000 km.')
    info = {
        'name': 'nshm23_cascadia_interface',
        'tectonic_region': 'subduction_interface',
        'nshm_component': 'Interface',
        'source': 'nshm-conus 6.2.0 subduction/interface (Cascadia: top, middle, bottom; full-rupture/cluster-out '
                  'and partial-rupture branches)',
        'cluster': False,
        'interface_events': True,
        'rupture_surface': {
            'type': 'ApproxGriddedSurface between joined upper and lower section traces',
            'surface_spacing_km': 5.0,
            'single_mfd_floating': 'OFF (full rupture)',
            'gr_mfd_floating': 'STRIKE_ONLY, full down-dip width, length 10^((M - 4.94)/1.39) km '
                               '(NSHM_SUB_GEOMAT_LENGTH), rate divided equally among floaters',
            'distance_cutoff_km': 1000.0,
            'distances': 'nshmp-lib Distance.compute: rRup and rJB are minimum distances from the site to the grid '
                         'points of the rupture (rJB = 0 if within 2.5 km and inside the perimeter); rX from the upper '
                         'edge. ucla_plha triangles of the strips between grid columns reproduce rRup within about 0.5 km.',
            'gmm_inputs': 'ruptures.npz: m, dip (average plunge of the first and last columns), ztor (mean depth of the '
                          'upper trace vertices for full ruptures, of the top row grid points for floating ruptures), '
                          'width ((rows - 1) x 5 km), zhyp = ztor + width/2 sin(dip), rake = 90; zbor = ztor + width sin(dip)',
        },
        'total_rate': float(data['rate'].sum()),
        'branch_total_rates': branch_rates,
        'notes': common_notes + ' The full-rupture cluster-in branch (weight 0.1) is in nshm23_cascadia_interface_cluster '
                 '(hazard service component FaultCluster). ruptures.npz extra keys: width (km, nshmp-lib down-dip width), '
                 'zhyp (km), rake (deg).',
    }
    json.dump(info, open(os.path.join(output_dir, 'source_info.json'), 'w', encoding='utf-8'), indent=2, ensure_ascii=False)

    info_cluster = {
        'name': 'nshm23_cascadia_interface_cluster',
        'tectonic_region': 'subduction_interface',
        'nshm_component': 'FaultCluster',
        'source': 'nshm-conus 6.2.0 subduction/interface/Cascadia/{top,middle,bottom}/full-rupture/cluster-in '
                  '(cluster-7a, cluster-7b, cluster-8)',
        'cluster': True,
        'interface_events': True,
        'clusters': [{'cluster_id': c[0], 'rate': c[1], 'weight': c[2], 'branch': c[3]} for c in clusters],
        'cluster_hazard': (
            'nshmp-lib ClusterRuptureSet (source type FAULT_CLUSTER, reported in the FaultCluster component) and '
            'calc.Transforms.ClusterGroundMotionsToCurves. Each cluster (clusters.npz cluster_id) is a sequence of '
            'rupture sets (cluster_section) that all rupture when the cluster occurs, with annual rate clusters.npz '
            '"rate" (0.0019) and logic tree weight clusters.npz "weight" (geometry weight x 0.1 x cluster-set weight). '
            'Each cluster section has magnitude variants (M1..M4) whose ruptures.npz "rate" holds the magnitude branch '
            'weight (0.25; equal magnitudes in a section are merged by adding weights). For each GMM g (and each of its '
            'epistemic branches, combined after), the probability that a section s causes exceedance is '
            'P_s(x) = sum_{r in s} rate_r * P_g(Y > x | r); the cluster exceedance probability is '
            'P_c(x) = 1 - prod_s (1 - P_s(x)); the annual rate of exceedance is '
            'lambda(x) = sum_c weight_c * rate_c * sum_g w_g * P_c,g(x). The product over sections is taken per GMM, '
            'not after combining GMMs. Ruptures are not independent Poissonian events.'),
        'notes': common_notes + ' ruptures.npz "rate" is the magnitude branch weight within a cluster section, not an '
                 'annual rate. ruptures.npz extra keys: cluster_id, cluster_section (unique over all clusters), width, '
                 'zhyp, rake.',
    }
    json.dump(info_cluster, open(os.path.join(output_dir_cluster, 'source_info.json'), 'w', encoding='utf-8'), indent=2, ensure_ascii=False)


if __name__ == '__main__':
    main()

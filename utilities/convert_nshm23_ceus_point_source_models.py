'''
Convert the stable-crust (central and eastern U.S.) zone and gridded seismicity sources of the USGS
NSHM conterminous U.S. model (nshm-conus release 6.2.0) into ucla_plha point source models:

    nshm23_ceus_zone         stable-crust/zone: AR Crowley's Ridge (south, west), Joiner Ridge,
                             Marianna, Saline River; IL Wabash Valley; ME Charlevoix; SC Charleston
                             (regional, local, narrow); VA Central Virginia (regional, local)
    nshm23_ceus_grid         stable-crust/grid/ceus-stable: smoothed seismicity, USGS and SSCn
                             Mmax-zone branches x 3 declustering x 2 smoothing spatial PDFs x
                             12 (Mmax, b-value) GR branches per Mmax zone
    nshm23_ceus_grid_system  stable-crust/grid/system-stable: the part of the WUS fault-system
                             branch-averaged gridded seismicity (active-crust/grid/grid-data/
                             branch-avg-grid.csv) that falls in the stable-crust polygon
                             (grid-system-stable.geojson; lon -116 to -104). nshmp-lib uses a mixed
                             NGA-East / NGA-West2 GMM tree for it (see its source_info.json).

Logic trees are collapsed to mean rates: the rate of each (node, magnitude) is the sum over source
tree leaves of leaf weight x MFD branch weight x MFD rate, exactly as nshmp-lib GridLoader and
ZoneRuptureSet build node MFDs. All three models are strike-slip only (focal-mech-tree STRIKE_SLIP
1.0) except the system grid, which uses the ss/r/n weights in branch-avg-grid.csv.

Run from the utilities directory: python convert_nshm23_ceus_point_source_models.py
Input files are downloaded into utilities/nshm23_ceus/ (git ignored).
'''
import json
import os

import numpy as np
import pandas as pd

import nshm23_ceus_common as nc

OUT = '../src/ucla_plha/source_models/point_source_models'
NSHMP_LIB_RESOURCES = 'https://code.usgs.gov/api/v4/projects/1356/repository/files/src%2Fmain%2Fresources%2Ffault%2Fsurface%2F{}/raw?ref=' + nc.NSHMP_LIB_COMMIT

# Grid MFD bins with a mean annual rate below this value are not written. The effect on the total
# rate of the grid model is reported when the script runs (it is < 1e-6 of the total).
GRID_RATE_FLOOR = 1e-12


def awt_contains(ring, lon, lat):
    '''
    java.awt.geom.Area.contains(x, y) for a simple polygon, as used by SourceFeature.Grid.contains().
    The Area is built from the polygon's straight edges (horizontal edges are dropped); a point is
    inside if it is inside the half-open bounding box [xmin, xmax) x [ymin, ymax) and the number of
    edges crossed by a ray towards -x (Curve.crossingsFor) is odd. This reproduces nshmp-lib's
    assignment of grid nodes that lie on, or within rounding error of, a polygon edge.
    '''
    ring = np.asarray(ring, dtype=float)
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    inside = ((lon >= ring[:, 0].min()) & (lat >= ring[:, 1].min()) &
              (lon < ring[:, 0].max()) & (lat < ring[:, 1].max()))
    crossings = np.zeros(len(lon), dtype=int)
    for (xa, ya), (xb, yb) in zip(ring[:-1], ring[1:]):
        if ya == yb:
            continue
        x0, y0, x1, y1 = (xa, ya, xb, yb) if ya < yb else (xb, yb, xa, ya)
        k = (lat >= y0) & (lat < y1) & (lon < max(x0, x1))
        if x0 == x1:
            xfy = np.full(len(lat), x0)
        else:
            xfy = np.where(lat <= y0, x0, np.where(lat >= y1, x1, x0 + (lat - y0) * (x1 - x0) / (y1 - y0)))
        k &= (lon < min(x0, x1)) | (lon < xfy)
        crossings += k
    return inside & (crossings % 2 == 1)


def polygon_ring(geometry):
    assert geometry['type'] == 'Polygon' and len(geometry['coordinates']) == 1
    ring = [tuple(map(float, p[:2])) for p in geometry['coordinates'][0]]
    if ring[0] != ring[-1]:
        ring.append(ring[0])
    return ring


def write_points(output_dir, lon, lat):
    os.makedirs(output_dir, exist_ok=True)
    xyz = nc.lat_lon_depth_to_xyz(lat, lon, np.zeros(len(lat)))
    np.save(os.path.join(output_dir, 'points.npy'), xyz)
    np.save(os.path.join(output_dir, 'node_index.npy'), np.arange(len(lat), dtype=np.int32))


def read_rjb_table(name):
    '''RuptureScaling.readRjb(): 26 magnitudes (6.05..8.55) x 1001 distances (0..1000 km).'''
    local = os.path.join(nc.INPUT_DIR, 'nshmp-lib', name)
    if not os.path.exists(local):
        os.makedirs(os.path.dirname(local), exist_ok=True)
        nc._curl(NSHMP_LIB_RESOURCES.format(name), local)
    tables = []
    rows = None
    for line in open(local):
        line = line.strip()
        if not line:
            continue
        if line.startswith('#Mag'):
            rows = []
            tables.append(rows)
            continue
        if line.startswith('#'):
            continue
        rows.append(float(line.split()[1]))
    rjb = np.array(tables)
    assert rjb.shape == (26, 1001)
    return rjb


def save_rjb_table(output_dir, name):
    np.savez_compressed(os.path.join(output_dir, 'rjb_correction.npz'),
                        m=np.round(6.05 + 0.1 * np.arange(26), 2), r=np.arange(1001.0), rjb=read_rjb_table(name))


# ------------------------------------------------------------------------------------------------
# Zones
# ------------------------------------------------------------------------------------------------

def zone_leaves(directory, weight, out):
    '''Walk a zone directory: source-tree.json branches (null = no source) or a single geojson.'''
    tree = nc.read_optional_json(directory, 'source-tree.json')
    if tree is not None:
        for b in tree:
            if b['id'] is not None:
                zone_leaves(os.path.join(directory, b['id']), weight * b['weight'], out)
        return
    geojsons = [n for n in sorted(os.listdir(directory)) if n.endswith('.geojson')]
    assert len(geojsons) == 1, directory
    out.append((os.path.join(directory, geojsons[0]), nc.jround(weight, 8)))


def convert_zones():
    root = nc.fetch_tree('stable-crust/zone')
    config = nc.read_json(os.path.join(root, 'zone-config.json'))
    leaves = []
    for state in sorted(os.listdir(root)):
        sdir = os.path.join(root, state)
        if not os.path.isdir(sdir):
            continue
        for name in sorted(os.listdir(sdir)):
            zdir = os.path.join(sdir, name)
            if os.path.isdir(zdir):
                start = len(leaves)
                zone_leaves(zdir, 1.0, leaves)
                for i in range(start, len(leaves)):
                    leaves[i] = leaves[i] + (zdir,)
    nodes = {}   # (lon, lat, strike) -> node index
    node_list = []
    rup = {}     # (node, m) -> rate
    zones = []
    for geojson, weight, zdir in leaves:
        feature = nc.read_json(geojson)
        props = feature['properties']
        mfd_map = nc.read_optional_json(zdir, 'mfd-map.json')
        tree = nc.mfd_tree(props['mfd-tree'], mfd_map)
        assert all(p['type'] == 'SINGLE' for _, p, _ in tree)
        strike = float(props['strike']) if 'strike' in props else np.nan
        rates = pd.read_csv(os.path.join(os.path.dirname(geojson), props['rate']))
        total = 0.0
        for lon, lat, node_rate in rates[['lon', 'lat', 'rate']].itertuples(index=False):
            key = (float(lon), float(lat), strike)
            if key not in nodes:
                nodes[key] = len(node_list)
                node_list.append(key)
            inode = nodes[key]
            # ZoneRuptureSet: node MFD = (weighted SINGLE MFD tree) scaled to the node rate
            for _, p, w in tree:
                r = node_rate * p['rate'] * w * weight
                rup[(inode, p['m'])] = rup.get((inode, p['m']), 0.0) + r
                total += r
        zones.append({'name': props['name'], 'id': feature['id'], 'path': os.path.relpath(geojson, root).replace(os.sep, '/'),
                      'leaf_weight': weight, 'strike': None if np.isnan(strike) else strike,
                      'point_source_type': 'FINITE' if np.isnan(strike) else 'FIXED_STRIKE',
                      'n_nodes': len(rates), 'sum_node_rates': float(rates['rate'].sum()),
                      'mfd_tree': [(i, p['m'], w) for i, p, w in tree], 'weighted_total_rate': total})
    out = os.path.join(OUT, 'nshm23_ceus_zone')
    lon = np.array([k[0] for k in node_list])
    lat = np.array([k[1] for k in node_list])
    write_points(out, lon, lat)
    np.save(os.path.join(out, 'strike.npy'), np.array([k[2] for k in node_list]))
    keys = sorted(rup)
    np.savez_compressed(os.path.join(out, 'ruptures.npz'),
                        node_index=np.array([k[0] for k in keys], dtype=np.int32),
                        m=np.array([k[1] for k in keys]),
                        style=np.full(len(keys), 3, dtype=np.int32),
                        rate=np.array([rup[k] for k in keys]))
    save_rjb_table(out, 'rjb_wc94length.dat')
    with open(os.path.join(out, 'zones.json'), 'w') as f:
        json.dump(zones, f, indent=1)
    print('nshm23_ceus_zone: %d nodes, %d ruptures, total rate %.5e' % (len(node_list), len(keys), sum(rup.values())))
    return config, zones


# ------------------------------------------------------------------------------------------------
# Gridded seismicity (ceus-stable)
# ------------------------------------------------------------------------------------------------

def grid_leaves(directory, weight, out):
    tree = nc.read_optional_json(directory, 'source-tree.json')
    if tree is not None:
        for b in tree:
            if b['id'] is not None:
                grid_leaves(os.path.join(directory, b['id']), weight * b['weight'], out)
        return
    # rupture-sets.json: a LogicGroup of rupture sets (one per Mmax zone), each with weight 1.0
    for rs in nc.read_json(os.path.join(directory, 'rupture-sets.json')):
        out.append((rs, nc.jround(weight, 8), directory))


def gr_node_shape(mfd_props_tree, mags):
    '''
    GridLoader.createNodeMfdTree() without a rate-tree: each GR branch is 10^(a - b m) on the common
    magnitude array of the tree (bins above the branch mMax are zero), scaled by the node pdf. Returns
    sum_b w_b MFD_b(m) for pdf = 1.
    '''
    shape_ = np.zeros(len(mags))
    for _, p, w in mfd_props_tree:
        rates = 10.0 ** (p['a'] - p['b'] * mags)
        rates[mags > p['mMax']] = 0.0
        shape_ += w * rates
    return shape_


def convert_grid():
    root = nc.fetch_tree('stable-crust/grid')
    config = nc.read_json(os.path.join(root, 'grid-config.json'))
    features = {}
    for name in os.listdir(os.path.join(root, 'features')):
        d = nc.read_json(os.path.join(root, 'features', name))
        features[d['id']] = (d['properties']['name'], polygon_ring(d['geometry']))
    ceus = os.path.join(root, 'ceus-stable')
    mfd_map = nc.read_json(os.path.join(ceus, 'mfd-map.json'))
    leaves = []
    grid_leaves(ceus, 1.0, leaves)

    # Common magnitude array: all ceus-stable MFD trees use mMin 4.7, Δm 0.1, max mMax 8.0
    props = [nc.mfd_tree(k, mfd_map) for k in mfd_map]
    mmin = {p['mMin'] for t in props for _, p, _ in t}
    dm = {p['Δm'] for t in props for _, p, _ in t}
    mmax = max(p['mMax'] for t in props for _, p, _ in t)
    assert len(mmin) == 1 and len(dm) == 1
    mags = nc.sequence(mmin.pop(), mmax, dm.pop(), centered=False)

    pdfs = {}
    for f in sorted({rs['spatial-pdf'] for rs, _, _ in leaves}):
        pdfs[f] = pd.read_csv(nc.fetch('stable-crust/grid/grid-data/' + f))
    # node set: union of all pdf nodes (keyed by the lon, lat strings rounded to 0.01 degree)
    allnodes = pd.concat([d[['lon', 'lat']] for d in pdfs.values()]).round(2).drop_duplicates()
    allnodes = allnodes.sort_values(['lat', 'lon']).reset_index(drop=True)
    index = {k: i for i, k in enumerate(zip(allnodes['lon'], allnodes['lat']))}
    n = len(allnodes)
    rate = np.zeros((n, len(mags)))
    summary = []
    membership = {}
    for rs, weight, directory in leaves:
        name, poly = features[rs['feature']]
        df = pdfs[rs['spatial-pdf']]
        key = (rs['feature'], rs['spatial-pdf'])
        if key not in membership:
            membership[key] = awt_contains(poly, df['lon'].values, df['lat'].values)
        inside = membership[key]
        idx = np.array([index[k] for k in zip(df['lon'].round(2)[inside], df['lat'].round(2)[inside])], dtype=int)
        node_shape = gr_node_shape(nc.mfd_tree(rs['mfd-tree'], mfd_map), mags)
        contrib = weight * np.outer(df['pdf'].values[inside], node_shape)
        np.add.at(rate, idx, contrib)
        summary.append({'name': rs['name'], 'id': rs['id'], 'feature': rs['feature'], 'feature_name': name,
                        'mfd_tree': rs['mfd-tree'], 'spatial_pdf': rs['spatial-pdf'], 'leaf_weight': weight,
                        'n_nodes': int(inside.sum()), 'pdf_sum': float(df['pdf'].values[inside].sum()),
                        'weighted_rate_m_ge_4p7': float(contrib.sum()),
                        'weighted_rate_m_ge_5': float(contrib[:, mags >= 5].sum())})
    # every pdf node is in exactly one Mmax zone of each group
    for f, df in pdfs.items():
        for group in ('USGS', 'SSCn'):
            count = sum(membership[(fid, f)].astype(int) for (fid, ff) in membership if ff == f
                        and features[fid][0].startswith(group))
            assert np.all(count == 1), (f, group, np.unique(count, return_counts=True))

    out = os.path.join(OUT, 'nshm23_ceus_grid')
    total = rate.sum()
    keep = rate >= GRID_RATE_FLOOR
    dropped = rate[~keep].sum()
    used = np.unique(np.nonzero(keep)[0])
    new_index = -np.ones(n, dtype=np.int64)
    new_index[used] = np.arange(len(used))
    write_points(out, allnodes['lon'].values[used], allnodes['lat'].values[used])
    inode, imag = np.nonzero(keep)
    np.savez_compressed(os.path.join(out, 'ruptures.npz'),
                        node_index=new_index[inode].astype(np.int32),
                        m=mags[imag],
                        style=np.full(len(inode), 3, dtype=np.int32),
                        rate=rate[inode, imag].astype(np.float32))
    save_rjb_table(out, 'rjb_somerville.dat')
    with open(os.path.join(out, 'rupture_sets.json'), 'w') as f:
        json.dump(summary, f, indent=1)
    print('nshm23_ceus_grid: %d pdf nodes, %d nodes written, %d ruptures, total rate %.6e (M>=5 %.6e); '
          'dropped %.3e (%.2e of total) below %.0e' % (n, len(used), len(inode), total, rate[:, mags >= 5].sum(),
                                                         dropped, dropped / total, GRID_RATE_FLOOR))
    return config, summary, {'total_rate': float(total), 'dropped_rate': float(dropped),
                             'n_pdf_nodes': n, 'n_nodes': int(len(used)), 'n_ruptures': int(len(inode)),
                             'float32_max_relative_error': float(np.max(np.abs(
                                 rate[inode, imag].astype(np.float32) / rate[inode, imag] - 1)))}


def convert_grid_system():
    root = nc.fetch_tree('stable-crust/grid')
    sysdir = os.path.join(root, 'system-stable')
    rs = nc.read_json(os.path.join(sysdir, 'branch-avg', 'rupture-sets.json'))
    assert len(rs) == 1
    rs = rs[0]
    feature = nc.read_json(os.path.join(root, 'features', 'grid-system-stable.geojson'))
    assert feature['id'] == rs['feature']
    poly = polygon_ring(feature['geometry'])
    mfd = nc.read_json(os.path.join(sysdir, 'mfd-map.json'))[rs['mfd-tree']]
    assert len(mfd) == 1 and mfd[0]['weight'] == 1.0
    df = pd.read_csv(nc.fetch('active-crust/grid/grid-data/branch-avg-grid.csv'))
    inside = awt_contains(poly, df['lon'].values, df['lat'].values)
    df = df[inside].reset_index(drop=True)
    m_columns = [repr(m) for m in mfd[0]['value']['magnitudes']]
    mags = np.array(mfd[0]['value']['magnitudes'])
    out = os.path.join(OUT, 'nshm23_ceus_grid_system')
    write_points(out, df['lon'].values, df['lat'].values)
    rates = df[m_columns].values
    node = np.repeat(np.arange(len(df)), len(mags))
    m = np.tile(mags, len(df))
    r = rates.reshape(-1)
    all_node, all_m, all_style, all_rate = [], [], [], []
    for style, col in [(3, 'ss_wt'), (1, 'r_wt'), (2, 'n_wt')]:  # 1 reverse, 2 normal, 3 strike slip
        rr = np.repeat(df[col].values, len(mags)) * r
        k = rr > 0
        all_node.append(node[k]); all_m.append(m[k]); all_style.append(np.full(k.sum(), style)); all_rate.append(rr[k])
    np.savez_compressed(os.path.join(out, 'ruptures.npz'), node_index=np.concatenate(all_node).astype(np.int32),
                        m=np.concatenate(all_m), style=np.concatenate(all_style).astype(np.int32),
                        rate=np.concatenate(all_rate))
    save_rjb_table(out, 'rjb_somerville.dat')
    total = np.concatenate(all_rate).sum()
    print('nshm23_ceus_grid_system: %d nodes, total rate %.6e' % (len(df), total))
    return {'n_nodes': int(len(df)), 'total_rate': float(total), 'rupture_set': rs,
            'gmm_tree': nc.read_json(os.path.join(sysdir, 'gmm-tree.json')),
            'gmm_config': nc.read_json(os.path.join(sysdir, 'gmm-config.json'))}


# ------------------------------------------------------------------------------------------------
# source_info.json
# ------------------------------------------------------------------------------------------------

STABLE_GMM_TREE = ('stable-crust/gmm-tree.json: NGA_EAST_2026 0.3333, NGA_EAST_2026_ADJUSTED 0.3333, '
                   'NGA_EAST_SEEDS_2026 0.1667, NGA_EAST_SEEDS_2026_ADJUSTED 0.1667; stable-crust/gmm-config.json '
                   'max-distance 1000 km')
RJB_NOTE = ('rJB = rjb_correction.npz["rjb"][i_m, i_r] if M >= 6.0, else rJB = r_epi, with '
            'i_m = min(round((M - 6.05) / 0.1), 25) and i_r = min(floor(r_epi), 1000) (no interpolation; '
            'RuptureScaling.correctedRjb). The table is nshmp-lib resources/fault/surface/{}: the mean rJB of '
            'a finite rupture of random strike centered on the node.')
EPI_NOTE = ('r_epi is the horizontal node-to-site distance from Locations.horzDistanceFast: '
            'R * sqrt(dlat^2 + (dlon * cos(mean lat))^2), R = 6371.0072 km, angles in radians. '
            '(ucla_plha points.npy are at the surface, so the chord distance used by plha.get_source_data '
            'differs from this by < 0.1% within 300 km.)')
FINITE_SS_NOTE = ('GridSourceFinite, strike-slip (dip 90, rake 0): rRup = hypot(rJB, zTor), rX = -rJB '
                  '(footwall), zTor = 5 km, down-dip width W = RuptureScaling.{}.dimensions(M, (22 - zTor) / '
                  'sin(dip)).width, zBot = zTor + W, zHyp = zTor + W / 2 (Faults.hypocentralDepth).')
RATE_SERVICE = 'USGS NSHM rate web service (https://earthquake.usgs.gov/ws/nshmp/conus-2023/dynamic/rate/lon/lat/r)'


def depth_info(config):
    gdm = config['grid-depth-map']
    return {'depth_map': {k: {'mMin': v['mMin'], 'mMax': v['mMax'],
                              'ztor_km_and_weight': [(b['value'], b['weight']) for b in v['depth-tree']]}
                          for k, v in gdm.items()},
            'max_depth_km': config['max-depth'],
            'depth_array': 'none: every rupture has zTor = 5 km (grid-depth-map "all": M 4.5-10, depth-tree 5 km, '
                           'weight 1.0), so no depth array is stored'}


def write_source_info(output_dir, info):
    with open(os.path.join(output_dir, 'source_info.json'), 'w', encoding='utf-8') as f:
        json.dump(info, f, indent=2, ensure_ascii=False)


def zone_source_info(config, zones):
    return {
        'name': 'nshm23_ceus_zone',
        'tectonic_region': 'stable_crust',
        'nshm_component': 'Zone',
        'source': 'nshm-conus %s stable-crust/zone (AR/Crowleys Ridge (south), AR/Crowleys Ridge (west), '
                  'AR/Joiner Ridge, AR/Marianna, AR/Saline River, IL/Wabash Valley, ME/Charlevoix, SC/Charleston, '
                  'VA/Central Virginia), zone-config.json' % nc.NSHM_REF,
        'nshm_release': nc.NSHM_REF,
        'nshmp_lib_commit': nc.NSHMP_LIB_COMMIT,
        'cluster': False,
        'gmm_tree': STABLE_GMM_TREE,
        'point_source': {
            'point_source_type': 'FIXED_STRIKE for nodes with a finite strike in strike.npy (all zones except '
                                 'Charlevoix; the strike is the zone feature "strike" property); FINITE for nodes '
                                 'with strike NaN (ME Charlevoix has no strike).',
            'rupture_scaling': config['rupture-scaling'],
            'rupture_dimensions': 'NSHM_POINT_WC94_LENGTH: L = 10^(-3.22 + 0.69 M) km, '
                                  'W = min((22 - zTor) / sin(dip), L / 1.5)',
            'focal_mech_tree': config['focal-mech-tree'],
            'style': 'all ruptures are strike-slip (style 3), dip 90, rake 0',
            'depth': depth_info(config),
            'magnitudes': 'the SINGLE magnitudes of each zone mfd-tree (zones.json); no aleatory or epistemic '
                          'magnitude variability',
            'distances_fixed_strike': 'GridSourceFixedStrike: a vertical rupture whose top trace runs from '
                                      'p1 to p2, the points L/2 from the node along azimuths strike and strike+180 '
                                      '(spherical, Locations.location), at depth zTor = 5 km, bottom zBot = zTor + W. '
                                      'rJB = Locations.distanceToSegmentFast(p1, p2, site) (longitudes scaled by the '
                                      'cosine of 0.5 lat_site + 0.25 lat_p1 + 0.25 lat_p2, then Line2D.ptSegDist, '
                                      'times R = 6371.0072 km), rRup = hypot(rJB, zTor), '
                                      'rX = Locations.distanceToLineFast(p1, p2, site). No distance correction.',
            'distances_finite': 'Charlevoix nodes (strike NaN): ' + RJB_NOTE.format('rjb_wc94length.dat') + ' ' +
                                FINITE_SS_NOTE.format('NSHM_POINT_WC94_LENGTH'),
            'r_epi': EPI_NOTE,
            'optimization': 'none: nshmp-lib builds distance-binned rate tables only for GRID rupture sets; zone '
                            'nodes are used individually.',
        },
        'files': {
            'ruptures.npz': 'node_index, m, style (3 = strike slip), rate (annual, mean over the logic tree)',
            'points.npy': 'node xyz (geometry.point_to_xyz convention) at the surface',
            'node_index.npy': '0..N-1',
            'strike.npy': 'per node (aligned with node_index.npy): fixed strike in degrees, NaN for FINITE nodes',
            'rjb_correction.npz': 'RuptureScaling RJB_WC94LENGTH table (m 6.05..8.55, r 0..1000 km, rjb), '
                                  'used only for the FINITE (Charlevoix) nodes',
            'zones.json': 'per zone source: path, NSHM id, leaf weight, strike, mfd tree, node count, rates',
        },
        'notes': 'Each zone is a list of nodes with annual rates (the zone csv files). nshmp-lib ZoneRuptureSet '
                 'scales the weighted SINGLE mfd-tree of the zone to each node rate, so rate(node, M_i) = '
                 'node_rate x w_i x leaf_weight, with leaf_weight the product of source-tree branch weights '
                 '("active" 0.5 vs null 0.5 for the AR zones and Central Virginia; Charleston regional 0.2, '
                 'local 0.5, narrow 0.3; Central Virginia regional 0.5, local 0.5; Wabash Valley and Charlevoix '
                 '1.0). Nodes shared by zones with the same strike are merged. Validated against the ' +
                 RATE_SERVICE + ' Zone component at Charleston, the AR zones, Wabash Valley, Charlevoix, Central '
                 'Virginia, Memphis, and within 1000 km of Ohio: total and per-0.1-magnitude-bin rates agree to '
                 '< 1e-5.',
        'zones': [{k: z[k] for k in ('name', 'id', 'leaf_weight', 'strike', 'point_source_type', 'n_nodes')}
                  for z in zones],
    }


def grid_source_info(config, summary, stats):
    return {
        'name': 'nshm23_ceus_grid',
        'tectonic_region': 'stable_crust',
        'nshm_component': 'Grid',
        'source': 'nshm-conus %s stable-crust/grid/ceus-stable (usgs and sscn Mmax-zone branches), '
                  'grid-data/pdf-{gk,nn,r85}-{adaptive,fixed}.csv, features/ (Mmax zone polygons), '
                  'ceus-stable/mfd-map.json, grid-config.json' % nc.NSHM_REF,
        'nshm_release': nc.NSHM_REF,
        'nshmp_lib_commit': nc.NSHMP_LIB_COMMIT,
        'cluster': False,
        'gmm_tree': STABLE_GMM_TREE,
        'point_source': {
            'point_source_type': config['point-source-type'] + ' (GridSourceFinite)',
            'rupture_scaling': config['rupture-scaling'],
            'rupture_dimensions': 'NSHM_SOMERVILLE: A = 10^(M - 4.366) km^2; W = sqrt(A) if sqrt(A) < maxWidth '
                                  'else maxWidth, with maxWidth = (22 - zTor) / sin(dip) = 17 km',
            'focal_mech_tree': config['focal-mech-tree'],
            'style': 'all ruptures are strike-slip (style 3), dip 90, rake 0',
            'depth': depth_info(config),
            'magnitudes': 'GR bin centers 4.75, 4.85, ..., 7.95 (mMin 4.7, dm 0.1, branch mMax 6.5-8.0; bins '
                          'above a branch mMax have zero rate in that branch)',
            'r_epi': EPI_NOTE,
            'rjb': RJB_NOTE.format('rjb_somerville.dat'),
            'rrup_rx_zhyp': FINITE_SS_NOTE.format('NSHM_SOMERVILLE'),
            'optimization': 'nshmp-lib defaults (calc-config optimizeGrids = true, smoothGrids = true) replace '
                            'the nodes, for each site, by a table of rates in 5 km distance bins '
                            '(opt-distance-bin %s km; rows are bin centers 2.5, 7.5, ... km out to 1000 km; '
                            'r_epi from horzDistanceFast). Nodes closer than smoothing-limit %s km are first split '
                            'into 4 x 4 sub-points at lat/lon offsets of -0.0375, -0.0125, 0.0125, 0.0375 degrees '
                            '(smoothing-density %d, grid spacing %s) with 1/16 of the rate each. Each table row '
                            'becomes a GridSourceFinite at the bin-center distance. ucla_plha can use the nodes '
                            'directly (exact) or emulate the table.' % (
                                config['opt-distance-bin'], config['smoothing-limit'], config['smoothing-density'],
                                config['grid-spacing']),
        },
        'files': {
            'ruptures.npz': 'node_index (int32), m (float64), style (int32, 3 = strike slip), rate (float32, '
                            'annual, mean over the logic tree; max relative rounding error %.1e)' %
                            stats['float32_max_relative_error'],
            'points.npy': 'node xyz (geometry.point_to_xyz convention) at the surface',
            'node_index.npy': '0..N-1',
            'rjb_correction.npz': 'RuptureScaling RJB_SOMERVILLE table (m 6.05..8.55, r 0..1000 km, rjb)',
            'rupture_sets.json': 'per source-tree leaf rupture set: Mmax zone, mfd tree, spatial pdf, leaf weight, '
                                 'node count, weighted rates',
        },
        'logic_tree': 'ceus-stable: usgs 0.5 / sscn 0.5 (Mmax zone polygons and Mmax trees) x declustering gk 0.4, '
                      'nn 0.4, r85 0.2 x smoothing fixed 0.6 / adaptive 0.4; each leaf is a group of rupture sets '
                      '(one per Mmax zone polygon, weight 1.0) sharing the leaf spatial pdf. Each Mmax zone has a '
                      '12-branch GR mfd-tree (4 mMax x 3 (a, b) pairs; a-values are total-catalog rates for '
                      'pdf = 1). nshmp-lib GridLoader node MFD: rate(M) = pdf_node x sum_b w_b 10^(a_b - b_b M) '
                      '[M <= mMax_b]. Mean rate(node, M) = sum over leaves of leaf_weight x that. A node is in a '
                      'polygon per java.awt.geom.Area.contains (reproduced exactly; every pdf node is in exactly '
                      'one polygon of each group).',
        'notes': 'Total annual rate M >= 4.7: %.6e. Bins with a mean rate below %.0e/yr are dropped: %.2e/yr in '
                 'total (%.1e of the model total; low-seismicity offshore and border nodes), leaving %d of %d '
                 'pdf nodes. Validated against the %s Grid component at Charleston, the AR zones, Wabash Valley, '
                 'Charlevoix, Central Virginia, Kansas, Memphis, West Texas, Montana, and within 1000 km of Ohio: '
                 'total and per-0.1-magnitude-bin rates agree to < 1e-5 (adding nshm23_ceus_grid_system and '
                 'nshm23_wus_grid where the service Grid component includes them).' % (
                     stats['total_rate'], GRID_RATE_FLOOR, stats['dropped_rate'],
                     stats['dropped_rate'] / stats['total_rate'], stats['n_nodes'], stats['n_pdf_nodes'],
                     RATE_SERVICE),
        'rupture_sets': [{k: r[k] for k in ('name', 'id', 'feature_name', 'mfd_tree', 'spatial_pdf', 'leaf_weight')}
                         for r in summary],
    }


def grid_system_source_info(config, stats):
    return {
        'name': 'nshm23_ceus_grid_system',
        'tectonic_region': 'stable_crust',
        'nshm_component': 'Grid',
        'source': 'nshm-conus %s stable-crust/grid/system-stable/branch-avg (tree 8059 "CEUS System Grid '
                  'Sources"): rupture set 8051 = active-crust/grid/grid-data/branch-avg-grid.csv nodes inside '
                  'stable-crust/grid/features/grid-system-stable.geojson (feature 8050)' % nc.NSHM_REF,
        'nshm_release': nc.NSHM_REF,
        'nshmp_lib_commit': nc.NSHMP_LIB_COMMIT,
        'cluster': False,
        'gmm_tree': stats['gmm_tree'],
        'gmm_max_distance_km': stats['gmm_config']['max-distance'],
        'gmm_note': 'This rupture set has its own gmm-tree.json and gmm-config.json (system-stable/): 1/3 '
                    'NGA-East (NGA_EAST_2026 0.1111, NGA_EAST_2026_ADJUSTED 0.1111, NGA_EAST_SEEDS_2026 0.05555, '
                    'NGA_EAST_SEEDS_2026_ADJUSTED 0.05555) and 2/3 NGA-West2 with basin terms (ASK_14_BASIN, '
                    'BSSA_14_BASIN, CB_14_BASIN 0.1667, CY_14_BASIN 0.1666), and a 300 km maximum distance '
                    '(the ceus-stable grid uses the stable-crust tree and 1000 km).',
        'point_source': {
            'point_source_type': config['point-source-type'] + ' (GridSourceFinite)',
            'rupture_scaling': config['rupture-scaling'],
            'rupture_dimensions': 'NSHM_SOMERVILLE: A = 10^(M - 4.366) km^2; W = sqrt(A) if sqrt(A) < maxWidth '
                                  'else maxWidth, with maxWidth = (22 - zTor) / sin(dip)',
            'focal_mech': 'from the ss_wt, r_wt, n_wt columns of branch-avg-grid.csv (stored as separate ruptures '
                          'with style 3, 1, 2); FocalMech dips: STRIKE_SLIP 90 (rake 0), REVERSE 50 (rake 90), '
                          'NORMAL 50 (rake -90). In 6.2.0 every node of this model has ss_wt = 1, so all '
                          'ruptures are strike-slip; the reverse/normal rules below are given for completeness.',
            'depth': depth_info(config),
            'magnitudes': 'INCR bins 5.05, 5.15, ..., 7.85 (system-stable/mfd-map.json)',
            'r_epi': EPI_NOTE,
            'rjb': RJB_NOTE.format('rjb_somerville.dat'),
            'rrup_rx_zhyp': 'GridSourceFinite: strike-slip as in nshm23_ceus_grid. Reverse and normal ruptures are '
                            'split into two equally weighted (0.5) representations: footwall (rX = -rJB, '
                            'rRup = hypot(rJB, zTor)) and hanging wall (rX = rJB + W cos(dip); with '
                            'rCut = zBot tan(dip): rRup = hypot(rJB, zBot) if rJB > rCut, else linear in rJB from '
                            'rRup0 = min(hypot(W cos(dip), zTor), zBot cos(dip)) at rJB = 0 to zBot / cos(dip) at '
                            'rCut), zBot = zTor + W sin(dip), zTor = 5 km.',
            'optimization': 'as nshm23_ceus_grid (5 km distance-bin tables with 4 x 4 smoothing within 40 km), with '
                            'separate strike-slip, reverse, and normal tables (gridFocalMechUpdate = true).',
        },
        'files': {
            'ruptures.npz': 'node_index, m, style (1 reverse, 2 normal, 3 strike slip), rate (annual)',
            'points.npy': 'node xyz (geometry.point_to_xyz convention) at the surface',
            'node_index.npy': '0..N-1',
            'rjb_correction.npz': 'RuptureScaling RJB_SOMERVILLE table (m 6.05..8.55, r 0..1000 km, rjb)',
        },
        'notes': 'The WUS fault-system (inversion) branch-averaged gridded seismicity covers lon -125 to -104; '
                 'nshm-conus applies the nodes inside grid-system-active.geojson as active-crust grid sources '
                 '(ucla_plha nshm23_wus_grid) and the %d nodes inside grid-system-stable.geojson (lon -116 to '
                 '-104: Colorado Plateau, Rocky Mountains, western Great Plains) as this stable-crust source with '
                 'the mixed GMM tree. Node membership follows java.awt.geom.Area.contains; with that rule the two '
                 'polygons partition the 68883 csv nodes exactly. Total annual rate M >= 5: %.6e. Validated '
                 'against the %s Grid component in Colorado, West Texas, and Montana.' % (
                     stats['n_nodes'], stats['total_rate'], RATE_SERVICE),
    }


if __name__ == '__main__':
    zone_config, zones = convert_zones()
    write_source_info(os.path.join(OUT, 'nshm23_ceus_zone'), zone_source_info(zone_config, zones))
    grid_config, summary, stats = convert_grid()
    write_source_info(os.path.join(OUT, 'nshm23_ceus_grid'), grid_source_info(grid_config, summary, stats))
    sys_stats = convert_grid_system()
    write_source_info(os.path.join(OUT, 'nshm23_ceus_grid_system'), grid_system_source_info(grid_config, sys_stats))

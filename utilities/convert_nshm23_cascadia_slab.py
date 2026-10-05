'''
Convert the Cascadia intraslab gridded seismicity of the USGS NSHM for the conterminous U.S.
(nshm-conus release 6.2.0, subduction/slab) into the Numpy files used by ucla_plha
(point_source_models/nshm23_cascadia_slab). The logic tree is collapsed to mean rates.

Input files are downloaded from https://code.usgs.gov/ghsc/nshmp/nshms/nshm-conus (tag 6.2.0) into
utilities/nshm23_cascadia_slab/ if they are not already there, together with the rJB correction
table rjb_geomatrix.dat from nshmp-lib (src/main/resources/fault/surface, commit 44728a7d).

The conversion follows nshmp-lib (commit 44728a7d) ModelLoader.Slab, GridLoader.SlabLoader,
SlabRuptureSet, and GridSourceFinite:
    - Source trees WA, OR (source-group: mfd-hi and mfd-lo, additive, scale 1) and CA (single-branch).
    - Each branch has a rate-tree (one branch, weight 1) whose value R scales the spatial PDF
      grid-data/pdf-{wa,or,ca}.csv (lon, lat, depth, pdf; the pdf sums to 1 in each state).
    - The node MFD on each mfd-tree branch (mfd-map.json) is a GR distribution with magnitudes at the
      bin centers from mMin + dm/2 to mMax - dm/2, scaled so that the incremental rate of the first bin is
      R * pdf * 10^(-b * (mMin + dm/2)), i.e. rate(M) = R * pdf * 10^(-b M). The MFD branches are combined
      with their weights (all branches share the magnitude bins up to the largest mMax, with zero rates
      above each branch's mMax).
    - slab-config.json: point-source-type FINITE, focal-mech-tree STRIKE_SLIP 1.0, rupture-scaling
      NSHM_SUB_GEOMAT_LENGTH, max-width 8 km, grid-spacing 0.1 deg. The depth to top of rupture is the
      node depth from the pdf file (GridSource.isSubductionGrid), with weight 1.
    - Nodes are kept separate for each state (the WA, OR, and CA grids overlap in a few places).

Distances used by nshmp-lib for these sources (GridSourceFinite.FiniteSurface, strike-slip so all
ruptures are "footwall" ruptures):
    r = Locations.horzDistanceFast(node, site) (spherical, R = 6371.0072 km)
    rJB = rjb_geomatrix[round((M - 6.05) / 0.1) clipped to 0..25, min(1000, floor(r))] for M >= 6, else r
    rRup = sqrt(rJB^2 + zTor^2), rX = -rJB
    dip = 90, width = 8 km, zTor = node depth, zHyp = zTor + 4 km, rake = 0
Only nodes within the GMM max-distance (gmm-config.json, 300 km) of the site are used.
'''
import json
import os
import urllib.parse
import urllib.request
import numpy as np
import pandas as pd

nshm_ref = '6.2.0'
nshm_project = 6063
nshm_path = 'subduction/slab'
lib_project = 1356
lib_ref = '44728a7d'
input_dir = 'nshm23_cascadia_slab'
output_dir = '../src/ucla_plha/source_models/point_source_models/nshm23_cascadia_slab'

STATES = ['WA', 'OR', 'CA']


def gitlab_get(url):
    # code.usgs.gov rejects the default Python user agent
    request = urllib.request.Request(url, headers={'User-Agent': 'curl/8.0'})
    return urllib.request.urlopen(request).read()


def download_inputs():
    '''
    Download the intraslab source model files and the rJB correction table unless they exist.
    '''
    if os.path.exists(os.path.join(input_dir, 'slab-config.json')):
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
    path = 'src/main/resources/fault/surface/rjb_geomatrix.dat'
    with open(os.path.join(input_dir, 'rjb_geomatrix.dat'), 'wb') as f:
        f.write(gitlab_get(f'https://code.usgs.gov/api/v4/projects/{lib_project}/repository/files/'
                           f'{urllib.parse.quote(path, safe="")}/raw?ref={lib_ref}'))


def read_json(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def read_rjb_table(filename):
    '''
    Port of nshmp-lib RuptureScaling.readRjb: 26 magnitudes (6.05 to 8.55) x 1001 distances (0 to 1000 km).
    '''
    table = np.zeros((26, 1001))
    mag_index = -1
    r_index = 0
    for line in open(filename):
        if not line.strip():
            continue
        if line.startswith('#Mag'):
            mag_index += 1
            r_index = 0
            continue
        if line.startswith('#'):
            continue
        table[mag_index, r_index] = float(line.split()[1])
        r_index += 1
    return table


def gr_magnitudes(mfd):
    dm = mfd['Δm']
    n = int(round((mfd['mMax'] - mfd['mMin']) / dm))
    return np.round(mfd['mMin'] + dm / 2.0 + dm * np.arange(n), 5)


def leaf_branches(state_dir):
    '''
    Return a list of (directory, weight) of the leaves of a state source tree. All branches of the
    current model are source-group branches with scale 1.
    '''
    tree = read_json(os.path.join(state_dir, 'source-group.json'))
    return [(os.path.join(state_dir, b['id']), b.get('scale', 1.0)) for b in tree]


def node_mfds(state_dir, mfd_map, grid_dir):
    '''
    Return the pdf DataFrame and an array of collapsed node rates [node, magnitude] with the magnitude bins.
    '''
    pdf = None
    rates = {}
    branch_info = []
    for leaf, scale in leaf_branches(state_dir):
        rupture_sets = read_json(os.path.join(leaf, 'rupture-sets.json'))
        rate_tree = read_json(os.path.join(leaf, 'rate-tree.json'))
        assert len(rupture_sets) == 1
        rs = rupture_sets[0]
        df = pd.read_csv(os.path.join(grid_dir, rs['spatial-pdf']))
        if pdf is None:
            pdf = df
        assert pdf.equals(df)
        for r_branch in rate_tree:
            for m_branch in mfd_map[rs['mfd-tree']]:
                mfd = m_branch['value']
                assert mfd['type'] == 'GR'
                w = scale * r_branch['weight'] * m_branch['weight']
                m = gr_magnitudes(mfd)
                for mi in m:
                    rates[mi] = rates.get(mi, 0.0) + w * r_branch['value'] * pdf['pdf'].values * 10.0 ** (-mfd['b'] * mi)
                branch_info.append({
                    'branch': os.path.relpath(leaf, input_dir).replace(os.sep, '/') + ':' + r_branch['id'] + ':' + m_branch['id'],
                    'weight': w, 'R': r_branch['value'], 'b': mfd['b'], 'mMin': mfd['mMin'], 'mMax': mfd['mMax'],
                    'total_rate': float(w * r_branch['value'] * np.sum(10.0 ** (-mfd['b'] * m)))})
    m = np.array(sorted(rates))
    return pdf, m, np.stack([rates[mi] for mi in m], axis=1), branch_info


def points_xyz(lat, lon):
    rad = np.pi/180.0
    a = 6378.1370 # Earth's equatorial radius in km
    b = 6356.7523 # Earth's polar radius in km
    r = np.sqrt(((a**2 * np.cos(lat * rad))**2 + (b**2 * np.sin(lat * rad))**2) / ((a * np.cos(lat * rad))**2 + (b * np.sin(lat * rad))**2))
    xyz = np.empty((len(lat), 3))
    xyz[:, 0] = r * np.cos(lat * rad) * np.cos(lon * rad)
    xyz[:, 1] = r * np.cos(lat * rad) * np.sin(lon * rad)
    xyz[:, 2] = r * np.sin(lat * rad)
    return xyz


def main():
    download_inputs()
    config = read_json(os.path.join(input_dir, 'slab-config.json'))
    assert config['point-source-type'] == 'FINITE'
    assert config['rupture-scaling'] == 'NSHM_SUB_GEOMAT_LENGTH'
    assert config['focal-mech-tree'] == [{'id': 'STRIKE_SLIP', 'weight': 1.0}]
    mfd_map = read_json(os.path.join(input_dir, 'mfd-map.json'))
    grid_dir = os.path.join(input_dir, 'grid-data')

    lat, lon, depth, state = [], [], [], []
    node_all, m_all, rate_all = [], [], []
    totals = {}
    branches = []
    n0 = 0
    for st in STATES:
        pdf, m, rates, info = node_mfds(os.path.join(input_dir, st), mfd_map, grid_dir)
        branches += info
        n = len(pdf)
        lat.append(pdf['lat'].values)
        lon.append(pdf['lon'].values)
        depth.append(pdf['depth'].values)
        state.append(np.full(n, st))
        node_all.append(np.repeat(np.arange(n0, n0 + n), len(m)))
        m_all.append(np.tile(m, n))
        rate_all.append(rates.reshape(-1))
        totals[st] = {f'M>={mm}': float(rates[:, m >= mm].sum()) for mm in [5.0, 6.0, 7.0, 7.5]}
        n0 += n
    lat, lon, depth, state = (np.concatenate(a) for a in (lat, lon, depth, state))
    node = np.concatenate(node_all).astype(np.int32)
    m = np.concatenate(m_all)
    rate = np.concatenate(rate_all)
    keep = rate > 0
    node, m, rate = node[keep], m[keep], rate[keep]

    os.makedirs(output_dir, exist_ok=True)
    np.save(os.path.join(output_dir, 'points.npy'), points_xyz(lat, lon))
    np.save(os.path.join(output_dir, 'node_index.npy'), np.arange(len(lat), dtype=np.int64))
    # Per-node data: depth (km) is the nshmp-lib depth to top of rupture; lon, lat (degrees)
    np.save(os.path.join(output_dir, 'depth.npy'), depth)
    np.save(os.path.join(output_dir, 'lonlat.npy'), np.stack((lon, lat), axis=1))
    # fault_type / style: 3 = strike slip (slab-config focal-mech-tree STRIKE_SLIP 1.0)
    np.savez_compressed(
        os.path.join(output_dir, 'ruptures.npz'),
        node_index=node, m=m, style=np.full(len(m), 3, dtype=np.int64), rate=rate, depth=depth[node])
    np.save(os.path.join(output_dir, 'rjb_geomatrix.npy'), read_rjb_table(os.path.join(input_dir, 'rjb_geomatrix.dat')))

    gmm = read_json(os.path.join(input_dir, 'gmm-config.json'))
    info = {
        'name': 'nshm23_cascadia_slab',
        'tectonic_region': 'subduction_slab',
        'nshm_component': 'Slab',
        'source': 'nshm-conus 6.2.0 subduction/slab (WA mfd-hi + mfd-lo, OR mfd-hi + mfd-lo, CA single-branch; '
                  'grid-data/pdf-{wa,or,ca}.csv, mfd-map.json, slab-config.json)',
        'cluster': False,
        'point_source': {
            'nshmp_lib_source_type': 'SLAB (SlabRuptureSet, GridSourceFinite)',
            'point_source_type': config['point-source-type'],
            'focal_mechanism': 'STRIKE_SLIP weight 1.0: dip 90, rake 0; style = 3 in ruptures.npz',
            'depth': 'zTor = node depth (km) from the pdf file, weight 1; per node in depth.npy and per rupture in '
                     'ruptures.npz "depth"',
            'width_km': config['max-width'],
            'zhyp': 'zTor + width / 2 = depth + 4 km (Faults.hypocentralDepth with dip 90)',
            'zbor': 'depth + 8 km',
            'rupture_scaling': config['rupture-scaling'],
            'horizontal_distance': 'r = Locations.horzDistanceFast(node, site): R = 6371.0072 km, '
                                   'sqrt(dlat^2 + (dlon cos(mean lat))^2) in radians; node lon, lat in lonlat.npy',
            'rjb': 'for M >= 6.0: rjb_geomatrix.npy[min(round((M - 6.05) / 0.1), 25), min(1000, floor(r))] '
                   '(Java Math.round, i.e. floor(x + 0.5)); for M < 6.0: rJB = r',
            'rrup': 'sqrt(rJB^2 + zTor^2) (footwall branch of GridSourceFinite.FiniteSurface.distanceTo for strike-slip)',
            'rx': '-rJB',
            'magnitudes': 'GR bin centers, dm = 0.1 (5.05..7.95 depending on branch)',
            'distance_cutoff_km': gmm['max-distance'],
            'distance_cutoff_metric': 'horizontal node-to-site distance r (Locations.distanceAndRectangleFilter)',
        },
        'total_rates_by_state': totals,
        'branches': branches,
        'notes': 'Converted from nshm-conus release 6.2.0 (identical to main) by utilities/convert_nshm23_cascadia_slab.py '
                 'following nshmp-lib commit 44728a7d. Mean rates: node rate(M) = sum over branches of weight * R * pdf * '
                 '10^(-b M), where R is the rate-tree value. Nodes are kept separate per state (overlapping grids). '
                 'Extra files: depth.npy (km, per node), lonlat.npy (per node), rjb_geomatrix.npy (nshmp-lib '
                 'rjb_geomatrix.dat, 26 magnitudes 6.05..8.55 x 1001 distances 0..1000 km); extra ruptures.npz key '
                 '"depth" (km, per rupture). GMM tree (subduction/slab/gmm-tree.json): ZHAO_06_SLAB_BASIN 0.25, '
                 'AG_20_CASCADIA_SLAB_BASIN 0.0825, AG_20_CASCADIA_SLAB_ADJUSTED_BASIN 0.1675, '
                 'KBCG_20_CASCADIA_SLAB_BASIN 0.25, PSBAH_20_CASCADIA_SLAB_BASIN 0.25.',
    }
    with open(os.path.join(output_dir, 'source_info.json'), 'w', encoding='utf-8') as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    print(f'{len(lat)} nodes, {len(m)} ruptures')
    for st in STATES:
        print(st, totals[st])


if __name__ == '__main__':
    main()

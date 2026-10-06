'''
Convert the stable-crust (central and eastern U.S.) fault sources of the USGS NSHM conterminous
U.S. model (nshm-conus release 6.2.0, stable-crust/fault) into two ucla_plha fault source models:

    nshm23_ceus_fault          independent fault ruptures: CO/Cheraw, MO/Commerce, OK/Meers,
                               TN/Eastern Rift Margin (north), TN/Eastern Rift Margin (south), and the
                               non-cluster ("cluster-out") branches of MO/New Madrid
    nshm23_ceus_fault_cluster  the New Madrid cluster ("cluster-in") branches, which nshmp-lib
                               treats as ClusterRuptureSets (see source_info.json for the hazard math)

The source logic trees are walked exactly as nshmp-lib ModelLoader.Fault does (source-tree.json,
source-group.json, rupture-set.json, cluster-set.json, rate-tree.json, mfd-config.json, mfd-map.json),
the MFDs are built as in FaultRuptureSet (rate-tree x SINGLE magnitude branches, moment balanced
epistemic +-0.2 magnitude branches, Gaussian aleatory magnitude variability, GR MFDs with NSHM rupture
floating), and fault surfaces are DefaultGriddedSurfaces with 1 km spacing. Logic trees are collapsed
to mean rates: every rupture's rate is the MFD rate times its MFD branch weight times its source tree
leaf weight, and identical ruptures (same surface, magnitude, and mechanism) are merged.

Each rupture surface is represented by 1 km wide strips between adjacent columns of the gridded
surface (the "segments" of ucla_plha), so that Rjb and Rrup computed by ucla_plha as the minimum over
segments match nshmp-lib's gridded-surface distances to within the grid spacing.

Run from the utilities directory: python convert_nshm23_ceus_fault_source_models.py
Input files are downloaded into utilities/nshm23_ceus/ (git ignored).
'''
import copy
import json
import os
from collections import OrderedDict

import numpy as np

import nshm23_ceus_common as nc

FAULT_ROOT = 'stable-crust/fault'
OUT_FAULT = '../src/ucla_plha/source_models/fault_source_models/nshm23_ceus_fault'
OUT_CLUSTER = '../src/ucla_plha/source_models/fault_source_models/nshm23_ceus_fault_cluster'
RATE_CUTOFF = 1e-14  # FaultRuptureSet skips MFD bins with rate < 1e-14


class Data:
    '''Mutable loader state copied at each directory level (model.ModelData).'''

    def __init__(self):
        self.mfd_config = None
        self.mfd_map = None
        self.rate_tree = None
        self.features = None
        self.cluster_model = False

    def copy(self):
        return copy.copy(self)


def read_features(directory):
    features = {}
    for name in sorted(os.listdir(directory)):
        if name.endswith('.geojson'):
            d = nc.read_json(os.path.join(directory, name))
            for f in (d['features'] if 'features' in d else [d]):
                features[int(f['id'])] = f
    return features


def create_feature(section_ids, features, rupture_props):
    '''FaultRuptureSet.Builder.createFeature(): stitch section traces, use properties of the 1st.'''
    props = dict(features[section_ids[0]]['properties'])
    if rupture_props:
        props.update(rupture_props)
    trace = []
    for sid in section_ids:
        for p in features[sid]['geometry']['coordinates']:
            p = (float(p[0]), float(p[1]))
            if p not in trace:  # LocationList distinct()
                trace.append(p)
    if 'magnitude-note' in props:
        raise NotImplementedError('magnitude-note not supported')
    dip = float(props['dip'])
    upper = float(props['upper-depth'])
    lower = float(props['lower-depth'])
    width = float(props['width']) if 'width' in props else (lower - upper) / np.sin(dip * nc.TO_RAD)
    return {'trace': np.array(trace), 'dip': dip, 'upper_depth': upper, 'width': width,
            'rake': float(props['rake']), 'name': props.get('name')}


def build_mfds(mfd_props_tree, data):
    '''
    FaultRuptureSet.Builder: (optional) rate-tree update, then createMfdTree(). Returns a list of
    (type, mags, rates) with rates already multiplied by the MFD branch weight.
    '''
    cfg = data.mfd_config
    tree = mfd_props_tree
    if data.rate_tree is not None and not data.cluster_model:
        # updateMfdsForRecurrence()
        new_tree = []
        for rid, rvalue, rweight in data.rate_tree:
            rate = 0.0 if rvalue <= 0.0 else rvalue if rvalue <= 1.0 else nc.jround_digits(1.0 / rvalue, 4)
            for mid, props, mweight in tree:
                assert props['type'] == 'SINGLE'
                new_tree.append(('%s:%s' % (mid, rid), {'type': 'SINGLE', 'm': props['m'], 'rate': rate},
                                 mweight * rweight))
        tree = new_tree
    out = []
    for _, props, weight in tree:
        if props['type'] == 'SINGLE':
            mo_rate = props['rate'] * nc.mag_to_moment(props['m'])
            branches = [(props['m'], weight)]
            if cfg.epistemic is not None:
                branches = [(nc.jround(props['m'] + v, 6), weight * w) for _, v, w in cfg.epistemic]
            for m, w in branches:
                if cfg.aleatory is not None:
                    a = cfg.aleatory
                    mags, rates = nc.gaussian_mfd(m, a['size'], a['σ'], a['nσ'])
                else:
                    mags, rates = np.array([m]), np.array([1.0])
                rates = nc.scale_to_moment_rate(mags, rates, mo_rate) if mo_rate > 0 else rates * 0.0
                out.append(('SINGLE', mags, rates * w))
        else:  # GR
            mo_rate = nc.gr_moment_rate(props)
            a, b, dm, mmin, mmax = props['a'], props['b'], props['Δm'], props['mMin'], props['mMax']
            assert nc.mag_count(mmin, mmax, dm) > 0
            epistemic = cfg.epistemic is not None and (mmax - dm / 2 + cfg.min_epi_offset) > cfg.min_magnitude
            if epistemic:
                for _, v, w in cfg.epistemic:
                    n = nc.mag_count(mmin, mmax + v, dm)
                    if n == 0:
                        mags, rates = nc.gr_mfd(1.0, b, dm, mmin, nc.jround(mmin + dm + 0.0001, 6))
                        out.append(('GR', mags, rates * 0.0))
                        continue
                    mags, rates = nc.gr_mfd(a, b, dm, mmin, nc.jround(mmin + n * dm + 0.0001, 6))
                    out.append(('GR', mags, nc.scale_to_moment_rate(mags, rates, mo_rate) * weight * w))
            else:
                mags, rates = nc.gr_mfd(a, b, dm, mmin, mmax)
                out.append(('GR', mags, nc.scale_to_moment_rate(mags, rates, mo_rate) * weight))
    return out


class Model:
    '''Accumulates surfaces and ruptures for one ucla_plha fault source model.'''

    def __init__(self):
        self.surfaces = []           # GriddedSurface objects
        self.surface_index = {}      # key -> index
        self.ruptures = OrderedDict()  # key -> rate

    def surface(self, feature):
        key = (tuple(map(tuple, feature['trace'])), feature['dip'], feature['upper_depth'], feature['width'])
        if key not in self.surface_index:
            self.surface_index[key] = len(self.surfaces)
            self.surfaces.append(nc.GriddedSurface(feature['trace'], feature['upper_depth'],
                                                   feature['dip'], feature['width']))
        return self.surface_index[key]

    def add(self, key, rate):
        self.ruptures[key] = self.ruptures.get(key, 0.0) + rate


def add_rupture_set(model, rs, data, leaf_weight, source_id, cluster=None):
    '''
    Build a FaultRuptureSet and add its ruptures. For cluster members (cluster = (cluster_id,
    section)) the rate stored is nshmp-lib's rupture rate (magnitude branch weight, leaf weight 1).
    Returns the total rate added.
    '''
    sections = rs.get('sections', [rs['id']])
    feature = create_feature(sections, data.features, rs.get('properties'))
    isurf = model.surface(feature)
    surf = model.surfaces[isurf]
    total = 0.0
    for mfd_type, mags, rates in build_mfds(nc.mfd_tree(rs['mfd-tree'], data.mfd_map), data):
        for m, rate in zip(mags, rates):
            if rate < RATE_CUTOFF:
                continue
            if mfd_type == 'GR':  # GR floats, SINGLE does not (no magnitude notes in CEUS)
                floaters = nc.float_nshm(surf, m)
            else:
                floaters = [(0, surf.rows, 0, surf.cols)]
            r = rate * leaf_weight / len(floaters)
            for fl in floaters:
                key = (isurf,) + fl + (round(float(m), 6), feature['rake'], source_id)
                if cluster is not None:
                    key = key + cluster
                model.add(key, r)
                total += r
    return total


def walk(directory, data, weight, ctx):
    '''ModelLoader.Fault.processBranch().'''
    if os.path.basename(directory) == 'null':
        return
    data = data.copy()
    path = os.path.join(directory, 'mfd-config.json')
    if os.path.exists(path):
        data.mfd_config = nc.MfdConfig(path)
    rt = nc.rate_tree(directory)
    if rt is not None:
        data.rate_tree = rt
    tree = nc.read_optional_json(directory, 'source-tree.json')
    if tree is not None:
        for b in tree:
            walk(os.path.join(directory, 'null' if b['id'] is None else b['id']), data, weight * b['weight'], ctx)
        return
    group = nc.read_optional_json(directory, 'source-group.json')
    if group is not None:
        for b in group:
            walk(os.path.join(directory, b['id']), data, weight * b.get('scale', 1.0), ctx)
        return
    rs = nc.read_optional_json(directory, 'rupture-set.json')
    rel = os.path.relpath(directory, ctx['root']).replace(os.sep, '/')
    if rs is not None:
        leaf = nc.jround(weight, 8)
        total = add_rupture_set(ctx['fault'], rs, data, leaf, ctx['source_id'])
        ctx['fault_leaves'].append({'path': rel, 'name': rs['name'], 'leaf_weight': leaf, 'total_rate': total})
        return
    cs = nc.read_json(os.path.join(directory, 'cluster-set.json'))
    data.cluster_model = True
    cluster_id = len(ctx['clusters'])
    sections = []
    for isec, member in enumerate(cs['rupture-sets']):
        total = add_rupture_set(ctx['cluster'], member, data, 1.0, ctx['source_id'], (cluster_id, isec))
        sections.append({'name': member['name'], 'id': member['id'], 'sum_of_rates': total})
    # ModelLoader.processClusterBranch(): one ClusterRuptureSet per rate-tree branch, all with the
    # same geometry. Mean hazard sum_b w_b * rate_b * P(x) = (sum_b w_b * rate_b) * P(x).
    branches = [(rid, value, nc.jround(weight * w, 8)) for rid, value, w in data.rate_tree]
    effective = sum(value * w for _, value, w in branches)
    ctx['clusters'].append({
        'cluster_id': cluster_id, 'path': rel, 'name': cs['name'], 'nshm_id': cs['id'],
        'weight': nc.jround(weight, 8), 'effective_rate': effective,
        'rate_branches': [{'id': rid, 'rate': value, 'leaf_weight': w} for rid, value, w in branches],
        'sections': sections, 'source_id': ctx['source_id']})


def load_source_tree(directory, data, ctx):
    '''ModelLoader.Fault.loadSourceTree() for a tree of rupture sets with a features directory.'''
    info = nc.read_json(os.path.join(directory, 'tree-info.json'))
    data = data.copy()
    data.features = read_features(os.path.join(directory, 'features'))
    ctx['root'] = directory
    ctx['source_id'] = info['id']
    ctx['tree_names'][info['id']] = os.path.relpath(directory, ctx['fault_root']).replace(os.sep, '/')
    for b in nc.read_json(os.path.join(directory, 'source-tree.json')):
        walk(os.path.join(directory, 'null' if b['id'] is None else b['id']), data, b['weight'], ctx)


def load_directory(directory, data, ctx):
    '''ModelLoader.loadSourceDirectory() for the fault source type.'''
    data = data.copy()
    path = os.path.join(directory, 'mfd-config.json')
    if os.path.exists(path):
        data.mfd_config = nc.MfdConfig(path)
    path = os.path.join(directory, 'mfd-map.json')
    if os.path.exists(path):
        data.mfd_map = nc.read_json(path)
    assert not os.path.exists(os.path.join(directory, 'rate-tree.json'))
    if os.path.exists(os.path.join(directory, 'source-tree.json')):
        load_source_tree(directory, data, ctx)
        return
    assert not any(n.endswith('.geojson') for n in os.listdir(directory)), directory
    for name in sorted(os.listdir(directory)):
        sub = os.path.join(directory, name)
        if os.path.isdir(sub) and name != 'features':
            load_directory(sub, data, ctx)


def write_model(model, output_dir, cluster):
    '''Write ruptures.npz, ruptures_segments.npz, and the tri_* and rect_* geometry arrays.'''
    os.makedirs(output_dir, exist_ok=True)
    keys = list(model.ruptures.keys())
    # Order ruptures by surface and position so that segments are grouped.
    order = sorted(range(len(keys)), key=lambda i: keys[i][:8] + (keys[i][8:] if cluster else ()))
    keys = [keys[i] for i in order]
    segment_ids = {}
    corners = []  # (p1, p2, p3, p4) each (lat, lon, depth)
    rupture_index = []
    segment_index = []
    for irup, key in enumerate(keys):
        isurf, r0, nr, c0, ncol = key[:5]
        grid = model.surfaces[isurf].grid
        r1 = r0 + nr - 1
        for c in range(c0, c0 + ncol - 1):
            skey = (isurf, r0, r1, c)
            if skey not in segment_ids:
                segment_ids[skey] = len(corners)
                p = [grid[r0, c], grid[r0, c + 1], grid[r1, c], grid[r1, c + 1]]
                corners.append([(q[1], q[0], q[2]) for q in p])
            rupture_index.append(irup)
            segment_index.append(segment_ids[skey])
    corners = np.array(corners)  # n x 4 x 3 (lat, lon, depth)
    lat, lon, dep = corners[..., 0], corners[..., 1], corners[..., 2]
    zero = np.zeros_like(dep)
    xyz = nc.lat_lon_depth_to_xyz(lat, lon, dep)
    xyz0 = nc.lat_lon_depth_to_xyz(lat, lon, zero)
    # Triangles (p1, p2, p4) and (p1, p3, p4), as in the nshm23_wus and UCERF3 conversions
    tri_rrup = np.concatenate((xyz[:, [0, 1, 3]], xyz[:, [0, 2, 3]]))
    tri_rjb = np.concatenate((xyz0[:, [0, 1, 3]], xyz0[:, [0, 2, 3]]))
    seg = np.arange(len(corners), dtype='intc')
    np.save(os.path.join(output_dir, 'tri_segment_id.npy'), np.concatenate((seg, seg)))
    np.save(os.path.join(output_dir, 'tri_rrup.npy'), tri_rrup)
    np.save(os.path.join(output_dir, 'tri_rjb.npy'), tri_rjb)
    np.save(os.path.join(output_dir, 'rect_segment_id.npy'), seg)
    np.save(os.path.join(output_dir, 'rect_rjb.npy'), xyz0)
    np.savez_compressed(os.path.join(output_dir, 'ruptures_segments.npz'),
                        rupture_index=np.asarray(rupture_index, dtype=np.int32),
                        segment_index=np.asarray(segment_index, dtype=np.int32))

    rate = np.array([model.ruptures[k] for k in keys])
    m = np.array([k[5] for k in keys])
    rake = np.array([k[6] for k in keys])
    source_id = np.array([k[7] for k in keys], dtype=np.int32)
    dip = np.array([model.surfaces[k[0]].dip for k in keys])
    ztor = np.array([model.surfaces[k[0]].grid[k[1], 0, 2] for k in keys])
    zbor = np.array([model.surfaces[k[0]].grid[k[1] + k[2] - 1, 0, 2] for k in keys])
    width = np.array([model.surfaces[k[0]].dip_spacing * (k[2] - 1) for k in keys])
    zhyp = ztor + 0.5 * width * np.sin(np.radians(dip))  # Faults.hypocentralDepth
    arrays = dict(m=m, rate=rate, fault_type=nc.fault_type_from_rake(rake), dip=dip, ztor=ztor, zbor=zbor,
                  rake=rake, width=width, zhyp=zhyp, source_id=source_id)
    if cluster:
        arrays['cluster_id'] = np.array([k[8] for k in keys], dtype=np.int32)
        # cluster_section is unique over all clusters (same convention as nshm23_cascadia_interface_cluster)
        sections = {s: i for i, s in enumerate(sorted({(k[8], k[9]) for k in keys}))}
        arrays['cluster_section'] = np.array([sections[(k[8], k[9])] for k in keys], dtype=np.int32)
    np.savez_compressed(os.path.join(output_dir, 'ruptures.npz'), **arrays)
    write_grid(model, keys, output_dir)
    return arrays


def write_grid(model, keys, output_dir):
    '''
    Write grid.npz: the nshmp-lib grid points of every surface and the grid subset of every rupture
    (in ruptures.npz order), used by ucla_plha (source_info "fault_distances": "nshmp_grid",
    geometry.gridded_rupture_distances) to compute rJB, rRup, and rX exactly as nshmp-lib
    Distance.compute does. Points are stored surface by surface and column by column (point (r, c) of
    surface s is surface_offset[s] + c * rows + r), as float32 longitude, latitude, and depth
    (rounding error < 1 m).
    '''
    surfaces = model.surfaces
    points = np.concatenate([s.grid.transpose(1, 0, 2).reshape(-1, 3) for s in surfaces])
    sizes = [s.rows * s.cols for s in surfaces]
    np.savez_compressed(
        os.path.join(output_dir, 'grid.npz'),
        lon=points[:, 0].astype(np.float32), lat=points[:, 1].astype(np.float32),
        depth=points[:, 2].astype(np.float32),
        surface_offset=np.r_[0, np.cumsum(sizes)[:-1]].astype(np.int64),
        surface_rows=np.array([s.rows for s in surfaces], dtype=np.int32),
        surface_cols=np.array([s.cols for s in surfaces], dtype=np.int32),
        surface_dip=np.array([s.dip for s in surfaces]),
        surface_spacing=np.array([(s.strike_spacing + s.dip_spacing) / 2 for s in surfaces]),
        rupture_surface=np.array([k[0] for k in keys], dtype=np.int32),
        rupture_row0=np.array([k[1] for k in keys], dtype=np.int32),
        rupture_rows=np.array([k[2] for k in keys], dtype=np.int32),
        rupture_col0=np.array([k[3] for k in keys], dtype=np.int32),
        rupture_cols=np.array([k[4] for k in keys], dtype=np.int32))


STABLE_GMM_TREE = ('stable-crust/gmm-tree.json: NGA_EAST_2026 0.3333, NGA_EAST_2026_ADJUSTED 0.3333, '
                   'NGA_EAST_SEEDS_2026 0.1667, NGA_EAST_SEEDS_2026_ADJUSTED 0.1667; stable-crust/gmm-config.json '
                   'max-distance 1000 km')
RATE_SERVICE = 'USGS NSHM rate web service (https://earthquake.usgs.gov/ws/nshmp/conus-2023/dynamic/rate/lon/lat/r)'
GEOMETRY = ('Rupture surfaces are nshmp-lib DefaultGriddedSurfaces (fault-config.json surface-spacing 1 km): the trace '
            '(sections stitched in rupture-set order, properties of the first section, rupture-set "properties" '
            'such as width override) is resampled at length/ceil(length) km, and each column is projected down-dip '
            'by the down-dip width (feature "width", else (lower - upper depth) / sin(dip)) in one direction, strike '
            'of the first-to-last trace point + 90 degrees. Each ucla_plha segment is the quadrilateral between two '
            'adjacent grid columns over the rows of a rupture, so Rjb, Rrup, Rx, Ry0 = min over a rupture\'s '
            'segments. nshmp-lib (model.Distance.compute) uses the minimum distance to the 1 km grid points '
            '(horizontal distance Locations.horzDistanceFast with R = 6371.0072 km, rRup = sqrt(horizontal^2 + '
            'depth^2), top row only for dip > 89), rJB = 0 if less than the grid spacing and inside the surface '
            'perimeter, and rX from the upper edge extended 1000 km along strike. ucla_plha computes these from '
            'grid.npz (the grid points and the grid subset of each rupture; source_info "fault_distances": '
            '"nshmp_grid", geometry.gridded_rupture_distances). With "constraints": {"fault_distances": '
            '"triangles"} it uses the segments instead, whose Cartesian distances differ from the nshmp-lib '
            'flat-earth ones by < 0.1 km within 50 km but by 2 to 3 km at 900 km (0.9% higher New Madrid cluster '
            'hazard at Charleston at 2475 years). GMM inputs as in nshmp-lib '
            '(calc.Transforms.IterableToInputs): ztor = depth of the rupture top row, dip, width (ruptures.npz '
            '"width" = down-dip width of the rupture), zhyp = ztor + width sin(dip) / 2, rake.')


def write_source_info(ctx, fault, cluster):
    names = {sid: name for sid, name in ctx['tree_names'].items()}
    fault_sources = {}
    for sid in sorted(set(fault['source_id'].tolist())):
        s = fault['source_id'] == sid
        fault_sources[str(sid)] = {'path': 'stable-crust/fault/' + names[sid], 'n_ruptures': int(s.sum()),
                                   'rate': float(fault['rate'][s].sum()),
                                   'rate_m_ge_7': float(fault['rate'][s & (fault['m'] >= 7.0)].sum())}
    info = {
        'name': 'nshm23_ceus_fault',
        'tectonic_region': 'stable_crust',
        'nshm_component': 'Fault',
        'source': 'nshm-conus %s stable-crust/fault: CO/Cheraw, MO/Commerce, OK/Meers, TN/Eastern Rift Margin '
                  '(north), TN/Eastern Rift Margin (south), and the non-cluster branches of MO/New Madrid '
                  '(sscn/cluster-out, usgs/*/cluster-out); fault-config.json, mfd-config.json files' % nc.NSHM_REF,
        'nshm_release': nc.NSHM_REF,
        'nshmp_lib_commit': nc.NSHMP_LIB_COMMIT,
        'cluster': False,
        'fault_distances': 'nshmp_grid',
        'gmm_tree': STABLE_GMM_TREE,
        'geometry': GEOMETRY,
        'files': {
            'ruptures.npz': 'm, rate (annual, mean over the logic tree), fault_type (1 reverse, 2 normal, 3 strike '
                            'slip, from rake), dip, ztor, zbor, plus rake, width (down-dip, km), zhyp, source_id '
                            '(NSHM source tree id, see "sources")',
            'ruptures_segments.npz, tri_*.npy, rect_*.npy': 'as for nshm23_wus',
        },
        'sources': fault_sources,
        'leaves': ctx['fault_leaves'],
        'notes': 'Built by utilities/convert_nshm23_ceus_fault_source_models.py, which follows nshmp-lib '
                 'ModelLoader.Fault and FaultRuptureSet: rate-tree x SINGLE magnitude branches (updateMfdsForRecurrence), '
                 'mfd-config.json epistemic +-0.2 magnitude branches (weights 0.2/0.6/0.2, moment balanced to the '
                 'SINGLE moment rate; stable-crust/fault/mfd-config.json, used by Cheraw (USGS) and Meers (USGS)) and '
                 'Gaussian aleatory variability (11 magnitudes over +-2 sigma, sigma 0.12, moment balanced; Cheraw, '
                 'Meers), GR MFDs (Cheraw partial ruptures) with epistemic mMax branches and NSHM rupture floating '
                 '(WC94 length; 1 to 4 down-dip positions, rate split evenly over floaters). SINGLE MFDs do not '
                 'float (full-surface ruptures). MO/ and TN/ mfd-config.json turn epistemic and aleatory terms off '
                 '(Commerce, Eastern Rift Margin, New Madrid). Rupture rate = MFD rate x MFD branch weight x leaf '
                 'weight (product of source-tree weights, rounded to 8 decimals as in SourceTree); identical '
                 'ruptures are merged. Validated against the %s Fault component: within 60 km of Cheraw, 50 km of '
                 'Meers, and 300 km of New Madrid (Commerce, Eastern Rift Margin north and south, New Madrid '
                 'cluster-out) the rates in every 0.1 magnitude bin agree to < 1e-6 (relative).' % RATE_SERVICE,
    }
    with open(os.path.join(OUT_FAULT, 'source_info.json'), 'w', encoding='utf-8') as f:
        json.dump(info, f, indent=2, ensure_ascii=False)

    clusters = ctx['clusters']
    cinfo = {
        'name': 'nshm23_ceus_fault_cluster',
        'tectonic_region': 'stable_crust',
        'nshm_component': 'FaultCluster',
        'source': 'nshm-conus %s stable-crust/fault/MO/New Madrid: sscn/cluster-in (8 cluster sets) and '
                  'usgs/{west,mid-west,center,mid-east,east}/cluster-in/{all,center-south} (10 cluster sets)'
                  % nc.NSHM_REF,
        'nshm_release': nc.NSHM_REF,
        'nshmp_lib_commit': nc.NSHMP_LIB_COMMIT,
        'cluster': True,
        'fault_distances': 'nshmp_grid',
        'gmm_tree': STABLE_GMM_TREE,
        'geometry': GEOMETRY,
        'files': {
            'ruptures.npz': 'm, rate (= magnitude branch WEIGHT of the rupture within its cluster section, not an '
                            'annual rate; equal magnitudes in a section are merged by adding weights; the weights of '
                            'a section sum to 1), fault_type, dip, ztor, zbor, rake, width, zhyp, source_id, '
                            'cluster_id (rupture -> cluster), cluster_section (rupture -> fault section of its '
                            'cluster; unique over all clusters, as in nshm23_cascadia_interface_cluster)',
            'clusters.npz': 'cluster_id, rate (annual cluster rate), weight (logic tree weight), n_sections, '
                            'source_id; the cluster contributes weight x rate x P_c(x)',
            'clusters.json': 'per cluster: source tree path, NSHM cluster-set name and id, weight, rate-tree '
                             'branches (rate, leaf weight), sections',
        },
        'cluster_hazard': 'nshmp-lib ClusterRuptureSet (source type FAULT_CLUSTER, FaultCluster component) and '
                          'calc.Transforms.ClusterGroundMotionsToCurves. When a cluster occurs, every one of its '
                          'sections (2 or 3 New Madrid faults) ruptures, each with one magnitude chosen from its '
                          'magnitude branches. For each GMM g (and, in nshmp-lib, each epistemic branch of the GMM, '
                          'combined after the cluster calculation): P_s(x) = sum_{r in section s} rate_r * '
                          'P_g(Y > x | r) (exceedance with the hazard exceedance model, truncated at +3 sigma); '
                          'P_c(x) = 1 - prod_s (1 - P_s(x)); annual rate of exceedance lambda(x) = sum_c weight_c * '
                          'rate_c * sum_g w_g * P_c,g(x). The product over sections is per GMM. Ruptures are not '
                          'independent Poisson events. A cluster is included if any section trace point is within '
                          'the GMM maximum distance (1000 km); individual ruptures are not distance filtered.',
        'cluster_weighting': 'Difference from nshm23_cascadia_interface_cluster: nshmp-lib creates one '
                             'ClusterRuptureSet per rate-tree branch of a cluster set (same geometry, rate rate_b, '
                             'leaf weight w_path x w_b). Because lambda(x) is linear in weight x rate for fixed '
                             'geometry, the rate branches are collapsed: clusters.npz weight = source-tree path '
                             'weight w_path of the cluster set, rate = sum_b round8(w_path w_b) rate_b / w_path '
                             '(the rate-tree mean). weight x rate is exact; clusters.json lists the branches. '
                             'USGS "all" and "center-south" cluster sets are a source-group (both apply, scale 1.0). '
                             'Null rate-tree values are 0.',
        'clusters': [{'cluster_id': c['cluster_id'], 'name': c['name'], 'path': c['path'],
                      'weight': c['weight'], 'rate': c['effective_rate'] / c['weight'],
                      'n_sections': len(c['sections'])} for c in clusters],
        'notes': 'Section rupture sets are FaultRuptureSets built with the cluster flag (no rate-tree update) from '
                 'the cluster-set "rupture-sets" (SINGLE magnitude trees in New Madrid/mfd-map.json; no epistemic or '
                 'aleatory terms, MO/mfd-config.json). Total cluster rate sum_c weight_c rate_c = %.6e/yr. Validated '
                 'against the %s FaultCluster component within 300 km of New Madrid (rate service: sum_c weight_c '
                 'rate_c sum_r rate_r per magnitude bin): all 0.1 magnitude bins agree to < 1e-6.' % (
                     sum(c['effective_rate'] for c in clusters), RATE_SERVICE),
    }
    with open(os.path.join(OUT_CLUSTER, 'source_info.json'), 'w', encoding='utf-8') as f:
        json.dump(cinfo, f, indent=2, ensure_ascii=False)


def main():
    root = nc.fetch_tree(FAULT_ROOT)
    data = Data()
    ctx = {'fault': Model(), 'cluster': Model(), 'clusters': [], 'fault_leaves': [], 'tree_names': {},
           'fault_root': root}
    load_directory(root, data, ctx)

    fault = write_model(ctx['fault'], OUT_FAULT, cluster=False)
    cluster = write_model(ctx['cluster'], OUT_CLUSTER, cluster=True)
    clusters = ctx['clusters']
    np.savez_compressed(os.path.join(OUT_CLUSTER, 'clusters.npz'),
                        cluster_id=np.array([c['cluster_id'] for c in clusters], dtype=np.int32),
                        rate=np.array([c['effective_rate'] / c['weight'] for c in clusters]),
                        weight=np.array([c['weight'] for c in clusters]),
                        n_sections=np.array([len(c['sections']) for c in clusters], dtype=np.int32),
                        source_id=np.array([c['source_id'] for c in clusters], dtype=np.int32))
    with open(os.path.join(OUT_CLUSTER, 'clusters.json'), 'w') as f:
        json.dump(clusters, f, indent=1)

    write_source_info(ctx, fault, cluster)

    used = sorted(set(fault['source_id'].tolist()))
    print('nshm23_ceus_fault: %d ruptures, %d sources' % (len(fault['m']), len(used)))
    for sid in used:
        s = fault['source_id'] == sid
        print('  %-40s %6d  rate %.4e  rate(M>=7) %.4e' % (ctx['tree_names'][sid], s.sum(), fault['rate'][s].sum(),
                                                          fault['rate'][s & (fault['m'] >= 7)].sum()))
    print('nshm23_ceus_fault_cluster: %d ruptures, %d clusters' % (len(cluster['m']), len(clusters)))
    return ctx


if __name__ == '__main__':
    main()

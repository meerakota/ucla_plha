'''
Shared helpers for converting the USGS NSHM conterminous U.S. model (nshm-conus release 6.2.0)
stable-crust sources into ucla_plha source models. Used by convert_nshm23_ceus_fault_source_models.py
and convert_nshm23_ceus_point_source_models.py.

The functions below re-implement, in Python, the parts of USGS nshmp-lib (commit 44728a7d) that turn
the model's JSON, CSV, and GeoJSON files into ruptures with annual rates:

    gov.usgs.earthquake.nshmp.mfd.Mfd                      MFD builders (SINGLE, Gaussian, GR)
    gov.usgs.earthquake.nshmp.data.Sequences               magnitude sequences
    gov.usgs.earthquake.nshmp.geo.Locations                spherical geodesy (R = 6371.0072 km)
    gov.usgs.earthquake.nshmp.fault.surface.DefaultGriddedSurface    fault surfaces
    gov.usgs.earthquake.nshmp.fault.surface.RuptureFloating.NSHM     floating ruptures
    gov.usgs.earthquake.nshmp.model.FaultRuptureSet         fault MFD logic trees

Model input files are downloaded with curl (the code.usgs.gov API refuses Python's urllib) into
utilities/nshm23_ceus/, mirroring the paths in the nshm-conus repository.
'''
import json
import os
import subprocess
import urllib.parse
from decimal import ROUND_DOWN, ROUND_HALF_UP, Decimal

import numpy as np

NSHM_REF = '6.2.0'
NSHM_PROJECT = 6063
NSHMP_LIB_COMMIT = '44728a7d'
API = 'https://code.usgs.gov/api/v4/projects/%d/repository' % NSHM_PROJECT
INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nshm23_ceus')

EARTH_RADIUS_MEAN = 6371.0072  # gov.usgs.earthquake.nshmp.geo.Coordinates
TO_RAD = np.pi / 180.0


# ------------------------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------------------------

def _curl(url, output=None):
    cmd = ['curl', '-s', '-f', '-L', url]
    if output is not None:
        cmd += ['-o', output]
    result = subprocess.run(cmd, capture_output=True, check=True)
    return result.stdout


def list_tree(path):
    '''Return the paths of all files below path in the nshm-conus repository.'''
    paths = []
    page = 1
    while True:
        url = '%s/tree?path=%s&ref=%s&per_page=100&recursive=true&page=%d' % (
            API, urllib.parse.quote(path), NSHM_REF, page)
        entries = json.loads(_curl(url))
        if not entries:
            break
        paths += [e['path'] for e in entries if e['type'] == 'blob']
        page += 1
    return paths


def fetch(path):
    '''Download one nshm-conus file into INPUT_DIR (if not already there) and return its local path.'''
    local = os.path.join(INPUT_DIR, *path.split('/'))
    if not os.path.exists(local) or os.path.getsize(local) == 0:
        os.makedirs(os.path.dirname(local), exist_ok=True)
        url = '%s/files/%s/raw?ref=%s' % (API, urllib.parse.quote(path, safe=''), NSHM_REF)
        _curl(url, local)
    return local


def fetch_tree(path):
    '''Download all files below path. Returns the local directory.'''
    marker = os.path.join(INPUT_DIR, *path.split('/'), '.complete')
    if not os.path.exists(marker):
        for p in list_tree(path):
            fetch(p)
        open(marker, 'w').close()
    return os.path.join(INPUT_DIR, *path.split('/'))


def read_json(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def read_optional_json(directory, name):
    path = os.path.join(directory, name)
    return read_json(path) if os.path.exists(path) else None


# ------------------------------------------------------------------------------------------------
# Java numerics
# ------------------------------------------------------------------------------------------------

def jround(value, scale, mode=ROUND_HALF_UP):
    '''nshmp Maths.round: BigDecimal.valueOf(value).setScale(scale, mode).'''
    return float(Decimal(repr(float(value))).quantize(Decimal(1).scaleb(-scale), rounding=mode))


def jround_down(value, scale):
    return jround(value, scale, ROUND_DOWN)


def jround_digits(value, digits):
    '''nshmp Maths.roundToDigits: BigDecimal.round(new MathContext(digits, HALF_UP)).'''
    d = Decimal(repr(float(value)))
    if d == 0:
        return 0.0
    exponent = d.adjusted() - digits + 1
    return float(d.quantize(Decimal(1).scaleb(exponent), rounding=ROUND_HALF_UP))


def jfloat_round(x):
    '''Java Math.round((float) x).'''
    return int(np.floor(np.float32(x) + np.float32(0.5)))


def jrint(x):
    '''Java Math.rint: round half to even.'''
    return int(np.rint(x))


def mag_to_moment(m):
    '''Earthquakes.magToMoment (N-m).'''
    return 10.0 ** (1.5 * np.asarray(m, dtype=float) + 9.05)


# ------------------------------------------------------------------------------------------------
# Magnitude-frequency distributions (gov.usgs.earthquake.nshmp.mfd.Mfd)
# ------------------------------------------------------------------------------------------------

def sequence(mmin, mmax, delta, centered, scale=5):
    '''Sequences.arrayBuilder(min, max, Δ)[.centered()].scale(scale).build().'''
    lo = mmin if centered else mmin + delta / 2.0
    hi = mmax if centered else mmax - delta / 2.0
    size = int(jround((hi - lo) / delta, 6)) + 1
    return np.array([jround(lo + delta * i, scale) for i in range(size)])


def gr_mfd(a, b, dm, mmin, mmax):
    '''GutenbergRichter.toBuilder(): bin centers between bin edges mmin and mmax, rates 10^(a - b m).'''
    m = sequence(mmin, mmax, dm, centered=False)
    return m, 10.0 ** (a - b * m)


def gaussian_mfd(m, size, sigma, nsigma):
    '''Single.toGaussianBuilder(): normal pdf between ±nσ, scaled to a cumulative rate of 1.'''
    lo = jround(m - nsigma * sigma, 6)
    hi = jround(m + nsigma * sigma, 6)
    dm = jround((hi - lo) / (size - 1), 6)
    mags = sequence(lo, hi, dm, centered=True)
    pdf = np.exp((m - mags) * (mags - m) / (2 * sigma * sigma)) / (sigma * np.sqrt(2 * np.pi))
    return mags, pdf / pdf.sum()


def scale_to_moment_rate(mags, rates, moment_rate):
    total = np.sum(rates * mag_to_moment(mags))
    return rates * (moment_rate / total)


def gr_moment_rate(gr):
    '''FaultRuptureSet.grMomentRate().'''
    m, r = gr_mfd(gr['a'], gr['b'], gr['Δm'], gr['mMin'], gr['mMax'] + 0.001)
    incr = 10.0 ** (gr['a'] - gr['b'] * (gr['mMin'] + gr['Δm'] / 2))
    r = r * incr / r[0]
    return float(np.sum(r * mag_to_moment(m)))


def mag_count(mmin, mmax, dm):
    '''FaultRuptureSet.magCount().'''
    return int((mmax - mmin - dm) / dm + 1.4)


def mfd_properties(value):
    '''Deserialize.mfdProperties(): Gson defaults are rate = 1.0 for SINGLE and a = 1.0 for GR.'''
    props = dict(value)
    if props['type'] == 'SINGLE':
        props.setdefault('rate', 1.0)
    elif props['type'] == 'GR':
        props.setdefault('a', 1.0)
    else:
        raise ValueError('Unsupported MFD type %s' % props['type'])
    return props


def mfd_tree(element, mfd_map):
    '''Deserialize.mfdTree(): an inline tree or a key into the mfd-map.'''
    if isinstance(element, str):
        element = mfd_map[element]
    return [(b['id'], mfd_properties(b['value']), b['weight']) for b in element]


def rate_tree(directory):
    '''Deserialize.rateTree(): null values become 0.'''
    tree = read_optional_json(directory, 'rate-tree.json')
    if tree is None:
        return None
    return [(b['id'], 0.0 if b['value'] is None else float(b['value']), b['weight']) for b in tree]


class MfdConfig:
    '''model.MfdConfig read from mfd-config.json.'''

    def __init__(self, path):
        d = read_json(path)
        self.epistemic = None
        if d['epistemic-tree'] is not None:
            self.epistemic = [(b['id'], b['value'], b['weight']) for b in d['epistemic-tree']]
        self.min_epi_offset = None if self.epistemic is None else min(v for _, v, _ in self.epistemic)
        self.aleatory = d['aleatory-properties']
        self.min_magnitude = d['minimum-magnitude']
        self.nshm_bin_model = True if d['nshm-bin-model'] is None else d['nshm-bin-model']

    def describe(self):
        return {'epistemic-tree': self.epistemic, 'aleatory-properties': self.aleatory,
                'minimum-magnitude': self.min_magnitude}


# ------------------------------------------------------------------------------------------------
# Spherical geodesy (gov.usgs.earthquake.nshmp.geo.Locations), angles in radians
# ------------------------------------------------------------------------------------------------

def azimuth_rad(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = (np.radians(v) for v in (lon1, lat1, lon2, lat2))
    dlon = lon2 - lon1
    az = np.arctan2(np.sin(dlon) * np.cos(lat2),
                    np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon))
    return (az + 2 * np.pi) % (2 * np.pi)


def horz_distance(lon1, lat1, lon2, lat2):
    '''Haversine distance (km).'''
    lon1, lat1, lon2, lat2 = (np.radians(v) for v in (lon1, lat1, lon2, lat2))
    s1 = np.sin((lat2 - lat1) / 2.0)
    s2 = np.sin((lon2 - lon1) / 2.0)
    c = s1 * s1 + np.cos(lat1) * np.cos(lat2) * s2 * s2
    return EARTH_RADIUS_MEAN * 2.0 * np.arctan2(np.sqrt(c), np.sqrt(1 - c))


def horz_distance_fast(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = (np.radians(v) for v in (lon1, lat1, lon2, lat2))
    dlat = lat1 - lat2
    dlon = (lon1 - lon2) * np.cos((lat1 + lat2) * 0.5)
    return EARTH_RADIUS_MEAN * np.sqrt(dlat * dlat + dlon * dlon)


def location(lon, lat, az, dh):
    '''Locations.location(p, LocationVector(az, Δh, Δv)); returns (lon, lat).'''
    lat1 = np.radians(lat)
    lon1 = np.radians(lon)
    ad = dh / EARTH_RADIUS_MEAN
    lat2 = np.arcsin(np.sin(lat1) * np.cos(ad) + np.cos(lat1) * np.sin(ad) * np.cos(az))
    lon2 = lon1 + np.arctan2(np.sin(az) * np.sin(ad) * np.cos(lat1), np.cos(ad) - np.sin(lat1) * np.sin(lat2))
    return np.degrees(lon2), np.degrees(lat2)


# ------------------------------------------------------------------------------------------------
# Fault surfaces (DefaultGriddedSurface) and floating ruptures (RuptureFloating.NSHM)
# ------------------------------------------------------------------------------------------------

class GriddedSurface:
    '''
    DefaultGriddedSurface built from a trace, upper depth, dip, and down-dip width with 1 km nominal
    spacing. The trace is resampled evenly along strike, and every column is projected down-dip in a
    single direction normal to the line joining the trace end points (right-hand rule).
    grid[row, col] = (lon, lat, depth).
    '''

    def __init__(self, trace, depth, dip, width, spacing=1.0):
        trace = np.asarray(trace, dtype=float)
        self.trace = trace
        self.dip = float(dip)
        self.depth = float(depth)
        dip_rad = self.dip * TO_RAD
        dip_dir = (azimuth_rad(trace[0, 0], trace[0, 1], trace[-1, 0], trace[-1, 1]) + np.pi / 2) % (2 * np.pi)
        self.dip_dir = dip_dir
        length = np.sum(horz_distance_fast(trace[:-1, 0], trace[:-1, 1], trace[1:, 0], trace[1:, 1]))
        self.strike_spacing = length / np.ceil(length / spacing)
        self.dip_spacing = width / np.ceil(width / spacing)
        seg_len = horz_distance(trace[:-1, 0], trace[:-1, 1], trace[1:, 0], trace[1:, 1])
        seg_az = azimuth_rad(trace[:-1, 0], trace[:-1, 1], trace[1:, 0], trace[1:, 1])
        cum = np.cumsum(seg_len)
        nseg = len(seg_len)
        rows = 1 + jfloat_round(width / self.dip_spacing)
        cols = 1 + jfloat_round(cum[-1] / self.strike_spacing)
        grid = np.empty((rows, cols, 3))
        for c in range(cols):
            along = c * self.strike_spacing
            s = 1
            while s <= nseg and along > cum[s - 1]:
                s += 1
            if s == nseg + 1:
                s -= 1
            d = along - cum[s - 2] if s > 1 else along
            lon, lat = location(trace[s - 1, 0], trace[s - 1, 1], seg_az[s - 1], d)
            grid[0, c] = (lon, lat, self.depth)
            for r in range(1, rows):
                h = r * self.dip_spacing * np.cos(dip_rad)
                v = r * self.dip_spacing * np.sin(dip_rad)
                lon2, lat2 = location(lon, lat, dip_dir, h)
                grid[r, c] = (lon2, lat2, self.depth + v)
        self.grid = grid
        self.rows = rows
        self.cols = cols
        self.width = self.dip_spacing * (rows - 1)


def wc94_length(m):
    '''RuptureScaling.lengthWc94.'''
    return 10.0 ** (-3.22 + 0.69 * m)


def float_nshm(surface, m):
    '''
    RuptureFloating.NSHM with RuptureScaling.NSHM_FAULT_WC94_LENGTH. Returns a list of
    (start_row, n_rows, start_col, n_cols) for the floating ruptures (each gets rate / count).
    '''
    ztor = surface.depth
    down_dip = 1 if (ztor > 1.0 or m > 7.0) else 2 if m > 6.75 else 3 if m > 6.5 else 4
    dz = 2.0 / np.sin(surface.dip * TO_RAD)
    out = []
    for i in range(down_dip):
        zw = i * dz
        length = wc94_length(m)
        width = min(surface.width - zw, length)
        start_row = jrint(zw / surface.dip_spacing)
        n_rows = jrint(width / surface.dip_spacing + 1)
        n_cols = jrint(length / surface.strike_spacing + 1)
        along = surface.cols - n_cols + 1
        if along <= 1:
            along = 1
            n_cols = surface.cols
        for c in range(along):
            out.append((start_row, n_rows, c, n_cols))
    return out


# ------------------------------------------------------------------------------------------------
# Cartesian conversion used by ucla_plha (same as geometry.point_to_xyz, depth positive down)
# ------------------------------------------------------------------------------------------------

def earth_radius(lat):
    rad = np.pi / 180.0
    a = 6378.1370
    b = 6356.7523
    return np.sqrt(((a**2 * np.cos(lat * rad))**2 + (b**2 * np.sin(lat * rad))**2) /
                   ((a * np.cos(lat * rad))**2 + (b * np.sin(lat * rad))**2))


def lat_lon_depth_to_xyz(lat, lon, depth):
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    depth = np.asarray(depth, dtype=float)
    rad = np.pi / 180.0
    r = earth_radius(lat) - depth
    return np.stack((r * np.cos(lat * rad) * np.cos(lon * rad),
                     r * np.cos(lat * rad) * np.sin(lon * rad),
                     r * np.sin(lat * rad)), axis=-1)


def fault_type_from_rake(rake):
    '''ucla_plha fault_type: 1 = reverse, 2 = normal, 3 = strike slip (same mapping as nshm23_wus).'''
    rake = np.asarray(rake, dtype=float)
    ft = np.full(rake.shape, 1)
    ft[(rake > -150) & (rake < -30)] = 2
    ft[(rake >= -180) & (rake <= -150)] = 3
    ft[(rake >= -30) & (rake <= 30)] = 3
    ft[(rake >= 150) & (rake <= 180)] = 3
    return ft

"""Point source distances following the USGS nshmp-lib conventions.

The functions in this module reproduce how nshmp-lib (commit 44728a7d,
``gov.usgs.earthquake.nshmp.model.GridSource``, ``GridSourceFinite``, and
``gov.usgs.earthquake.nshmp.fault.surface.RuptureScaling``) turns a point source (grid node
or zone node) into distances for the ground motion models:

* The horizontal distance from the site to the node is computed with nshmp-lib's
  ``Locations.horzDistanceFast`` (flat-Earth approximation with a mean Earth radius of
  6371.0072 km).
* ``RuptureScaling.pointSourceDistance``: for M >= 6, the horizontal distance is replaced by
  the average Joyner-Boore distance of finite ruptures of random strike centered on the
  node, read from the nshmp-lib tables (``nshmp_rjb_tables.npz``, magnitudes 6.05 to 8.55,
  distances 0 to 1000 km; magnitude index ``round((M - 6.05) / 0.1)`` capped at the last
  row, distance index ``floor(r)`` capped at 1000). The NONE scaling uses the horizontal
  distance, and PEER uses ``max(0.5, r)``.
* ``POINT`` sources: rRup = hypot(rJB, zTor) and rX = rJB.
* ``FINITE`` sources: each magnitude is a finite rupture with the focal mechanism dip (90
  degrees for strike slip, 50 degrees for reverse and normal) and a down-dip width from the
  rupture scaling relation, limited by the maximum depth (crustal grids) or a maximum width
  (slab grids). Strike-slip ruptures are on the footwall (rX = -rJB, rRup = hypot(rJB,
  zTor)). Reverse and normal ruptures are split into a footwall and a hanging wall rupture
  with half of the rate each. For the hanging wall rupture rX = rJB + horizontal width and
  rRup is hypot(rJB, zBot) beyond rCut = zBot tan(dip), and otherwise interpolated linearly
  in rJB between min(hypot(widthH, zTor), zBot cos(dip)) at rJB = 0 and zBot / cos(dip) at
  rCut.
* Grid optimization (``GridRuptureSet.Table``): when a distance bin is given, the horizontal
  distances are replaced by the centers of bins of that width. Nodes closer than the
  smoothing limit are replaced by ``density`` x ``density`` nodes offset by fractions of the
  grid spacing (in degrees), each with 1 / density**2 of the rate.

The depth to the top of rupture is the node depth for subduction (slab) grids and comes
from the magnitude-depth model for crustal grids (see :func:`expand_depths`).
"""

from importlib.resources import files

import numpy as np

EARTH_RADIUS_MEAN = 6371.0072  # km, nshmp-lib Coordinates.EARTH_RADIUS_MEAN

# nshmp-lib FocalMech: (dip in degrees) by ucla_plha style (1 reverse, 2 normal, 3 strike slip)
MECH_DIP = {1: 50.0, 2: 50.0, 3: 90.0}

RUPTURE_SCALINGS = [
    "nshm_point_wc94_length",
    "nshm_sub_geomat_length",
    "nshm_somerville",
    "peer",
    "none",
]

_TABLES = {}


def _table(scaling):
    if not _TABLES:
        data = np.load(str(files("ucla_plha").joinpath("geometry/nshmp_rjb_tables.npz")))
        for key in data.files:
            _TABLES[key] = data[key]
    return _TABLES[scaling]


def xyz_to_latlon(xyz):
    """Return latitude and longitude (degrees) of points from geometry.point_to_xyz.

    point_to_xyz scales the unit vector (cos(lat) cos(lon), cos(lat) sin(lon), sin(lat)) by
    a radius, so the latitude and longitude are recovered exactly from the direction.
    """
    xyz = np.atleast_2d(np.asarray(xyz, dtype=float))
    lat = np.degrees(np.arctan2(xyz[:, 2], np.hypot(xyz[:, 0], xyz[:, 1])))
    lon = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    return lat, lon


def horz_distance_fast(lat1, lon1, lat2, lon2):
    """nshmp-lib Locations.horzDistanceFast (km); latitudes and longitudes in degrees."""
    lat1, lon1, lat2, lon2 = (np.radians(np.asarray(v, dtype=float)) for v in (lat1, lon1, lat2, lon2))
    d_lat = lat1 - lat2
    d_lon = (lon1 - lon2) * np.cos((lat1 + lat2) * 0.5)
    return EARTH_RADIUS_MEAN * np.sqrt(d_lat * d_lat + d_lon * d_lon)


def point_source_distance(m, r, scaling):
    """nshmp-lib RuptureScaling.pointSourceDistance, vectorized.

    Args:
        m (array): magnitudes
        r (array): horizontal distances from the site to the point sources (km)
        scaling (str): rupture scaling name (lower case), one of RUPTURE_SCALINGS

    Returns:
        array: Joyner-Boore distance (km)
    """
    m = np.asarray(m, dtype=float)
    r = np.asarray(r, dtype=float)
    scaling = scaling.lower()
    if scaling == "none":
        return r.copy()
    if scaling == "peer":
        return np.maximum(0.5, r)
    table = _table(scaling)
    out = np.array(r, dtype=float, copy=True)
    use = m >= 6.0
    # Java Math.round(x) = floor(x + 0.5)
    i_m = np.minimum(np.floor((m[use] - 6.05) / 0.1 + 0.5).astype(int), 25)
    i_r = np.minimum(1000, np.floor(r[use]).astype(int))
    out[use] = table[i_m, i_r]
    return out


def rupture_width(m, max_width, scaling):
    """Down-dip width (km) from nshmp-lib RuptureScaling.dimensions(mag, maxWidth).width."""
    m = np.asarray(m, dtype=float)
    max_width = np.asarray(max_width, dtype=float)
    scaling = scaling.lower()
    if scaling == "nshm_point_wc94_length":
        length = 10.0 ** (-3.22 + 0.69 * m)
        return np.minimum(max_width, length / 1.5)
    if scaling == "nshm_fault_wc94_length":
        length = 10.0 ** (-3.22 + 0.69 * m)
        return np.minimum(max_width, length)
    if scaling == "nshm_somerville":
        width = np.sqrt(10.0 ** (m - 4.366))
        return np.where(width < max_width, width, max_width)
    if scaling == "nshm_sub_geomat_length":
        return np.broadcast_to(max_width, m.shape).astype(float)
    if scaling == "peer":
        return np.minimum(max_width, np.sqrt(10.0 ** (m - 4.0) / 2.0))
    raise ValueError(f'rupture scaling "{scaling}" has no rupture dimensions')


def smooth_nodes(site_lat, site_lon, node_lat, node_lon, density, limit, spacing):
    """nshmp-lib grid smoothing (GridRuptureSet.addSmoothed) for nodes near a site.

    Returns (index, lat, lon, scale): for every output node, the index of the input node,
    its location, and the factor applied to its rates. Nodes at a horizontal distance of
    ``limit`` km or more are returned unchanged (scale 1); closer nodes are replaced by
    density**2 nodes offset by the nshmp-lib offsets times the grid spacing (degrees).
    """
    if density == 10:
        offsets = np.arange(-0.45, 0.45 + 1e-9, 0.1)
    elif density == 4:
        offsets = np.array([-0.375, -0.125, 0.125, 0.375])
    else:
        raise ValueError("nshmp-lib grid smoothing density must be 4 or 10")
    offsets = np.round(spacing * offsets, 5)
    node_lat = np.asarray(node_lat, dtype=float)
    node_lon = np.asarray(node_lon, dtype=float)
    r = horz_distance_fast(site_lat, site_lon, node_lat, node_lon)
    near = np.flatnonzero(r < limit)
    far = np.flatnonzero(r >= limit)
    if len(near) == 0:
        return far, node_lat[far], node_lon[far], np.ones(len(far))
    d_lat, d_lon = np.meshgrid(offsets, offsets, indexing="ij")
    d_lat = d_lat.ravel()
    d_lon = d_lon.ravel()
    n = len(d_lat)
    index = np.concatenate([far, np.repeat(near, n)])
    lat = np.concatenate([node_lat[far], (node_lat[near][:, None] + d_lat).ravel()])
    lon = np.concatenate([node_lon[far], (node_lon[near][:, None] + d_lon).ravel()])
    scale = np.concatenate([np.ones(len(far)), np.full(len(near) * n, 1.0 / n)])
    return index, lat, lon, scale


def expand_depths(m, depth_map):
    """Expand ruptures over a magnitude-depth model (nshmp-lib grid-depth-map).

    Args:
        m (array): magnitudes, length N
        depth_map (list): entries {"m_min": .., "m_max": .., "depths": [...], "weights": [...]};
            a rupture uses the entry with m_min <= m < m_max (nshmp-lib uses the magnitude
            cutoffs as m < mMax)

    Returns:
        (index, ztor, weight): index into the input ruptures, depth to top of rupture, and
        depth weight of every expanded rupture
    """
    m = np.asarray(m, dtype=float)
    index, ztor, weight = [], [], []
    assigned = np.zeros(len(m), dtype=bool)
    for entry in depth_map:
        sel = np.flatnonzero((m >= entry["m_min"]) & (m < entry["m_max"]) & ~assigned)
        assigned[sel] = True
        for depth, w in zip(entry["depths"], entry["weights"]):
            index.append(sel)
            ztor.append(np.full(len(sel), float(depth)))
            weight.append(np.full(len(sel), float(w)))
    if not np.all(assigned):
        raise ValueError("grid depth map does not cover all rupture magnitudes")
    index = np.concatenate(index)
    order = np.argsort(index, kind="stable")
    return index[order], np.concatenate(ztor)[order], np.concatenate(weight)[order]


def finite_point_source_distances(
    m, style, r_horizontal, ztor, scaling, source_type="finite", max_depth=None, max_width=None
):
    """Distances and geometry of nshmp-lib POINT or FINITE point sources.

    Args:
        m (array): magnitudes, length N
        style (array): 1 reverse, 2 normal, 3 strike slip, length N
        r_horizontal (array): horizontal site-to-node distances (km), length N
        ztor (array): depth to top of rupture (km), length N
        scaling (str): rupture scaling (RUPTURE_SCALINGS)
        source_type (str): "point" or "finite" (nshmp-lib point-source-type)
        max_depth (float): maximum depth of finite ruptures (crustal grids), or
        max_width (float): maximum horizontal-equivalent width (slab grids); exactly one of
            max_depth and max_width is used for "finite" sources

    Returns:
        dict of arrays of length K (K >= N; reverse and normal finite sources are split into
        footwall and hanging wall ruptures): "index" (input rupture of each output rupture),
        "rate_scale" (0.5 for the footwall and hanging wall halves, 1 otherwise), "rjb",
        "rrup", "rx", "rx1", "ry0", "dip", "ztor", "zbor", "width"
    """
    m = np.asarray(m, dtype=float)
    style = np.asarray(style)
    r_horizontal = np.asarray(r_horizontal, dtype=float)
    ztor = np.asarray(ztor, dtype=float)
    n = len(m)
    dip = np.vectorize(MECH_DIP.get, otypes=[float])(style) if n else np.empty(0)
    dip_rad = np.radians(dip)
    rjb = point_source_distance(m, r_horizontal, scaling)
    source_type = source_type.lower()

    if source_type == "point":
        # GridSource.PointSurface: no finiteness; width is a generic 10 km in nshmp-lib
        rrup = np.hypot(rjb, ztor)
        width = np.full(n, 10.0)
        zbor = ztor + width * np.sin(dip_rad)
        return {
            "index": np.arange(n),
            "rate_scale": np.ones(n),
            "rjb": rjb,
            "rrup": rrup,
            "rx": rjb.copy(),
            "rx1": rjb - width * np.cos(dip_rad),
            "ry0": np.zeros(n),
            "dip": dip,
            "ztor": ztor,
            "zbor": zbor,
            "width": width,
        }
    if source_type != "finite":
        raise ValueError(f'point source type "{source_type}" is not supported')

    if (max_depth is None) == (max_width is None):
        raise ValueError("finite point sources need exactly one of max_depth and max_width")
    sin_dip = np.sin(dip_rad)
    cos_dip = np.cos(dip_rad)
    if max_depth is not None:
        max_width_dd = (max_depth - ztor) / sin_dip
    else:
        max_width_dd = max_width / sin_dip
    width_dd = rupture_width(m, max_width_dd, scaling)
    width_h = width_dd * cos_dip
    zbot = ztor + width_dd * sin_dip

    # Strike-slip ruptures (footwall), then footwall and hanging wall halves of the others
    dip_slip = np.flatnonzero(style != 3)
    index = np.concatenate([np.arange(n), dip_slip])
    hanging = np.concatenate([np.zeros(n, dtype=bool), np.ones(len(dip_slip), dtype=bool)])
    rate_scale = np.where(style[index] != 3, 0.5, 1.0)

    rjb_k = rjb[index]
    ztor_k = ztor[index]
    zbot_k = zbot[index]
    width_h_k = width_h[index]
    dip_k = dip_rad[index]
    rx = np.where(hanging, rjb_k + width_h_k, -rjb_k)
    rrup = np.hypot(rjb_k, ztor_k)
    if np.any(hanging):
        h = hanging
        r_cut = zbot_k[h] * np.tan(dip_k[h])
        rrup0 = np.minimum(np.hypot(width_h_k[h], ztor_k[h]), zbot_k[h] * np.cos(dip_k[h]))
        rrup_c = zbot_k[h] / np.cos(dip_k[h])
        with np.errstate(divide="ignore", invalid="ignore"):
            inner = (rrup_c - rrup0) * rjb_k[h] / r_cut + rrup0
        rrup[h] = np.where(rjb_k[h] > r_cut, np.hypot(rjb_k[h], zbot_k[h]), inner)
    return {
        "index": index,
        "rate_scale": rate_scale,
        "rjb": rjb_k,
        "rrup": rrup,
        "rx": rx,
        "rx1": rx - width_h_k,
        "ry0": np.zeros(len(index)),
        "dip": np.degrees(dip_k),
        "ztor": ztor_k,
        "zbor": zbot_k,
        "width": width_dd[index],
    }


def location_at(lat, lon, azimuth_deg, distance):
    """nshmp-lib Locations.location: point at a distance (km) along an azimuth (degrees)."""
    lat1 = np.radians(np.asarray(lat, dtype=float))
    lon1 = np.radians(np.asarray(lon, dtype=float))
    az = np.radians(np.asarray(azimuth_deg, dtype=float))
    ad = np.asarray(distance, dtype=float) / EARTH_RADIUS_MEAN
    lat2 = np.arcsin(np.sin(lat1) * np.cos(ad) + np.cos(lat1) * np.sin(ad) * np.cos(az))
    lon2 = lon1 + np.arctan2(
        np.sin(az) * np.sin(ad) * np.cos(lat1), np.cos(ad) - np.sin(lat1) * np.sin(lat2)
    )
    return np.degrees(lat2), np.degrees(lon2)


def _segment_frame(lat1, lon1, lat2, lon2, site_lat, site_lon):
    p1, l1, p2, l2, p3, l3 = (
        np.radians(np.asarray(v, dtype=float)) for v in (lat1, lon1, lat2, lon2, site_lat, site_lon)
    )
    scale = np.cos(0.5 * p3 + 0.25 * p1 + 0.25 * p2)
    x2 = (l2 - l1) * scale
    y2 = p2 - p1
    x3 = (l3 - l1) * scale
    y3 = p3 - p1
    return x2, y2, x3, y3


def distance_to_line_fast(lat1, lon1, lat2, lon2, site_lat, site_lon):
    """nshmp-lib Locations.distanceToLineFast (signed, km)."""
    x2, y2, x3, y3 = _segment_frame(lat1, lon1, lat2, lon2, site_lat, site_lon)
    return (x3 * y2 - x2 * y3) / np.sqrt(x2 * x2 + y2 * y2) * EARTH_RADIUS_MEAN


def distance_to_segment_fast(lat1, lon1, lat2, lon2, site_lat, site_lon):
    """nshmp-lib Locations.distanceToSegmentFast (java.awt.geom.Line2D.ptSegDist, km)."""
    x2, y2, x3, y3 = _segment_frame(lat1, lon1, lat2, lon2, site_lat, site_lon)
    length2 = x2 * x2 + y2 * y2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(length2 > 0, (x3 * x2 + y3 * y2) / length2, 0.0)
    t = np.clip(t, 0.0, 1.0)
    return np.hypot(x3 - t * x2, y3 - t * y2) * EARTH_RADIUS_MEAN


def fixed_strike_distances(m, style, node_lat, node_lon, strike, site_lat, site_lon, ztor,
                           scaling, max_depth=None, max_width=None):
    """Distances of nshmp-lib FIXED_STRIKE point sources (GridSourceFixedStrike), strike slip.

    The rupture is vertical with a top trace of length L (rupture scaling relation) centered
    on the node along the strike; rJB = distanceToSegmentFast, rRup = hypot(rJB, zTor),
    rX = distanceToLineFast. No point-source distance correction is applied. Only
    strike-slip (vertical) ruptures are supported, which is what the NSHM zones use.

    Returns a dict like finite_point_source_distances.
    """
    m = np.asarray(m, dtype=float)
    style = np.asarray(style)
    if np.any(style != 3):
        raise ValueError("FIXED_STRIKE point sources are implemented for strike-slip ruptures only")
    ztor = np.asarray(ztor, dtype=float)
    n = len(m)
    if max_depth is not None:
        max_width_dd = max_depth - ztor
    else:
        max_width_dd = np.full(n, float(max_width))
    scaling = scaling.lower()
    if scaling not in ("nshm_point_wc94_length", "nshm_fault_wc94_length"):
        raise ValueError(f'FIXED_STRIKE rupture length for "{scaling}" is not implemented')
    length = 10.0 ** (-3.22 + 0.69 * m)
    width = rupture_width(m, max_width_dd, scaling)
    lat1, lon1 = location_at(node_lat, node_lon, strike, length / 2.0)
    lat2, lon2 = location_at(node_lat, node_lon, np.asarray(strike) + 180.0, length / 2.0)
    rjb = distance_to_segment_fast(lat1, lon1, lat2, lon2, site_lat, site_lon)
    rx = distance_to_line_fast(lat1, lon1, lat2, lon2, site_lat, site_lon)
    return {
        "index": np.arange(n),
        "rate_scale": np.ones(n),
        "rjb": rjb,
        "rrup": np.hypot(rjb, ztor),
        "rx": rx,
        "rx1": rx,
        "ry0": np.zeros(n),
        "dip": np.full(n, 90.0),
        "ztor": ztor,
        "zbor": ztor + width,
        "width": width,
    }

import numpy as np


def point_to_xyz(p):
    """
    Convert point latitude, longitude, elevation values to Cartesian coordinates.
    p = N x 3 Numpy array, where N is the number of points.
    Each row of p has [latitude, longitude, elevation], where latitude and longitude are in degrees and elevation is in km.
    r = Earth's radius in km.
    """
    p_xyz = np.empty(3)
    rad = np.pi / 180.0
    lat = p[0]
    lon = p[1]
    elevation = p[2]
    a = 6378.1370  # Earth's equatorial radius in km
    b = 6356.7523  # Earth's polar radius in km
    r = np.sqrt(
        ((a**2 * np.cos(lat * rad)) ** 2 + (b**2 * np.sin(lat * rad)) ** 2)
        / ((a * np.cos(lat * rad)) ** 2 + (b * np.sin(lat * rad)) ** 2)
    )
    p_xyz[0] = (r + elevation) * np.cos(lat * rad) * np.cos(lon * rad)
    p_xyz[1] = (r + elevation) * np.cos(lat * rad) * np.sin(lon * rad)
    p_xyz[2] = (r + elevation) * np.sin(lat * rad)
    return np.asarray(p_xyz)


def point_triangle_distance(tri_xyz, p_xyz, tri_segment_id):
    """
    Compute distance between all of the triangles for a source model (tri_xyz) and the point p_xyz.
    Return a Numpy array of distances and associated fault_id's.
    tri_xyz = N x 3 x 3 Numpy array of points defining triangle in Cartesian coordinates
    p_xyz = Numpy array defining point in Cartesian coordinates
    fault_id = N x 1 Numpy array of integers defining the fault id for each triangle
    """
    # Perform checks that all three points are different. Add 1m to elevation for points that have the same coordinates.
    filt1 = (tri_xyz[:, 0] == tri_xyz[:, 1]).all(axis=1)
    filt2 = (tri_xyz[:, 0] == tri_xyz[:, 2]).all(axis=1)
    filt3 = (tri_xyz[:, 1] == tri_xyz[:, 2]).all(axis=1)
    tri_xyz[filt1, 0, 2] += 0.001
    tri_xyz[filt2, 0, 2] += 0.001
    tri_xyz[filt3, 1, 2] += 0.001

    e0 = tri_xyz[:, 1] - tri_xyz[:, 0]
    e1 = tri_xyz[:, 2] - tri_xyz[:, 0]
    a = np.sum(np.multiply(e0, e0), axis=1)
    b = np.sum(np.multiply(e0, e1), axis=1)
    c = np.sum(np.multiply(e1, e1), axis=1)

    det = a * c - b * b
    det[det == 0] = 1.0e-8

    d0 = tri_xyz[:, 0] - p_xyz

    d = np.sum(np.multiply(e0, d0), axis=1)
    e = np.sum(np.multiply(e1, d0), axis=1)
    f = np.sum(np.multiply(d0, d0), axis=1)

    s = b * e - c * d
    t = b * d - a * e

    sqrdistance = np.empty(len(tri_xyz), dtype=np.float64)

    # Region 4
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d < 0.0) & (-d >= a)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d < 0.0) & (-d < a)
    sqrdistance[cond] = (-d * d / a + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d >= 0.0) & (e >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d >= 0.0) & (e < 0.0) & (-e >= c)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d >= 0.0) & (e < 0.0) & (-e < c)
    sqrdistance[cond] = (-e * e / c + f)[cond]

    # Region 3
    cond = (s + t <= det) & (s < 0.0) & (t >= 0.0) & (e >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t >= 0.0) & (e < 0.0) & (-e >= c)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t >= 0.0) & (e < 0.0) & (-e < c)
    sqrdistance[cond] = (-e * e / c + f)[cond]

    # Region 5
    cond = (s + t <= det) & (s >= 0.0) & (t < 0.0) & (d >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t <= det) & (s >= 0.0) & (t < 0.0) & (d < 0.0) & (-d >= a)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t <= det) & (s >= 0.0) & (t < 0.0) & (d < 0.0) & (-d < a)
    sqrdistance[cond] = (-d * d / a + f)[cond]

    # Region 0
    invDet = 1.0 / det
    stemp = s * invDet
    ttemp = t * invDet
    cond = (s + t <= det) & (s >= 0) & (t >= 0)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]

    # Region 2
    tmp0 = b + d
    tmp1 = c + e
    numer = tmp1 - tmp0
    denom = a - 2.0 * b + c
    denom[denom == 0] = 1.0e-8
    stemp = numer / denom
    ttemp = 1.0 - stemp
    cond = (s + t > det) & (s < 0.0) & (tmp1 > tmp0) & (numer >= denom)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 > tmp0) & (numer < denom)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 <= tmp0) & (tmp1 <= 0.0)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 <= tmp0) & (tmp1 > 0.0) & (e >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 <= tmp0) & (tmp1 > 0.0) & (e < 0.0)
    sqrdistance[cond] = (-e * e / c + f)[cond]

    # Region 6
    tmp0 = b + e
    tmp1 = a + d
    numer = tmp1 - tmp0
    denom = a - 2.0 * b + c
    denom[denom == 0] = 1.0e-8
    ttemp = numer / denom
    stemp = 1.0 - ttemp
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 > tmp0) & (numer >= denom)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 > tmp0) & (numer < denom)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 <= tmp0) & (tmp1 <= 0)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 <= tmp0) & (tmp1 > 0) & (d >= 0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 <= tmp0) & (tmp1 > 0) & (d < 0)
    sqrdistance[cond] = (-d * d / a + f)[cond]

    # Region 1
    numer = c + e - b - d
    denom = a - 2.0 * b + c
    denom[denom == 0] = 1.0e-8
    stemp = numer / denom
    ttemp = 1.0 - stemp
    cond = (s + t > det) & (s >= 0) & (t >= 0) & (numer <= 0)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t >= 0) & (numer > 0) & (numer >= denom)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t >= 0) & (numer > 0) & (numer < denom)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]

    # account for numerical round-off error
    sqrdistance[sqrdistance <= 0] = 0

    sort_idx = np.argsort(tri_segment_id)
    tri_segment_id_sorted = tri_segment_id[sort_idx]
    sqrdistance_sorted = sqrdistance[sort_idx]

    split_indices = np.diff(tri_segment_id_sorted) != 0
    boundaries = np.r_[0, np.where(split_indices)[0] + 1]
    sqrdistance_out = np.minimum.reduceat(sqrdistance_sorted, boundaries)

    # return distance array
    return np.sqrt(sqrdistance_out).T


def get_Rx_Rx1_Ry0(rect_points, point, rect_segment_id):
    """
    rect_points is an Nx4x3 Numpy array with the x, y, z coordinates of four points on the surface projection of the fault, and N is the number of faults
    point is a 1x3 numpy array with the x, y, z coordinates of the point
    top edge of fault is defined by (x1,y1,z1) (x2,y2,z2), and bottom edge by (x3,y3,z3) (x4,y4,z4)

                         Rx            point
              |<---------------------->o(x0,y0,z0)
              |                        ^
              |               Rx1      |
              |        |<------------->| Ry0
              |        |               |
              |        |               v
    (x2,y2,z2)o--------o(x4,y4,z4)-------
              |////////|
              |////////|<--Surface projection of fault
        top-->|////////|<--bottom
              |////////|
              |////////|
    (x1,y1,z1)o--------o(x3,y3,z3)

    """

    width = np.sqrt(np.sum((rect_points[:, 2] - rect_points[:, 0]) ** 2, axis=1))
    width[width < 0.001] = 0.001
    length = np.sqrt(np.sum((rect_points[:, 1] - rect_points[:, 0]) ** 2, axis=1))
    Rx = (
        np.sqrt(
            np.sum(
                (
                    np.cross(
                        rect_points[:, 1] - rect_points[:, 0], rect_points[:, 0] - point
                    )
                )
                ** 2,
                axis=1,
            )
        )
        / length
    )
    # Rx is positive on the hanging wall side of the top edge, which is the side containing the
    # bottom edge. Sites on the top edge line and sites near vertical faults have positive Rx.
    # Rx1 = Rx - W cos(dip), where W cos(dip) is the horizontal distance between the top and
    # bottom edges, so Rx1 is negative for sites above the rupture.
    strike = rect_points[:, 1] - rect_points[:, 0]
    site_side = np.cross(strike, point - rect_points[:, 0])
    bottom_side = np.cross(strike, rect_points[:, 2] - rect_points[:, 0])
    footwall = np.sum(site_side * bottom_side, axis=1) < 0.0
    Rx[footwall] = -Rx[footwall]
    Rx1 = Rx - np.sqrt(np.sum(bottom_side**2, axis=1)) / length
    Ry0a = (
        np.sqrt(
            np.sum(
                (
                    np.cross(
                        rect_points[:, 1] - rect_points[:, 3], rect_points[:, 1] - point
                    )
                )
                ** 2,
                axis=1,
            )
        )
        / width
    )
    Ry0b = (
        np.sqrt(
            np.sum(
                (
                    np.cross(
                        rect_points[:, 0] - rect_points[:, 2], rect_points[:, 0] - point
                    )
                )
                ** 2,
                axis=1,
            )
        )
        / width
    )
    Ry0 = np.empty(len(rect_points), dtype=np.float64)
    Ry0[Ry0a < Ry0b] = Ry0a[Ry0a < Ry0b]
    Ry0[Ry0b <= Ry0a] = Ry0b[Ry0b <= Ry0a]
    Ry0[(Ry0a < length) & (Ry0b < length)] = 0

    # A "segment" in the UCERF3 model sometimes consists of multiple geometric objects, so we need to
    # compute the shortest distance between each geometric object and the point to find the shortest distance
    # to the "segment"

    split_indices = np.diff(rect_segment_id) != 0
    boundaries = np.r_[0, np.where(split_indices)[0] + 1]
    Rx_out = np.minimum.reduceat(Rx, boundaries)
    Rx1_out = np.minimum.reduceat(Rx1, boundaries)
    Ry0_out = np.minimum.reduceat(Ry0, boundaries)

    return (Rx_out, Rx1_out, Ry0_out)


# ------------------------------------------------------------------------------------------------
# nshmp-lib distance conventions
# ------------------------------------------------------------------------------------------------

NSHMP_EARTH_RADIUS = 6371.0072  # km, nshmp-lib Coordinates.EARTH_RADIUS_MEAN


def _group_boundaries(group_id):
    """First index of each run of equal values in group_id (which must be grouped)."""
    return np.r_[0, np.flatnonzero(np.diff(group_id) != 0) + 1]


def xyz_to_lon_lat(xyz):
    """Longitude and latitude (degrees) of Cartesian points (inverse of the angles of
    point_to_xyz)."""
    lon = np.degrees(np.arctan2(xyz[..., 1], xyz[..., 0]))
    lat = np.degrees(np.arctan2(xyz[..., 2], np.hypot(xyz[..., 0], xyz[..., 1])))
    return lon, lat


def _nshmp_azimuth(lon1, lat1, lon2, lat2):
    """Locations.azimuthRad (radians, degrees in)."""
    lon1, lat1, lon2, lat2 = (np.radians(v) for v in (lon1, lat1, lon2, lat2))
    dlon = lon2 - lon1
    az = np.arctan2(np.sin(dlon) * np.cos(lat2),
                    np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon))
    return (az + 2 * np.pi) % (2 * np.pi)


def _nshmp_location(lon, lat, azimuth, distance):
    """Locations.location: the point distance (km) from (lon, lat) along azimuth (radians)."""
    lat1 = np.radians(lat)
    lon1 = np.radians(lon)
    ad = distance / NSHMP_EARTH_RADIUS
    lat2 = np.arcsin(np.sin(lat1) * np.cos(ad) + np.cos(lat1) * np.sin(ad) * np.cos(azimuth))
    lon2 = lon1 + np.arctan2(np.sin(azimuth) * np.sin(ad) * np.cos(lat1),
                             np.cos(ad) - np.sin(lat1) * np.sin(lat2))
    return np.degrees(lon2), np.degrees(lat2)


def _nshmp_segment_distance(lon1, lat1, lon2, lat2, site_lon, site_lat):
    """Locations.distanceToSegmentFast (km) from the site to the segments (degrees in)."""
    rad = np.pi / 180.0
    lon1, lat1, lon2, lat2 = lon1 * rad, lat1 * rad, lon2 * rad, lat2 * rad
    lon3, lat3 = site_lon * rad, site_lat * rad
    scale = np.cos(0.5 * lat3 + 0.25 * lat1 + 0.25 * lat2)
    x2 = (lon2 - lon1) * scale
    y2 = lat2 - lat1
    x3 = (lon3 - lon1) * scale
    y3 = lat3 - lat1
    l2 = x2 * x2 + y2 * y2
    with np.errstate(invalid="ignore", divide="ignore"):
        t = np.where(l2 > 0, (x3 * x2 + y3 * y2) / l2, 0.0)
    t = np.clip(t, 0.0, 1.0)
    dx = x3 - t * x2
    dy = y3 - t * y2
    return np.sqrt(dx * dx + dy * dy) * NSHMP_EARTH_RADIUS


def _mercator_y(lat):
    return np.degrees(np.log(np.tan(np.pi / 4 + np.radians(lat) / 2)))


def extended_trace_rx(lon1, lat1, lon2, lat2, boundaries, site_lon, site_lat, extension=1000.0):
    """Rx of polylines as nshmp-lib computes it (model.Distance.getDistanceX).

    Each polyline (e.g. the upper edge of a fault section or rupture) is given by its pieces,
    piece k going from (lon1[k], lat1[k]) to (lon2[k], lat2[k]) (degrees); the pieces of
    polyline i are boundaries[i]:boundaries[i + 1], in order along the polyline. With strike
    the azimuth from the first to the last point, nshmp-lib extends the polyline by P2 and P3,
    `extension` (1000) km from its ends along -strike and strike, and closes it by P1 and P4,
    1000 km from P2 and P3 along the dip direction (strike + 90 degrees). |Rx| is the
    distance from the site to the extended polyline P2, trace, P3 (Locations.minDistanceToLine,
    i.e. distanceToSegmentFast to each piece), positive if the site is inside the polygon P1,
    P2, trace, P3, P4 (straight edges in Mercator coordinates, even-odd rule) or on the
    polyline, and negative otherwise.
    """
    n = len(lon1)
    boundaries = np.asarray(boundaries, dtype=np.int64)
    last = np.r_[boundaries[1:], n] - 1
    flon, flat = lon1[boundaries], lat1[boundaries]
    llon, llat = lon2[last], lat2[last]
    strike = _nshmp_azimuth(flon, flat, llon, llat)
    dip_dir = strike + np.pi / 2
    p3 = _nshmp_location(llon, llat, strike, extension)
    p4 = _nshmp_location(p3[0], p3[1], dip_dir, extension)
    p2 = _nshmp_location(flon, flat, strike + np.pi, extension)
    p1 = _nshmp_location(p2[0], p2[1], dip_dir, extension)

    d = np.minimum.reduceat(_nshmp_segment_distance(lon1, lat1, lon2, lat2, site_lon, site_lat),
                            boundaries)
    d = np.minimum(d, _nshmp_segment_distance(p2[0], p2[1], flon, flat, site_lon, site_lat))
    d = np.minimum(d, _nshmp_segment_distance(llon, llat, p3[0], p3[1], site_lon, site_lat))

    sy = _mercator_y(site_lat)

    def crossings(ax, ay, bx, by):
        # horizontal ray from the site toward +x (Mercator coordinates)
        ay, by = _mercator_y(ay), _mercator_y(by)
        cond = (ay > sy) != (by > sy)
        with np.errstate(invalid="ignore", divide="ignore"):
            xint = ax + (sy - ay) * (bx - ax) / (by - ay)
        return (cond & (site_lon < xint)).astype(np.int64)

    c = np.add.reduceat(crossings(lon1, lat1, lon2, lat2), boundaries)
    for a, b in [(p1, p2), (p2, (flon, flat)), ((llon, llat), p3), (p3, p4), (p4, p1)]:
        c += crossings(a[0], a[1], b[0], b[1])
    inside = (c % 2) == 1
    return np.where(inside | (d == 0.0), d, -d)


def section_rx_extended_trace(rect_points, p_xyz, rect_segment_id):
    """Rx and Rx1 of each fault segment (section) with the nshmp-lib convention.

    rect_points and rect_segment_id are as in get_Rx_Rx1_Ry0 (the pieces of a segment must be
    consecutive and in order along its upper edge). Rx is that of extended_trace_rx for the
    upper edge of the segment (all of its pieces), and Rx1 = Rx - the mean horizontal width of
    its pieces. Returns (Rx, Rx1) per segment, in the order of the segment ids.
    """
    lon, lat = xyz_to_lon_lat(rect_points[:, :2])
    site_lon, site_lat = xyz_to_lon_lat(np.asarray(p_xyz, dtype=float))
    boundaries = _group_boundaries(rect_segment_id)
    rx = extended_trace_rx(lon[:, 0], lat[:, 0], lon[:, 1], lat[:, 1], boundaries,
                           site_lon, site_lat)
    strike = rect_points[:, 1] - rect_points[:, 0]
    length = np.sqrt(np.sum(strike**2, axis=1))
    length[length == 0] = 1.0
    wh = np.sqrt(np.sum(np.cross(strike, rect_points[:, 2] - rect_points[:, 0]) ** 2, axis=1)) / length
    counts = np.diff(np.r_[boundaries, len(rect_segment_id)])
    wh = np.add.reduceat(wh, boundaries) / counts
    return rx, rx - wh


def _range_reduce(ufunc, values, start, stop, fill):
    """ufunc.reduce of values[start[i]:stop[i]] for each i (fill for empty ranges)."""
    padded = np.r_[values, np.asarray([fill], dtype=values.dtype)]
    idx = np.empty(2 * len(start), dtype=np.int64)
    idx[0::2] = start
    idx[1::2] = np.minimum(stop, len(values))
    out = ufunc.reduceat(padded, idx)[0::2]
    return np.where(stop > start, out, fill)


def _inside_quads(px, py, qx, qy, rx_, ry_, sx_, sy_, x, y):
    """Even-odd test of the point (x, y) in the quadrilaterals p, q, r, s (arrays of vertices)."""
    inside = np.zeros(len(px), dtype=bool)
    verts = [(px, py), (qx, qy), (rx_, ry_), (sx_, sy_)]
    for k in range(4):
        ax, ay = verts[k]
        bx, by = verts[(k + 1) % 4]
        cond = (ay > y) != (by > y)
        with np.errstate(invalid="ignore", divide="ignore"):
            xint = ax + (y - ay) * (bx - ax) / (by - ay)
        inside ^= cond & (x < xint)
    return inside


def gridded_rupture_distances(grid, site_lat, site_lon, need_rx=True):
    """rJB, rRup, and rX of ruptures on gridded surfaces as nshmp-lib computes them
    (model.Distance.compute).

    grid: dict-like with the arrays of grid.npz: lon, lat, depth (grid points of all surfaces,
    surface by surface, column by column: point (r, c) of surface s is
    surface_offset[s] + c * surface_rows[s] + r), surface_offset, surface_rows, surface_cols,
    surface_dip (degrees), surface_spacing (average grid spacing, km), and per rupture
    rupture_surface, rupture_row0, rupture_rows, rupture_col0, rupture_cols (the subset of the
    grid of its surface, nshmp-lib GriddedSubsetSurface).

    rJB and rRup are the minimum horizontal (Locations.horzDistanceFast, R = 6371.0072 km) and
    the minimum sqrt(horizontal^2 + depth^2) distances to the grid points of the rupture (top
    row only if dip > 89 degrees). rJB is 0 if it is less than the grid spacing and the site is
    inside the perimeter of the rupture (top and bottom rows, longitude-latitude polygon). rX
    is extended_trace_rx of the top row of the rupture. Returns (rjb, rrup, rx); rx is None if
    need_rx is False.
    """
    rad = np.pi / 180.0
    lon = np.asarray(grid["lon"], dtype=np.float64)
    lat = np.asarray(grid["lat"], dtype=np.float64)
    depth = np.asarray(grid["depth"], dtype=np.float64)
    dlat = (lat - site_lat) * rad
    dlon = (lon - site_lon) * rad * np.cos(0.5 * (lat + site_lat) * rad)
    h = NSHMP_EARTH_RADIUS * np.sqrt(dlat * dlat + dlon * dlon)
    r2 = h * h + depth * depth

    offset = np.asarray(grid["surface_offset"], dtype=np.int64)
    rows = np.asarray(grid["surface_rows"], dtype=np.int64)
    cols = np.asarray(grid["surface_cols"], dtype=np.int64)
    dip = np.asarray(grid["surface_dip"], dtype=np.float64)
    spacing = np.asarray(grid["surface_spacing"], dtype=np.float64)
    rs = np.asarray(grid["rupture_surface"], dtype=np.int64)
    r0 = np.asarray(grid["rupture_row0"], dtype=np.int64)
    nr = np.asarray(grid["rupture_rows"], dtype=np.int64)
    c0 = np.asarray(grid["rupture_col0"], dtype=np.int64)
    nc = np.asarray(grid["rupture_cols"], dtype=np.int64)

    # column minima for each (surface, first row, number of rows) used by the ruptures
    groups, group_index = np.unique(np.c_[rs, r0, nr], axis=0, return_inverse=True)
    group_index = np.ravel(group_index)
    gstart = np.zeros(len(groups), dtype=np.int64)
    hmin, r2min, quad_in = [], [], []
    pos = 0
    for g, (s, a, n) in enumerate(groups):
        gstart[g] = pos
        block = slice(offset[s], offset[s] + rows[s] * cols[s])
        hb = h[block].reshape(cols[s], rows[s])
        rb = r2[block].reshape(cols[s], rows[s])
        n_eff = 1 if dip[s] > 89.0 else n
        hmin.append(hb[:, a:a + n_eff].min(axis=1))
        r2min.append(rb[:, a:a + n_eff].min(axis=1))
        # the perimeter of a rupture encloses the quadrilaterals between adjacent columns
        lo = lon[block].reshape(cols[s], rows[s])
        la = lat[block].reshape(cols[s], rows[s])
        tx, ty = lo[:, a], la[:, a]
        bx, by = lo[:, a + n - 1], la[:, a + n - 1]
        q = _inside_quads(tx[:-1], ty[:-1], tx[1:], ty[1:], bx[1:], by[1:], bx[:-1], by[:-1],
                          site_lon, site_lat)
        quad_in.append(np.r_[q, False])  # one entry per column
        pos += cols[s]
    hmin = np.concatenate(hmin)
    r2min = np.concatenate(r2min)
    quad_in = np.concatenate(quad_in).astype(np.int8)
    start = gstart[group_index] + c0
    rjb = _range_reduce(np.minimum, hmin, start, start + nc, np.inf)
    rrup = np.sqrt(_range_reduce(np.minimum, r2min, start, start + nc, np.inf))
    inside = _range_reduce(np.maximum, quad_in, start, start + nc - 1, 0) > 0
    rjb = np.where(inside & (rjb < spacing[rs]), 0.0, rjb)

    rx = None
    if need_rx:
        n_pieces = nc - 1
        rup_of_piece = np.repeat(np.arange(len(nc)), n_pieces)
        first = np.r_[0, np.cumsum(n_pieces)[:-1]]
        col = c0[rup_of_piece] + (np.arange(n_pieces.sum()) - first[rup_of_piece])
        s = rs[rup_of_piece]
        i1 = offset[s] + col * rows[s] + r0[rup_of_piece]
        i2 = i1 + rows[s]
        rx = extended_trace_rx(lon[i1], lat[i1], lon[i2], lat[i2], first, site_lon, site_lat)
    return rjb, rrup, rx


def _nshmp_horz_distance(lon1, lat1, lon2, lat2):
    """Locations.horzDistance (haversine, km; degrees in)."""
    lon1, lat1, lon2, lat2 = (np.radians(v) for v in (lon1, lat1, lon2, lat2))
    s1 = np.sin((lat2 - lat1) / 2.0)
    s2 = np.sin((lon2 - lon1) / 2.0)
    c = s1 * s1 + np.cos(lat1) * np.cos(lat2) * s2 * s2
    return NSHMP_EARTH_RADIUS * 2.0 * np.arctan2(np.sqrt(c), np.sqrt(1 - c))


def _nshmp_horz_distance_fast(lon1, lat1, lon2, lat2):
    """Locations.horzDistanceFast (km; degrees in)."""
    rad = np.pi / 180.0
    dlat = (lat1 - lat2) * rad
    dlon = (lon1 - lon2) * rad * np.cos(0.5 * (lat1 + lat2) * rad)
    return NSHMP_EARTH_RADIUS * np.sqrt(dlat * dlat + dlon * dlon)


def _java_round(x):
    return np.floor(x + 0.5).astype(np.int64)


def nshmp_section_grids(sections, spacing=1.0):
    """nshmp-lib DefaultGriddedSurfaces (with aseismicity) of fault sections.

    sections: dict-like with the arrays of sections.npz: trace_lon, trace_lat (the trace points
    of all sections, concatenated), trace_offset (first trace point of each section, with a last
    entry equal to the number of points), upper_depth, lower_depth, aseismicity, dip, and
    dip_direction (degrees) of each section.

    As DefaultGriddedSurface.createEvenlyGriddedSurface: the down-dip width (lower - upper) /
    sin(dip) is reduced by the aseismicity, which removes the top of the plane (the top row is
    moved down-dip from the trace by aseismicity x width cos(dip) horizontally and aseismicity x
    width sin(dip) vertically); the grid spacings are length / ceil(length / spacing) along strike
    (length from Locations.horzDistanceFast) and width / ceil(width / spacing) down dip; the
    columns are at multiples of the strike spacing along the trace (Locations.location along the
    trace segments, haversine segment lengths), and the rows go down dip along the dip direction.

    Returns a dict with lon, lat, depth (grid points, section by section, column by column, as in
    grid.npz), surface_offset, surface_rows, surface_cols, surface_dip, surface_spacing (average of
    the strike and dip spacings), and the perimeter of each section (DefaultGriddedSurface.perimeter,
    which uses the trace, not the top row: the trace and, in reverse order, the trace moved by the
    reduced width cos(dip) along the dip direction): perimeter_lon, perimeter_lat, perimeter_offset.
    """
    tlon = np.asarray(sections["trace_lon"], dtype=np.float64)
    tlat = np.asarray(sections["trace_lat"], dtype=np.float64)
    toff = np.asarray(sections["trace_offset"], dtype=np.int64)
    upper = np.asarray(sections["upper_depth"], dtype=np.float64)
    lower = np.asarray(sections["lower_depth"], dtype=np.float64)
    aseis = np.asarray(sections["aseismicity"], dtype=np.float64)
    dip = np.asarray(sections["dip"], dtype=np.float64)
    dip_dir = np.radians(np.asarray(sections["dip_direction"], dtype=np.float64))
    n_sec = len(upper)
    npts = np.diff(toff)
    nseg = npts - 1
    sec_of_pt = np.repeat(np.arange(n_sec), npts)

    dip_rad = np.radians(dip)
    width0 = (lower - upper) / np.sin(dip_rad)
    width = width0 * (1 - aseis)
    rv = np.sin(dip_rad) * width0 * aseis
    rh = np.cos(dip_rad) * width0 * aseis

    # trace segments (point k to k + 1 of the same section)
    is_seg = np.ones(len(tlon), dtype=bool)
    is_seg[toff[1:] - 1] = False
    a = np.flatnonzero(is_seg)
    seg_sec = sec_of_pt[a]
    seg_fast = _nshmp_horz_distance_fast(tlon[a], tlat[a], tlon[a + 1], tlat[a + 1])
    seg_len = _nshmp_horz_distance(tlon[a], tlat[a], tlon[a + 1], tlat[a + 1])
    seg_az = _nshmp_azimuth(tlon[a], tlat[a], tlon[a + 1], tlat[a + 1])
    seg_first = np.r_[0, np.cumsum(nseg)[:-1]]
    length = np.add.reduceat(seg_fast, seg_first)
    total = np.add.reduceat(seg_len, seg_first)
    cum = np.cumsum(seg_len)
    cum_before = cum - seg_len  # cumulative length at the start of each segment (all sections)
    sec_start_cum = cum_before[seg_first]

    strike_spacing = length / np.ceil(length / spacing)
    dip_spacing = width / np.ceil(width / spacing)
    rows = 1 + _java_round(width / dip_spacing)
    cols = 1 + _java_round(total / strike_spacing)

    # columns: point along the trace
    col_sec = np.repeat(np.arange(n_sec), cols)
    col_first = np.r_[0, np.cumsum(cols)[:-1]]
    c = np.arange(cols.sum()) - col_first[col_sec]
    along = c * strike_spacing[col_sec]
    # segment index: the first segment whose end is at or beyond `along` (the last if none)
    local_cum = cum - sec_start_cum[seg_sec]
    key = col_sec * 1.0e7 + along
    seg_key = seg_sec * 1.0e7 + local_cum
    k = np.searchsorted(seg_key, key, side="left")
    k = np.minimum(k, seg_first[col_sec] + nseg[col_sec] - 1)
    d = along - (local_cum[k] - seg_len[k])
    clon, clat = _nshmp_location(tlon[a[k]], tlat[a[k]], seg_az[k], d)
    shift = (rv > 0) | (rh > 0)
    s_ = shift[col_sec]
    slon, slat = _nshmp_location(clon, clat, dip_dir[col_sec], rh[col_sec])
    clon = np.where(s_, slon, clon)
    clat = np.where(s_, slat, clat)
    ztop = upper[col_sec] + np.where(s_, rv[col_sec], 0.0)

    # rows: points down dip from each column
    pt_col = np.repeat(np.arange(len(c)), rows[col_sec])
    pt_first = np.r_[0, np.cumsum(rows[col_sec])[:-1]]
    r = np.arange(pt_col.size) - pt_first[pt_col]
    ps = col_sec[pt_col]
    h = r * dip_spacing[ps] * np.cos(dip_rad[ps])
    v = r * dip_spacing[ps] * np.sin(dip_rad[ps])
    lon, lat = _nshmp_location(clon[pt_col], clat[pt_col], dip_dir[ps], h)
    lon = np.where(r == 0, clon[pt_col], lon)
    lat = np.where(r == 0, clat[pt_col], lat)
    depth = ztop[pt_col] + v

    # perimeter: trace at the upper depth and the trace moved by the reduced width cos(dip)
    hd = (width * np.sin(dip_rad)) / np.tan(dip_rad)
    blon, blat = _nshmp_location(tlon, tlat, dip_dir[sec_of_pt], hd[sec_of_pt])
    per_lon, per_lat = [], []
    for s in range(n_sec):
        sl = slice(toff[s], toff[s + 1])
        per_lon.append(np.r_[tlon[sl], blon[sl][::-1], tlon[toff[s]]])
        per_lat.append(np.r_[tlat[sl], blat[sl][::-1], tlat[toff[s]]])
    sizes = rows * cols
    return {
        "lon": lon, "lat": lat, "depth": depth,
        "surface_offset": np.r_[0, np.cumsum(sizes)[:-1]],
        "surface_rows": rows, "surface_cols": cols, "surface_dip": dip,
        "surface_spacing": 0.5 * (strike_spacing + dip_spacing),
        "perimeter_lon": np.concatenate(per_lon), "perimeter_lat": np.concatenate(per_lat),
        "perimeter_offset": np.r_[0, np.cumsum(2 * npts + 1)],
    }


def _inside_polygons(plon, plat, poff, site_lon, site_lat):
    """Even-odd test of the site in closed polygons (vertices plon, plat; polygon i is
    poff[i]:poff[i + 1], last vertex equal to the first)."""
    n = len(plon)
    is_edge = np.ones(n, dtype=bool)
    is_edge[poff[1:] - 1] = False
    a = np.flatnonzero(is_edge)
    ay, by = plat[a], plat[a + 1]
    ax, bx = plon[a], plon[a + 1]
    cond = (ay > site_lat) != (by > site_lat)
    with np.errstate(invalid="ignore", divide="ignore"):
        xint = ax + (site_lat - ay) * (bx - ax) / (by - ay)
    cross = (cond & (site_lon < xint)).astype(np.int64)
    first = np.r_[0, np.cumsum(np.diff(poff) - 1)[:-1]]
    return (np.add.reduceat(cross, first) % 2) == 1


def gridded_section_distances(grids, site_lat, site_lon, need_rx=True):
    """rJB, rRup, and rX of fault sections as nshmp-lib computes them (model.Distance.compute)
    from their gridded surfaces (nshmp_section_grids): minimum horizontal
    (Locations.horzDistanceFast) and sqrt(horizontal^2 + depth^2) distances to the grid points
    (top row only if dip > 89 degrees), rJB = 0 if it is less than the average grid spacing and
    the site is inside the perimeter, and rX = extended_trace_rx of the top row. Returns (rjb,
    rrup, rx) per section (rx None if need_rx is False)."""
    lon, lat, depth = grids["lon"], grids["lat"], grids["depth"]
    offset = np.asarray(grids["surface_offset"], dtype=np.int64)
    rows = np.asarray(grids["surface_rows"], dtype=np.int64)
    cols = np.asarray(grids["surface_cols"], dtype=np.int64)
    h = _nshmp_horz_distance_fast(site_lon, site_lat, lon, lat)
    r2 = h * h + depth * depth
    if "_top_only" not in grids:
        sec = np.repeat(np.arange(len(rows)), rows * cols)
        row = (np.arange(len(lon)) - offset[sec]) % rows[sec]
        grids["_top_only"] = (np.asarray(grids["surface_dip"])[sec] > 89.0) & (row > 0)
    top_only = grids["_top_only"]
    h_use = np.where(top_only, np.inf, h)
    rjb = np.minimum.reduceat(h_use, offset)
    rrup = np.sqrt(np.minimum.reduceat(np.where(top_only, np.inf, r2), offset))
    near = rjb < grids["surface_spacing"]
    if np.any(near):
        inside = _inside_polygons(grids["perimeter_lon"], grids["perimeter_lat"],
                                  np.asarray(grids["perimeter_offset"]), site_lon, site_lat)
        rjb = np.where(near & inside, 0.0, rjb)
    rx = None
    if need_rx:
        n_pieces = cols - 1
        sec_of_piece = np.repeat(np.arange(len(cols)), n_pieces)
        first = np.r_[0, np.cumsum(n_pieces)[:-1]]
        col = np.arange(n_pieces.sum()) - first[sec_of_piece]
        i1 = offset[sec_of_piece] + col * rows[sec_of_piece]
        i2 = i1 + rows[sec_of_piece]
        rx = extended_trace_rx(lon[i1], lat[i1], lon[i2], lat[i2], first, site_lon, site_lat)
    return rjb, rrup, rx


# ------------------------------------------------------------------------------------------------
# Site-centred local projection and the GC2 coordinates of Spudich and Chiou (2015)
# ------------------------------------------------------------------------------------------------

LOCAL_EARTH_RADIUS = NSHMP_EARTH_RADIUS  # km, sphere of the local projection
GC2_RUN_GAP = 1.0  # km, consecutive pieces farther apart than this start a new GC2 trace


def spheroid_radius(lat):
    """Radius (km) of the oblate spheroid of point_to_xyz at latitude lat (degrees)."""
    rad = np.pi / 180.0
    a = 6378.1370
    b = 6356.7523
    return np.sqrt(
        ((a**2 * np.cos(lat * rad)) ** 2 + (b**2 * np.sin(lat * rad)) ** 2)
        / ((a * np.cos(lat * rad)) ** 2 + (b * np.sin(lat * rad)) ** 2)
    )


def xyz_to_lat_lon_depth(xyz):
    """Latitude, longitude (degrees), and depth (km) of Cartesian points made by point_to_xyz
    or the source model conversion scripts (spheroid radius minus depth)."""
    lon, lat = xyz_to_lon_lat(xyz)
    depth = spheroid_radius(lat) - np.linalg.norm(xyz, axis=-1)
    return lat, lon, depth


def local_projection(lat, lon, site_lat, site_lon, radius=LOCAL_EARTH_RADIUS):
    """Azimuthal equidistant projection on a sphere centred on the site.

    Returns x (east) and y (north) in km: the great-circle distance from the site times the
    sine and cosine of the azimuth from the site. Distances from the site are surface
    distances, and the vertical is the same at the site and at the fault.
    """
    rad = np.pi / 180.0
    p1 = site_lat * rad
    p2 = np.asarray(lat) * rad
    dl = (np.asarray(lon) - site_lon) * rad
    hav = np.sin((p2 - p1) / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2.0) ** 2
    c = 2.0 * np.arcsin(np.sqrt(np.clip(hav, 0.0, 1.0)))
    az = np.arctan2(np.sin(dl) * np.cos(p2),
                    np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl))
    d = radius * c
    return d * np.sin(az), d * np.cos(az)


def local_triangle_distance(tri_lld, site_lat, site_lon, site_z, tri_segment_id):
    """point_triangle_distance in the site-centred local projection.

    tri_lld: N x 3 x 3 array of triangle vertices as [latitude, longitude, depth (km)]. The
    vertices are projected with local_projection and put at z = -depth; the site is at
    (0, 0, site_z) (site_z = elevation in km for rrup, 0 for rjb). Unlike distances in
    Earth-centred Cartesian coordinates (chords, with a vertical that tilts with the distance
    from the site), these are surface distances, and rrup >= rjb. Returns the minimum distance
    per segment, as point_triangle_distance.
    """
    x, y = local_projection(tri_lld[..., 0], tri_lld[..., 1], site_lat, site_lon)
    tri = np.stack([x, y, -tri_lld[..., 2]], axis=-1)
    return point_triangle_distance(tri, np.array([0.0, 0.0, site_z]), tri_segment_id)


def _first_where(mask, starts, n):
    """Index of the first True of mask in each group (groups start at starts)."""
    big = np.iinfo(np.int64).max
    return np.minimum.reduceat(np.where(mask, np.arange(n), big), starts)


def gc2_coordinates(x1, y1, x2, y2, piece_first, piece_count, entry_segment, group_start,
                    run_gap=GC2_RUN_GAP):
    """GC2 T and U (Spudich and Chiou, 2015, USGS Open-File Report 2015-1028) of a site at the
    origin for groups of fault segments, vectorized over the groups.

    x1, y1, x2, y2: start and end (km, planar coordinates centred on the site) of the upper edge
        of each planar piece, oriented with the fault dipping to the right (Aki and Richards).
    piece_first, piece_count: first piece and number of pieces of each segment (the pieces of a
        segment are consecutive and in order along its upper edge).
    entry_segment: the segments of the groups, group after group, each group in order along its
        rupture (as in ruptures_segments.npz).
    group_start: first entry of each group (increasing, starting with 0).

    For each group: a_hat joins the two most distant piece end points (found by two farthest-
    point sweeps), oriented eastward; pieces whose direction projects on a_hat with the
    opposite sign of the total are reversed (discordance); b_hat, the nominal strike, is the
    direction of the sum of the oriented pieces, and the origin p0 is the end of a_hat that
    comes first along it. Consecutive pieces closer than run_gap (km) form a trace; the offset
    of a piece is the position of the start of its trace along b_hat relative to p0 plus the
    length of the pieces before it in the trace (in the direction of b_hat). With t and u the
    site coordinates relative to a piece (t positive to its right) and
    w = [atan((l - u)/t) - atan(-u/t)]/t its weight, T = sum(w t)/sum(w) and
    U = sum(w (u + offset))/sum(w); a site on a piece gets T = 0 and U = u + offset. U_max, the
    largest U of the piece ends, is used for Ry0 (gc2_ry0).

    Returns (T, U, U_max), one value per group.
    """
    eps = 1e-9
    entry_segment = np.asarray(entry_segment, dtype=np.int64)
    group_start = np.asarray(group_start, dtype=np.int64)
    n_entry = len(entry_segment)
    cnt = piece_count[entry_segment]
    n = int(cnt.sum())
    entry_row0 = np.r_[0, np.cumsum(cnt)[:-1]]
    piece = np.repeat(piece_first[entry_segment] - entry_row0, cnt) + np.arange(n)
    is_start = np.zeros(n_entry, dtype=bool)
    is_start[group_start] = True
    g = np.repeat(np.cumsum(is_start) - 1, cnt)
    gs = entry_row0[group_start]  # first row of each group
    gcount = np.diff(np.r_[gs, n])

    ax, ay, bx, by = x1[piece], y1[piece], x2[piece], y2[piece]
    px1, py1, px2, py2 = ax, ay, bx, by  # as stored, for the trace gaps

    # a_hat: two farthest-point sweeps, starting from the centroid of the end points
    cx = np.add.reduceat(ax + bx, gs) / (2 * gcount)
    cy = np.add.reduceat(ay + by, gs) / (2 * gcount)

    def farthest(px, py):
        da = np.hypot(ax - px[g], ay - py[g])
        db = np.hypot(bx - px[g], by - py[g])
        dm = np.maximum(da, db)
        i = _first_where(dm == np.maximum.reduceat(dm, gs)[g], gs, n)
        use_a = da[i] >= db[i]
        return np.where(use_a, ax[i], bx[i]), np.where(use_a, ay[i], by[i])

    f1x, f1y = farthest(cx, cy)
    f2x, f2y = farthest(f1x, f1y)
    east = f1x <= f2x
    begx, begy = np.where(east, f1x, f2x), np.where(east, f1y, f2y)
    endx, endy = np.where(east, f2x, f1x), np.where(east, f2y, f1y)
    avx, avy = endx - begx, endy - begy
    degenerate = np.hypot(avx, avy) <= eps  # all end points coincide: use the first piece
    avx = np.where(degenerate, (bx - ax)[gs], avx)
    avy = np.where(degenerate, (by - ay)[gs], avy)
    an = np.maximum(np.hypot(avx, avy), eps)
    ahx, ahy = avx / an, avy / an

    # discordance
    e = (bx - ax) * ahx[g] + (by - ay) * ahy[g]
    sgn = np.sign(np.add.reduceat(e, gs))
    sgn[sgn == 0] = 1.0
    flip = np.sign(e) * sgn[g] < 0
    ax, bx = np.where(flip, bx, ax), np.where(flip, ax, bx)
    ay, by = np.where(flip, by, ay), np.where(flip, ay, by)
    dx, dy = bx - ax, by - ay
    length = np.hypot(dx, dy)

    bvx, bvy = np.add.reduceat(dx, gs), np.add.reduceat(dy, gs)
    bn = np.maximum(np.hypot(bvx, bvy), eps)
    bhx, bhy = bvx / bn, bvy / bn
    first_beg = avx * bhx + avy * bhy >= 0
    p0x, p0y = np.where(first_beg, begx, endx), np.where(first_beg, begy, endy)

    # traces: runs of consecutive pieces closer than run_gap
    new_run = np.ones(n, dtype=bool)
    if n > 1:
        gap = np.minimum(
            np.minimum(np.hypot(px1[1:] - px2[:-1], py1[1:] - py2[:-1]),
                       np.hypot(px1[1:] - px1[:-1], py1[1:] - py1[:-1])),
            np.minimum(np.hypot(px2[1:] - px2[:-1], py2[1:] - py2[:-1]),
                       np.hypot(px2[1:] - px1[:-1], py2[1:] - py1[:-1])),
        )
        new_run[1:] = gap > run_gap
    new_run[gs] = True
    rs = np.flatnonzero(new_run)
    run = np.cumsum(new_run) - 1
    re = np.r_[rs[1:], n] - 1
    proj_a = ax * bhx[g] + ay * bhy[g]
    proj_b = bx * bhx[g] + by * bhy[g]
    mid = 0.5 * (proj_a + proj_b)
    forward = mid[re] >= mid[rs]
    csum = np.cumsum(length)
    before = csum - length - (csum[rs] - length[rs])[run]  # length of earlier pieces of the run
    total = np.add.reduceat(length, rs)
    p0b = p0x * bhx + p0y * bhy
    base = np.minimum.reduceat(np.minimum(proj_a, proj_b), rs) - p0b[g[rs]]
    s = base[run] + np.where(forward[run], before, total[run] - before - length)
    u_max = np.maximum.reduceat(s + length, gs)

    # site (origin) coordinates relative to each piece
    ln = np.maximum(length, eps)
    uhx, uhy = dx / ln, dy / ln
    t = -ax * uhy + ay * uhx  # (site - start) x u_hat: positive to the right of the piece
    u = -ax * uhx - ay * uhy
    w = np.zeros(n)
    ok = length > eps
    nz = ok & (np.abs(t) > eps)
    off = ok & ~nz & ((u < -eps) | (u > length + eps))
    on = ok & ~nz & ~off
    with np.errstate(divide="ignore", invalid="ignore"):
        w[nz] = (np.arctan((length[nz] - u[nz]) / t[nz]) - np.arctan(-u[nz] / t[nz])) / t[nz]
        w[off] = 1.0 / (u[off] - length[off]) - 1.0 / u[off]
        sw = np.add.reduceat(w, gs)
        T = np.add.reduceat(w * t, gs) / sw
        U = np.add.reduceat(w * (u + s), gs) / sw
    if on.any():
        hit = np.add.reduceat(on.astype(np.int64), gs) > 0
        i = _first_where(on, gs, n)[hit]
        T[hit] = 0.0
        U[hit] = (u + s)[i]
    return T, U, u_max


def haversine_distance(lat1, lon1, lat2, lon2, radius=LOCAL_EARTH_RADIUS):
    """Great-circle distance (km) between points (degrees)."""
    rad = np.pi / 180.0
    p1, p2 = lat1 * rad, lat2 * rad
    hav = np.sin((p2 - p1) / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin((lon2 - lon1) * rad / 2.0) ** 2
    return 2.0 * radius * np.arcsin(np.sqrt(np.clip(hav, 0.0, 1.0)))


def gc2_rupture_setup(seg_lat, seg_lon, seg_length, entry_segment, rupture_start,
                      run_gap=GC2_RUN_GAP):
    """Site-independent part of GC2 for ruptures made of fault segments (vectorized).

    seg_lat, seg_lon: n_seg x 2 arrays with the start and end of the upper edge of each
        segment (start of its first piece and end of its last piece, with the fault dipping to
        the right of the pieces).
    seg_length: length of the upper edge of each segment (km, sum over its pieces).
    entry_segment: the segments of the ruptures, rupture after rupture, each in order along its
        rupture (ruptures_segments.npz).
    rupture_start: first entry of each rupture.

    The steps are those of gc2_coordinates with each segment as one piece of trace (a segment
    is reversed as a whole if it is discordant, as Spudich and Chiou reverse whole traces), in
    an azimuthal equidistant frame centred on each rupture: a_hat from two farthest-point
    sweeps over the segment end points, oriented eastward; discordance; b_hat and the origin
    p0; traces of segments closer than run_gap (km); and the offset of each entry along b_hat.

    Returns (flip, offset, u_max): whether each entry is reversed, the offset (km) of the start
    of each entry (after reversal) along the GC2 U axis, and U_max of each rupture.
    """
    eps = 1e-9
    es = np.asarray(entry_segment, dtype=np.int64)
    rs0 = np.asarray(rupture_start, dtype=np.int64)
    n = len(es)
    is_start = np.zeros(n, dtype=bool)
    is_start[rs0] = True
    g = np.cumsum(is_start) - 1
    gcount = np.diff(np.r_[rs0, n])

    # azimuthal equidistant frame centred on each rupture (mean of its segment end points)
    lat0 = (np.add.reduceat(seg_lat[es].sum(axis=1), rs0) / (2 * gcount))[g]
    lon0 = (np.add.reduceat(seg_lon[es].sum(axis=1), rs0) / (2 * gcount))[g]
    ax, ay = local_projection(seg_lat[es, 0], seg_lon[es, 0], lat0, lon0)
    bx, by = local_projection(seg_lat[es, 1], seg_lon[es, 1], lat0, lon0)
    gap_ax, gap_ay, gap_bx, gap_by = ax, ay, bx, by
    length = seg_length[es]

    cx = np.add.reduceat(ax + bx, rs0) / (2 * gcount)
    cy = np.add.reduceat(ay + by, rs0) / (2 * gcount)

    def farthest(px, py):
        da = np.hypot(ax - px[g], ay - py[g])
        db = np.hypot(bx - px[g], by - py[g])
        dm = np.maximum(da, db)
        i = _first_where(dm == np.maximum.reduceat(dm, rs0)[g], rs0, n)
        use_a = da[i] >= db[i]
        return np.where(use_a, ax[i], bx[i]), np.where(use_a, ay[i], by[i])

    f1x, f1y = farthest(cx, cy)
    f2x, f2y = farthest(f1x, f1y)
    east = f1x <= f2x
    begx, begy = np.where(east, f1x, f2x), np.where(east, f1y, f2y)
    endx, endy = np.where(east, f2x, f1x), np.where(east, f2y, f1y)
    avx, avy = endx - begx, endy - begy
    degenerate = np.hypot(avx, avy) <= eps
    avx = np.where(degenerate, (bx - ax)[rs0], avx)
    avy = np.where(degenerate, (by - ay)[rs0], avy)
    an = np.maximum(np.hypot(avx, avy), eps)
    ahx, ahy = avx / an, avy / an

    e = (bx - ax) * ahx[g] + (by - ay) * ahy[g]
    sgn = np.sign(np.add.reduceat(e, rs0))
    sgn[sgn == 0] = 1.0
    flip = np.sign(e) * sgn[g] < 0
    ax, bx = np.where(flip, bx, ax), np.where(flip, ax, bx)
    ay, by = np.where(flip, by, ay), np.where(flip, ay, by)
    dx, dy = bx - ax, by - ay
    bvx, bvy = np.add.reduceat(dx, rs0), np.add.reduceat(dy, rs0)
    bn = np.maximum(np.hypot(bvx, bvy), eps)
    bhx, bhy = bvx / bn, bvy / bn
    first_beg = avx * bhx + avy * bhy >= 0
    p0x, p0y = np.where(first_beg, begx, endx), np.where(first_beg, begy, endy)

    new_run = np.ones(n, dtype=bool)
    if n > 1:
        gap = np.minimum(
            np.minimum(np.hypot(gap_ax[1:] - gap_bx[:-1], gap_ay[1:] - gap_by[:-1]),
                       np.hypot(gap_ax[1:] - gap_ax[:-1], gap_ay[1:] - gap_ay[:-1])),
            np.minimum(np.hypot(gap_bx[1:] - gap_bx[:-1], gap_by[1:] - gap_by[:-1]),
                       np.hypot(gap_bx[1:] - gap_ax[:-1], gap_by[1:] - gap_ay[:-1])),
        )
        new_run[1:] = gap > run_gap
    new_run[rs0] = True
    rs = np.flatnonzero(new_run)
    run = np.cumsum(new_run) - 1
    re = np.r_[rs[1:], n] - 1
    proj_a = ax * bhx[g] + ay * bhy[g]
    proj_b = bx * bhx[g] + by * bhy[g]
    mid = 0.5 * (proj_a + proj_b)
    forward = mid[re] >= mid[rs]
    csum = np.cumsum(length)
    before = csum - length - (csum[rs] - length[rs])[run]
    total = np.add.reduceat(length, rs)
    p0b = p0x * bhx + p0y * bhy
    base = np.minimum.reduceat(np.minimum(proj_a, proj_b), rs) - p0b[g[rs]]
    offset = base[run] + np.where(forward[run], before, total[run] - before - length)
    u_max = np.maximum.reduceat(offset + length, rs0)
    return flip, offset, u_max


def gc2_segment_sums(x1, y1, x2, y2, piece_first, piece_before):
    """Site-dependent part of GC2: the weights of the pieces for a site at the origin, summed
    over the pieces of each segment.

    x1, y1, x2, y2: upper edges of the pieces (km, centred on the site), oriented with the
        fault dipping to their right; piece_first: first piece of each segment (the pieces of a
        segment are consecutive and in order); piece_before: length (km) of the pieces of the
        segment before each piece.

    Returns (W, A, B, on, on_u): per segment W = sum(w), A = sum(w t), and B = sum(w (u + c))
    with c = piece_before, as in gc2_coordinates; on is True for segments with the site on one
    of their pieces (T = 0 there) and on_u the u + c of that piece.
    """
    eps = 1e-9
    dx, dy = x2 - x1, y2 - y1
    length = np.hypot(dx, dy)
    ln = np.maximum(length, eps)
    uhx, uhy = dx / ln, dy / ln
    t = -x1 * uhy + y1 * uhx
    u = -x1 * uhx - y1 * uhy
    w = np.zeros(len(x1))
    ok = length > eps
    nz = ok & (np.abs(t) > eps)
    off = ok & ~nz & ((u < -eps) | (u > length + eps))
    on = ok & ~nz & ~off
    with np.errstate(divide="ignore", invalid="ignore"):
        w[nz] = (np.arctan((length[nz] - u[nz]) / t[nz]) - np.arctan(-u[nz] / t[nz])) / t[nz]
        w[off] = 1.0 / (u[off] - length[off]) - 1.0 / u[off]
    uc = u + piece_before
    W = np.add.reduceat(w, piece_first)
    A = np.add.reduceat(w * t, piece_first)
    B = np.add.reduceat(w * uc, piece_first)
    on_seg = np.add.reduceat(on.astype(np.int64), piece_first) > 0
    on_u = np.zeros(len(piece_first))
    if on.any():
        i = _first_where(on, piece_first, len(x1))
        on_u[on_seg] = uc[i[on_seg]]
    return W, A, B, on_seg, on_u


def gc2_rupture_coordinates(W, A, B, on_seg, on_u, seg_length, entry_segment, rupture_start,
                            flip, offset):
    """GC2 T and U of ruptures from the segment sums of gc2_segment_sums and the set-up of
    gc2_rupture_setup (entries and ruptures as there; flip and offset per entry).

    A reversed segment contributes -A to sum(w t) and (offset + L) W - B to sum(w (u + s));
    one in its own direction contributes A and B + offset W. A site on a segment of the
    rupture gets T = 0 and U = u + offset of that segment.
    """
    es = np.asarray(entry_segment, dtype=np.int64)
    rs0 = np.asarray(rupture_start, dtype=np.int64)
    w = W[es]
    a = A[es]
    b = B[es]
    L = seg_length[es]
    sw = np.add.reduceat(w, rs0)
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.add.reduceat(np.where(flip, -a, a), rs0) / sw
        U = np.add.reduceat(np.where(flip, (offset + L) * w - b, b + offset * w), rs0) / sw
    on = on_seg[es]
    if on.any():
        hit = np.add.reduceat(on.astype(np.int64), rs0) > 0
        i = _first_where(on, rs0, len(es))[hit]
        T[hit] = 0.0
        ou = on_u[es[i]]
        U[hit] = np.where(flip[i], offset[i] + L[i] - ou, offset[i] + ou)
    return T, U


def gc2_ry0(U, u_max):
    """Ry0 from the GC2 U coordinate: -U before the start, U - U_max after the end, else 0."""
    return np.where(U < 0.0, -U, np.where(U > u_max, U - u_max, 0.0))

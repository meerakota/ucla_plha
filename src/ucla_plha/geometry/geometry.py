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


def xyz_to_lon_lat(xyz):
    """Longitude and latitude (degrees) of Cartesian points (inverse of the angles of
    point_to_xyz)."""
    lon = np.degrees(np.arctan2(xyz[..., 1], xyz[..., 0]))
    lat = np.degrees(np.arctan2(xyz[..., 2], np.hypot(xyz[..., 0], xyz[..., 1])))
    return lon, lat


# ------------------------------------------------------------------------------------------------
# Site-centred local projection and the GC2 coordinates of Spudich and Chiou (2015)
# ------------------------------------------------------------------------------------------------

LOCAL_EARTH_RADIUS = 6371.0072  # km, mean Earth radius, sphere of the local projection
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

"""Tests of the site-centred local projection (Rjb, Rrup) and of the GC2 coordinates (Rx, Ry0)
of Spudich and Chiou (2015)."""

import numpy as np
import pytest

from ucla_plha import plha
from ucla_plha.geometry import geometry


def _one_rupture(starts, ends):
    """gc2_coordinates for one rupture whose pieces are each a segment (site at the origin)."""
    starts, ends = np.asarray(starts, float), np.asarray(ends, float)
    n = len(starts)
    T, U, u_max = geometry.gc2_coordinates(
        starts[:, 0], starts[:, 1], ends[:, 0], ends[:, 1],
        np.arange(n), np.ones(n, dtype=np.int64), np.arange(n), np.array([0]),
    )
    return T[0], U[0], u_max[0]


def _shift(points, site):
    return np.asarray(points, float) - np.asarray(site, float)


# ------------------------------------------------------------------------------------------------
# GC2 on planar traces
# ------------------------------------------------------------------------------------------------

@pytest.mark.parametrize(
    "site,T,U,ry0",
    [
        ((5.0, -3.0), 3.0, 5.0, 0.0),  # right of the trace (hanging wall), above the middle
        ((5.0, 4.0), -4.0, 5.0, 0.0),  # left (footwall)
        ((-6.0, 2.0), None, None, 6.0),  # before the start
        ((13.0, -1.0), None, None, 3.0),  # beyond the end
        ((4.0, 0.0), 0.0, 4.0, 0.0),  # on the trace
    ],
)
def test_gc2_single_piece(site, T, U, ry0):
    # trace from (0, 0) to (10, 0): strike east, fault dipping south (to the right)
    t, u, u_max = _one_rupture(_shift([[0, 0]], site), _shift([[10, 0]], site))
    assert u_max == pytest.approx(10.0)
    if T is not None:
        assert t == pytest.approx(T, abs=1e-9)
        assert u == pytest.approx(U, abs=1e-9)
    assert geometry.gc2_ry0(np.array([u]), np.array([u_max]))[0] == pytest.approx(ry0, abs=1e-9)


def test_gc2_collinear_split_and_reversed_pieces():
    # a straight trace split into three pieces, one of them reversed (discordant), gives the
    # coordinates of the single piece
    rng = np.random.default_rng(3)
    for site in rng.uniform(-30, 40, size=(20, 2)):
        ref = _one_rupture(_shift([[0, 0]], site), _shift([[10, 0]], site))
        split = _one_rupture(_shift([[0, 0], [7, 0], [3, 0]], site),
                             _shift([[3, 0], [10, 0], [7, 0]], site))
        np.testing.assert_allclose(split, ref, atol=1e-9)


def test_gc2_kinked_trace_matches_spudich_chiou_weights():
    # two pieces at a right angle: T and U are the w-weighted averages of the local t and u
    site = np.array([12.0, 5.0])
    a, b, c = np.array([0.0, 0.0]), np.array([10.0, 0.0]), np.array([10.0, -10.0])
    t, u, u_max = _one_rupture(_shift([a, b], site), _shift([b, c], site))

    def local(p, q, s):
        uh = (q - p) / np.linalg.norm(q - p)
        th = np.array([uh[1], -uh[0]])
        tt, uu, ll = (site - p) @ th, (site - p) @ uh, np.linalg.norm(q - p)
        w = (np.arctan((ll - uu) / tt) - np.arctan(-uu / tt)) / tt
        return w, tt, uu + s

    w1, t1, u1 = local(a, b, 0.0)
    w2, t2, u2 = local(b, c, 10.0)
    assert t == pytest.approx((w1 * t1 + w2 * t2) / (w1 + w2), abs=1e-9)
    assert u == pytest.approx((w1 * u1 + w2 * u2) / (w1 + w2), abs=1e-9)
    assert u_max == pytest.approx(20.0)


def test_gc2_segment_form_matches_piece_form():
    # gc2_rupture_setup + gc2_segment_sums + gc2_rupture_coordinates (the cached form used by
    # get_source_data) against gc2_coordinates, for random concordant kinked ruptures of several
    # multi-piece segments
    rng = np.random.default_rng(7)
    lat0, lon0 = 36.0, -120.0
    km = geometry.LOCAL_EARTH_RADIUS * np.pi / 180.0
    for _ in range(25):
        n_seg = rng.integers(1, 5)
        counts = rng.integers(1, 4, n_seg)
        n = counts.sum()
        az = np.radians(rng.uniform(100, 140)) + np.cumsum(rng.normal(0, 0.15, n))
        step = rng.uniform(2, 8, n)
        xy = np.vstack([[0, 0], np.cumsum(np.c_[step * np.sin(az), step * np.cos(az)], 0)])
        lat = lat0 + xy[:, 1] / km
        lon = lon0 + xy[:, 0] / (km * np.cos(np.radians(lat0)))
        top_lat, top_lon = np.c_[lat[:-1], lat[1:]], np.c_[lon[:-1], lon[1:]]
        first = np.r_[0, np.cumsum(counts)[:-1]]
        plen = geometry.haversine_distance(top_lat[:, 0], top_lon[:, 0], top_lat[:, 1], top_lon[:, 1])
        seg_len = np.add.reduceat(plen, first)
        seg = np.repeat(np.arange(n_seg), counts)
        before = np.cumsum(plen) - plen - (np.cumsum(plen)[first] - plen[first])[seg]
        last = first + counts - 1
        flip, offset, u_max = geometry.gc2_rupture_setup(
            np.c_[top_lat[first, 0], top_lat[last, 1]], np.c_[top_lon[first, 0], top_lon[last, 1]],
            seg_len, np.arange(n_seg), np.array([0]))
        for slat, slon in np.c_[lat0 + rng.uniform(-0.4, 0.4, 5), lon0 + rng.uniform(-0.4, 0.4, 5)]:
            x1, y1 = geometry.local_projection(top_lat[:, 0], top_lon[:, 0], slat, slon)
            x2, y2 = geometry.local_projection(top_lat[:, 1], top_lon[:, 1], slat, slon)
            sums = geometry.gc2_segment_sums(x1, y1, x2, y2, first, before)
            T, U = geometry.gc2_rupture_coordinates(*sums, seg_len, np.arange(n_seg), np.array([0]),
                                                    flip, offset)
            T_ref, U_ref, um_ref = geometry.gc2_coordinates(
                x1, y1, x2, y2, first, counts, np.arange(n_seg), np.array([0]))
            assert T[0] == pytest.approx(T_ref[0], abs=1e-6)
            # U: rupture-centred vs site-centred frames differ by well under 0.1%
            assert U[0] == pytest.approx(U_ref[0], abs=0.05)
            assert u_max[0] == pytest.approx(um_ref[0], abs=0.05)


# ------------------------------------------------------------------------------------------------
# Local projection
# ------------------------------------------------------------------------------------------------

def test_local_projection_distance_and_azimuth():
    x, y = geometry.local_projection(np.array([35.0]), np.array([-118.0]), 34.0, -118.0)
    km = geometry.LOCAL_EARTH_RADIUS * np.pi / 180.0
    assert x[0] == pytest.approx(0.0, abs=1e-9)
    assert y[0] == pytest.approx(km, rel=1e-12)  # due north, great-circle distance
    lat, lon = np.array([36.3]), np.array([-115.1])
    x, y = geometry.local_projection(lat, lon, 34.0, -118.0)
    d = geometry.haversine_distance(34.0, -118.0, lat, lon)
    assert np.hypot(x, y)[0] == pytest.approx(d[0], rel=1e-12)


def test_local_triangle_distance_vertical_and_horizontal():
    # a triangle 5 km below the site and a triangle 100 km north at 10 km depth
    km = geometry.LOCAL_EARTH_RADIUS * np.pi / 180.0
    d = 0.01
    below = [[34.0 - d, -118.0 - d, 5.0], [34.0 - d, -118.0 + d, 5.0], [34.0 + d, -118.0, 5.0]]
    far_lat = 34.0 + 100.0 / km
    far = [[far_lat, -118.1, 10.0], [far_lat, -117.9, 10.0], [far_lat + 0.1, -118.0, 10.0]]
    tri = np.array([below, far])
    rrup = geometry.local_triangle_distance(tri, 34.0, -118.0, 0.0, np.array([0, 1]))
    surf = tri.copy()
    surf[..., 2] = 0.0
    rjb = geometry.local_triangle_distance(surf, 34.0, -118.0, 0.0, np.array([0, 1]))
    assert rrup[0] == pytest.approx(5.0, abs=1e-6)
    assert rjb[0] == pytest.approx(0.0, abs=1e-6)
    # the closest edge follows a parallel of latitude, which bends 4.5 m toward the pole over
    # 0.1 degree of longitude in the site-centred projection
    assert rjb[1] == pytest.approx(100.0, abs=0.01)
    assert rrup[1] == pytest.approx(np.hypot(100.0, 10.0), abs=0.01)


def test_xyz_to_lat_lon_depth_inverts_point_to_xyz():
    for lat, lon, elev in [(34.0, -118.0, 0.0), (47.6, -122.3, -12.5), (-33.4, 151.2, 1.2)]:
        out = geometry.xyz_to_lat_lon_depth(geometry.point_to_xyz(np.array([lat, lon, elev])))
        np.testing.assert_allclose(out, [lat, lon, -elev], atol=1e-9)


@pytest.mark.parametrize("model", ["ucerf3_fm31", "nshm23_wus"])
def test_rrup_not_less_than_rjb(model):
    # Los Angeles: with the local projection rrup >= rjb for every rupture; the Earth-centred
    # Cartesian triangles give rrup < rjb for most ruptures beyond ~100 km
    p_xyz = geometry.point_to_xyz(np.array([34.05, -118.25, 0.0]))
    # (no distance cutoff, so that both frames keep the same ruptures)
    kw = dict(fault_distances="triangles")
    local = plha.get_source_data("fault_source_models", model, p_xyz, None, 7.0, ["cy14"], **kw)
    geo = plha.get_source_data("fault_source_models", model, p_xyz, None, 7.0, ["cy14"],
                               triangle_frame="geocentric", **kw)
    rjb, rrup = local[3], local[4]
    assert np.all(rrup >= rjb - 1e-9)
    far = (geo[3] > 200.0) & (geo[3] < 500.0)
    assert np.mean(geo[4][far] < geo[3][far]) > 0.5
    # the frames differ by the chord-arc difference and the tilt of the vertical: under 0.1%
    near = geo[3] < 500.0
    np.testing.assert_allclose(local[3][near], geo[3][near], rtol=1e-3, atol=0.05)
    np.testing.assert_allclose(local[4][near], geo[4][near], rtol=2e-3, atol=0.05)


# ------------------------------------------------------------------------------------------------
# NGA flatfile: South Napa (nc72282711), single segment
# ------------------------------------------------------------------------------------------------

# finite fault (nga-gm.org finiteFaultSegments) and flatfile distances (km)
NAPA = dict(lat=38.32726, lon=-122.35902, depth=0.0973, length=15.0, width=10.0,
            strike=155.0, dip=82.0)
NAPA_STATIONS = [
    # station, latitude, longitude, rjb, rrup, rx, ry0
    ("N019B", 38.30095, -122.40154, 3.19981, 4.56040, 4.59158, 0.0000),
    ("N016", 38.29875, -122.28484, 4.53752, 4.53855, -4.53752, 0.0000),
    ("1759", 38.10767, -122.25623, 10.93730, 11.11540, 2.11791, 10.9131),
    ("NMI", 38.07601, -122.25850, 14.22980, 14.52240, 3.78085, 14.0278),
    ("C012", 37.92636, -122.31248, 30.38430, 30.97120, 15.06720, 27.1328),
    ("N006", 38.40629, -122.70915, 30.72850, 31.59900, 23.94660, 20.8692),
    ("GASB", 39.65471, -122.71595, 150.79500, 150.79500, -34.44870, 146.8070),
    ("SCK", 37.18321, -121.07116, 156.01000, 156.01000, -50.25280, 147.6950),
]


def _napa_geometry():
    f = NAPA
    stk = np.radians(f["strike"])
    h = f["width"] * np.cos(np.radians(f["dip"]))
    lon2, lat2 = geometry._nshmp_location(f["lon"], f["lat"], stk, f["length"])
    lon3, lat3 = geometry._nshmp_location(f["lon"], f["lat"], stk + np.pi / 2, h)
    lon4, lat4 = geometry._nshmp_location(lon2, lat2, stk + np.pi / 2, h)
    zb = f["depth"] + f["width"] * np.sin(np.radians(f["dip"]))
    p1, p2 = [f["lat"], f["lon"], f["depth"]], [lat2, lon2, f["depth"]]
    p3, p4 = [lat3, lon3, zb], [lat4, lon4, zb]
    return np.array([[p1, p2, p4], [p1, p3, p4]]), np.array([[f["lat"], lat2]]), np.array([[f["lon"], lon2]])


@pytest.mark.parametrize("station", NAPA_STATIONS, ids=[s[0] for s in NAPA_STATIONS])
def test_south_napa_matches_nga_flatfile(station):
    _, slat, slon, rjb, rrup, rx, ry0 = station
    tri, top_lat, top_lon = _napa_geometry()
    surf = tri.copy()
    surf[..., 2] = 0.0
    seg = np.array([0, 0])
    assert geometry.local_triangle_distance(surf, slat, slon, 0.0, seg)[0] == pytest.approx(rjb, abs=0.05)
    assert geometry.local_triangle_distance(tri, slat, slon, 0.0, seg)[0] == pytest.approx(rrup, abs=0.05)
    length = geometry.haversine_distance(top_lat[:, 0], top_lon[:, 0], top_lat[:, 1], top_lon[:, 1])
    flip, offset, u_max = geometry.gc2_rupture_setup(top_lat, top_lon, length, np.array([0]), np.array([0]))
    x1, y1 = geometry.local_projection(top_lat[:, 0], top_lon[:, 0], slat, slon)
    x2, y2 = geometry.local_projection(top_lat[:, 1], top_lon[:, 1], slat, slon)
    sums = geometry.gc2_segment_sums(x1, y1, x2, y2, np.array([0]), np.array([0.0]))
    T, U = geometry.gc2_rupture_coordinates(*sums, length, np.array([0]), np.array([0]), flip, offset)
    # Rx and Ry0 to 0.1 km near the fault, 0.5% of the distance farther away
    tol = max(0.1, 0.005 * rjb)
    assert T[0] == pytest.approx(rx, abs=tol)
    assert geometry.gc2_ry0(U, u_max)[0] == pytest.approx(ry0, abs=tol)


def test_get_source_data_gc2_default_and_options():
    p_xyz = geometry.point_to_xyz(np.array([34.05, -118.25, 0.0]))
    default = plha.get_source_data("fault_source_models", "ucerf3_fm31", p_xyz, 100.0, 6.0, ["ask14"])
    gc2 = plha.get_source_data("fault_source_models", "ucerf3_fm31", p_xyz, 100.0, 6.0, ["ask14"],
                               rupture_rx="gc2")
    closest = plha.get_source_data("fault_source_models", "ucerf3_fm31", p_xyz, 100.0, 6.0, ["ask14"],
                                   rupture_rx="closest_section")
    for a, b in zip(default, gc2):
        np.testing.assert_array_equal(a, b)
    rx, rx1, ry0 = gc2[5], gc2[6], gc2[7]
    assert np.all(np.isfinite(rx)) and np.all(ry0 >= 0.0)
    # rx1 = rx - the horizontal width of the closest segment, as for "closest_section"
    np.testing.assert_allclose(rx - rx1, closest[5] - closest[6], atol=1e-9)
    with pytest.raises(ValueError):
        plha.get_source_data("fault_source_models", "ucerf3_fm31", p_xyz, 100.0, 6.0, ["ask14"],
                             rupture_rx="gc3")
    with pytest.raises(ValueError):
        plha.get_source_data("fault_source_models", "ucerf3_fm31", p_xyz, 100.0, 6.0, ["ask14"],
                             triangle_frame="flat")


@pytest.mark.parametrize("model", ["ucerf3_fm31", "nshm23_wus"])
def test_gc2_matches_segment_rx_alongside_dipping_segments(model):
    # For a dipping single-piece segment with the site alongside it (Ry0 = 0), GC2 T is the
    # signed distance from the line of the upper edge, i.e. the Rx of get_Rx_Rx1_Ry0
    lat, lon = 34.05, -118.25
    p_xyz = geometry.point_to_xyz(np.array([lat, lon, 0.0]))
    path = plha._source_model_path("fault_source_models", model)
    rect = np.load(str(path.joinpath("rect_rjb.npy")))
    rid = np.load(str(path.joinpath("rect_segment_id.npy")))
    rx_seg, rx1_seg, _ = geometry.get_Rx_Rx1_Ry0(rect, p_xyz, rid)
    local = plha._local_geometry(path)
    n_seg = len(local["piece_first"])
    x1, y1 = geometry.local_projection(local["top_lat"][:, 0], local["top_lon"][:, 0], lat, lon)
    x2, y2 = geometry.local_projection(local["top_lat"][:, 1], local["top_lon"][:, 1], lat, lon)
    T, U, u_max = geometry.gc2_coordinates(x1, y1, x2, y2, local["piece_first"],
                                           local["piece_count"], np.arange(n_seg), np.arange(n_seg))
    ry0 = geometry.gc2_ry0(U, u_max)
    m = ((rx_seg - rx1_seg) > 0.5) & (ry0 == 0.0) & (np.abs(rx_seg) < 50.0) & (local["piece_count"] == 1)
    assert m.sum() >= 2
    np.testing.assert_allclose(T[m], rx_seg[m], atol=0.01)


def test_gc2_ry0_of_vertical_segments():
    # get_Rx_Rx1_Ry0 gives Ry0 = 0 for vertical rectangles wherever the site is (the
    # down-dip edge has no length); GC2 Ry0 is the distance beyond the ends along strike
    rect = np.array([[[0, 0, 0], [10, 0, 0], [0, 0, 0], [10, 0, 0]]], dtype=float)
    _, _, ry0 = geometry.get_Rx_Rx1_Ry0(rect, np.array([25.0, 3.0, 0.0]), np.array([0]))
    assert ry0[0] == 0.0
    t, u, u_max = _one_rupture(_shift([[0, 0]], (25.0, 3.0)), _shift([[10, 0]], (25.0, 3.0)))
    assert geometry.gc2_ry0(np.array([u]), np.array([u_max]))[0] == pytest.approx(15.0)

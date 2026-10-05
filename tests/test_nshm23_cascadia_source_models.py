"""
Tests for the Cascadia subduction source models converted from nshm-conus 6.2.0
(utilities/convert_nshm23_cascadia_interface.py and utilities/convert_nshm23_cascadia_slab.py).
Expected rates are computed by hand from the NSHM logic tree weights and MFD parameters.
"""

import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import pytest

from ucla_plha import plha
from ucla_plha.geometry import geometry

SOURCE_MODELS = files("ucla_plha").joinpath("source_models")
UTILITIES = Path(__file__).parents[1] / "utilities"
INTERFACE = "nshm23_cascadia_interface"
CLUSTER = "nshm23_cascadia_interface_cluster"
SLAB = "nshm23_cascadia_slab"

# Geometry (depth) branch weights in Cascadia/source-tree.json
DEPTH_WEIGHTS = {"bottom": 0.3, "middle": 0.5, "top": 0.2}
# Scale of the unsegmented GEA12-A branch (source-group.json) on each depth branch
GEA12_A_SCALE = {"bottom": 1.856, "middle": 1.8528, "top": 1.856}
# Rates (mfd-map.json) of the segmented ruptures; GEA12 northern has scale 0.25
GEA12_RATES = [0.0002174, 0.0004888, 0.0005438, 0.25 * 0.001]
GEA17_RATES = [0.000156, 0.000364, 0.000104, 0.000468, 0.000208, 0.0001]
FULL_RATE = 0.0019


def _load(source_type, model, name):
    return np.load(str(SOURCE_MODELS.joinpath(source_type, model, name)))


def _info(source_type, model):
    return json.loads(SOURCE_MODELS.joinpath(source_type, model, "source_info.json").read_text())


def _gr_total(a, b):
    # GR mMin 7.95, mMax 8.75, dm 0.1: bin centers 8.0 ... 8.7
    m = 8.0 + 0.1 * np.arange(8)
    return np.sum(10.0 ** (a - b * m))


@pytest.fixture(scope="module")
def interface():
    return dict(_load("fault_source_models", INTERFACE, "ruptures.npz"))


@pytest.fixture(scope="module")
def cluster():
    return dict(_load("fault_source_models", CLUSTER, "ruptures.npz"))


@pytest.fixture(scope="module")
def slab():
    return dict(_load("point_source_models", SLAB, "ruptures.npz"))


@pytest.mark.parametrize(
    "source_type, model, region, component, is_cluster",
    [
        ("fault_source_models", INTERFACE, "subduction_interface", "Interface", False),
        ("fault_source_models", CLUSTER, "subduction_interface", "FaultCluster", True),
        ("point_source_models", SLAB, "subduction_slab", "Slab", False),
    ],
)
def test_source_info(source_type, model, region, component, is_cluster):
    info = _info(source_type, model)
    assert info["name"] == model
    assert info["tectonic_region"] == region
    assert info["nshm_component"] == component
    assert info["cluster"] is is_cluster
    assert "nshm-conus 6.2.0" in info["source"]
    assert "6.2.0" in info["notes"]
    if source_type == "point_source_models":
        assert info["point_source"]["point_source_type"] == "FINITE"
    else:
        assert info["interface_events"] is True


def test_interface_total_rate(interface):
    # Full rupture, cluster-out branch (weight 0.9); the geometry weights sum to 1
    full = 0.9 * FULL_RATE
    # Partial rupture, segmented branch (0.5) with GEA12 and GEA17 branches (0.5 each)
    segmented = 0.5 * (0.5 * sum(GEA12_RATES) + 0.5 * sum(GEA17_RATES))
    # Partial rupture, unsegmented branch (0.5): GR b=0 and b=1 branches (0.5 each), GEA12-B
    # (weight 0.75, scale 1.2) and GEA12-A (weight 0.25, scale per depth branch)
    gr = 0.5 * _gr_total(-3.806, 0.0) + 0.5 * _gr_total(4.485, 1.0)
    scale_a = sum(DEPTH_WEIGHTS[d] * GEA12_A_SCALE[d] for d in DEPTH_WEIGHTS)
    unsegmented = 0.5 * gr * (0.75 * 1.2 + 0.25 * scale_a)
    expected = full + segmented + unsegmented
    assert interface["rate"].sum() == pytest.approx(expected, rel=1e-6)
    assert expected == pytest.approx(0.0032874, rel=1e-4)


def test_interface_rate_of_m9_events(interface):
    # M >= 9 SINGLE magnitudes in mfd-map.json, each with magnitude branch weight 0.25:
    # full rupture bottom 9.01, 9.34, 9.21; middle 9.12, 9.03; top 9.01 (cluster-out weight 0.9),
    # segmented (0.5 x 0.5) bottom GEA17-B 9.20, 9.10, GEA12-B 9.07, GEA17-C 9.02.
    full = 0.9 * FULL_RATE * 0.25 * (3 * 0.3 + 2 * 0.5 + 1 * 0.2)
    seg = 0.3 * 0.25 * 0.25 * (2 * 0.000156 + 0.0002174 + 0.000364)
    rate = interface["rate"][interface["m"] >= 9.0].sum()
    assert rate == pytest.approx(full + seg, rel=1e-6)
    # The largest event is the bottom-geometry full rupture with the M3 (Papazachos) magnitude
    assert interface["m"].max() == pytest.approx(9.34)
    assert interface["rate"][interface["m"] == interface["m"].max()].sum() == pytest.approx(
        0.9 * FULL_RATE * 0.25 * 0.3
    )


def test_interface_rate_of_m8_events(interface):
    # All ruptures are M >= 8 except the GEA17-E segmented ruptures (M 7.52 to 8.14)
    e_mags = {
        "bottom": [7.94, 7.93, 8.09, 8.14],
        "middle": [7.68, 7.63, 7.73, 7.83],
        "top": [7.59, 7.52, 7.61, 7.73],
    }
    below = sum(
        DEPTH_WEIGHTS[d] * 0.25 * 0.25 * 0.000208 * sum(m < 8.0 for m in e_mags[d])
        for d in e_mags
    )
    rate = interface["rate"][interface["m"] >= 8.0 - 1e-9].sum()
    assert rate == pytest.approx(interface["rate"].sum() - below, rel=1e-6)


def test_interface_rupture_properties(interface):
    assert set(np.unique(interface["fault_type"])) == {1}
    np.testing.assert_allclose(interface["rake"], 90.0)
    # Rupture widths are multiples of the 5 km grid spacing, and zbor follows from the width
    np.testing.assert_allclose(interface["width"] / 5.0, np.round(interface["width"] / 5.0))
    np.testing.assert_allclose(
        interface["zbor"],
        interface["ztor"] + interface["width"] * np.sin(np.radians(interface["dip"])),
    )
    assert np.all((interface["dip"] > 5.0) & (interface["dip"] < 15.0))
    assert np.all((interface["ztor"] >= 5.0) & (interface["ztor"] < 7.0))
    assert np.all(interface["rate"] > 0.0)


def test_interface_ruptures_use_consecutive_segments():
    rs = dict(_load("fault_source_models", INTERFACE, "ruptures_segments.npz"))
    r = dict(_load("fault_source_models", INTERFACE, "ruptures.npz"))
    breaks = np.diff(rs["rupture_index"]) != 0
    assert np.all(np.diff(rs["segment_index"])[~breaks] == 1)
    # Floating ruptures of the unsegmented branches are shorter than the full rupture
    n_seg = np.bincount(rs["rupture_index"])
    assert len(n_seg) == len(r["m"])
    full = n_seg[np.isclose(r["m"], 9.34)]
    assert len(full) == 1 and full[0] == n_seg.max()
    # Full rupture length: 10^((M - 4.94) / 1.39) km at 5 km per segment for M 8.0 floaters
    m8 = n_seg[np.isclose(r["m"], 8.0)]
    assert np.all(np.abs(m8 * 5.0 - 10.0 ** ((8.0 - 4.94) / 1.39)) < 5.0)


def test_interface_full_rupture_spans_the_margin():
    r = dict(_load("fault_source_models", INTERFACE, "ruptures.npz"))
    rs = dict(_load("fault_source_models", INTERFACE, "ruptures_segments.npz"))
    tri = _load("fault_source_models", INTERFACE, "tri_rrup.npy")
    tri_id = _load("fault_source_models", INTERFACE, "tri_segment_id.npy")
    i = np.flatnonzero(np.isclose(r["m"], 9.34))[0]
    seg = rs["segment_index"][rs["rupture_index"] == i]
    xyz = tri[np.isin(tri_id, seg)].reshape(-1, 3)
    lat = np.degrees(np.arcsin(xyz[:, 2] / np.linalg.norm(xyz, axis=1)))
    assert lat.min() == pytest.approx(40.35, abs=0.05)
    assert lat.max() > 49.5


@pytest.mark.parametrize(
    "site, rrup_min, rrup_max",
    [
        ((47.6, -122.3), 85.0, 95.0),  # Seattle
        ((45.52, -122.68), 72.0, 85.0),  # Portland
        ((44.63, -124.05), 20.0, 30.0),  # Newport, Oregon
        ((40.8, -124.16), 14.0, 22.0),  # Eureka, California
    ],
)
def test_interface_distances(site, rrup_min, rrup_max):
    p_xyz = geometry.point_to_xyz(np.array([site[0], site[1], 0.0]))
    out = plha.get_source_data(
        "fault_source_models", INTERFACE, p_xyz, None, None, ["ask14", "bssa14"]
    )
    rrup = out[4]
    assert rrup_min < rrup.min() < rrup_max
    assert np.all(np.isfinite(rrup))


def test_cluster_model(cluster):
    c = dict(_load("fault_source_models", CLUSTER, "clusters.npz"))
    # 3 geometry branches x 3 cluster sets (7a, 7b, 8) with cluster rate 0.0019 and total weight 0.1
    assert len(c["cluster_id"]) == 9
    np.testing.assert_allclose(c["rate"], 0.0019)
    assert c["weight"].sum() == pytest.approx(0.1, abs=1e-8)
    assert set(np.unique(cluster["cluster_id"])) == set(c["cluster_id"])
    # Magnitude branch weights of each section sum to one
    weights = np.bincount(cluster["cluster_section"], weights=cluster["rate"])
    np.testing.assert_allclose(weights, 1.0)
    sections = [
        len(np.unique(cluster["cluster_section"][cluster["cluster_id"] == i]))
        for i in c["cluster_id"]
    ]
    assert sorted(sections) == [7] * 6 + [8] * 3
    assert np.all((cluster["m"] > 7.5) & (cluster["m"] < 8.6))
    assert set(np.unique(cluster["fault_type"])) == {1}


def test_slab_total_rates(slab):
    # rate(M) = R * pdf * 10^(-b M) on GR bin centers; pdfs sum to one in each state
    def gr(R, b, m_min, m_max):
        m = m_min + 0.05 + 0.1 * np.arange(int(round((m_max - m_min) / 0.1)))
        return R * 10.0 ** (-b * m), m

    def total(branches, m_lo):
        return sum(w * np.sum(r[m >= m_lo]) for w, (r, m) in branches)

    hi = lambda R: [(0.9, gr(R, 0.8, 7.2, 7.5)), (0.1, gr(R, 0.8, 7.2, 8.0))]
    expected = {
        "WA": hi(90.38385106241607) + [(1.0, gr(1.20769448734618, 0.4, 5.0, 7.2))],
        "OR": hi(10.801015830699988) + [(1.0, gr(0.1443214369750003, 0.4, 6.5, 7.2))],
        "CA": [(0.9, gr(51.28063780017224, 0.8, 5.0, 7.5)), (0.1, gr(51.28063780017224, 0.8, 5.0, 8.0))],
    }
    info = _info("point_source_models", SLAB)
    lonlat = _load("point_source_models", SLAB, "lonlat.npy")
    n = {"WA": 2490, "OR": 821, "CA": 2306}
    start = {"WA": 0, "OR": 2490, "CA": 2490 + 821}
    assert len(lonlat) == sum(n.values())
    for state, branches in expected.items():
        sel = (slab["node_index"] >= start[state]) & (slab["node_index"] < start[state] + n[state])
        for m_lo in [5.0, 6.0, 7.0]:
            rate = slab["rate"][sel & (slab["m"] >= m_lo)].sum()
            assert rate == pytest.approx(total(branches, m_lo), rel=1e-6), (state, m_lo)
            assert info["total_rates_by_state"][state][f"M>={m_lo}"] == pytest.approx(rate)
    # Known values: WA M>=5 about 0.114/yr, CA M>=5 about 0.0275/yr, OR M>=6.5 about 0.0019/yr
    assert total(expected["WA"], 5.0) == pytest.approx(0.1142, rel=1e-3)
    assert total(expected["CA"], 5.0) == pytest.approx(0.02754, rel=1e-3)
    assert total(expected["OR"], 5.0) == pytest.approx(0.001916, rel=1e-3)


def test_slab_point_source_data(slab):
    depth = _load("point_source_models", SLAB, "depth.npy")
    np.testing.assert_allclose(slab["depth"], depth[slab["node_index"]])
    assert np.all((depth > 4.0) & (depth < 240.0))
    assert set(np.unique(slab["style"])) == {3}
    assert np.all((slab["m"] > 5.0) & (slab["m"] < 8.0))
    table = _load("point_source_models", SLAB, "rjb_geomatrix.npy")
    assert table.shape == (26, 1001)
    # First rows of rjb_geomatrix.dat (M 6.05): 0 -> 0.00, 1 -> 0.64, 2 -> 1.27 km
    np.testing.assert_allclose(table[0, :3], [0.0, 0.64, 1.27])
    assert np.all(np.diff(table, axis=1) >= 0.0)
    assert np.all(table <= np.arange(1001) + 1e-9)


def test_slab_nodes_match_pdf_files(slab):
    pd = pytest.importorskip("pandas")
    grid = UTILITIES / "nshm23_cascadia_slab" / "grid-data"
    if not grid.exists():
        pytest.skip("slab input files not available")
    df = pd.concat([pd.read_csv(grid / f"pdf-{s}.csv") for s in ["wa", "or", "ca"]])
    lonlat = _load("point_source_models", SLAB, "lonlat.npy")
    depth = _load("point_source_models", SLAB, "depth.npy")
    np.testing.assert_allclose(lonlat[:, 0], df["lon"].values)
    np.testing.assert_allclose(lonlat[:, 1], df["lat"].values)
    np.testing.assert_allclose(depth, df["depth"].values)
    # Spot check one node: the WA lo branch rate at M 5.05 is R * pdf * 10^(-0.4 * 5.05)
    i = int(np.argmax(df["pdf"].values[:2490]))
    sel = (slab["node_index"] == i) & np.isclose(slab["m"], 5.05)
    assert slab["rate"][sel].sum() == pytest.approx(
        1.20769448734618 * df["pdf"].values[i] * 10.0 ** (-0.4 * 5.05), rel=1e-9
    )

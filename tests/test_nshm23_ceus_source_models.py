"""
Tests for the NSHM 2023 (nshm-conus 6.2.0) stable-crust source models produced by
utilities/convert_nshm23_ceus_fault_source_models.py and convert_nshm23_ceus_point_source_models.py.

Expected values come from (1) hand calculations from the NSHM logic tree files (rate trees, source
tree weights, zone rate files) and (2) the USGS NSHM rate web service
(https://earthquake.usgs.gov/ws/nshmp/conus-2023/dynamic/rate/<lon>/<lat>/<distance>), which reports
the annual rate of earthquakes within a distance of a site by source component.
"""

import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import pytest

from ucla_plha import plha
from ucla_plha.geometry import geometry

SOURCE_MODELS = Path(str(files("ucla_plha").joinpath("source_models")))
FAULT = SOURCE_MODELS / "fault_source_models"
POINT = SOURCE_MODELS / "point_source_models"
MODELS = {
    "nshm23_ceus_fault": (FAULT, "Fault", False),
    "nshm23_ceus_fault_cluster": (FAULT, "FaultCluster", True),
    "nshm23_ceus_zone": (POINT, "Zone", False),
    "nshm23_ceus_grid": (POINT, "Grid", False),
    "nshm23_ceus_grid_system": (POINT, "Grid", False),
}
SOURCE_IDS = {
    "Cheraw": 2189,
    "Commerce": 10007,
    "Meers": 3619,
    "ERM north": 3800,
    "ERM south": 3809,
    "New Madrid": 3099,
}


def _npz(model, name):
    return dict(np.load(MODELS[model][0] / model / name))


def _mean(tree):
    return sum(w * r for w, r in tree)


def horz_distance_fast(lon1, lat1, lon2, lat2):
    """nshmp-lib Locations.horzDistanceFast (km)."""
    lon1, lat1, lon2, lat2 = (np.radians(v) for v in (lon1, lat1, lon2, lat2))
    dlat = lat1 - lat2
    dlon = (lon1 - lon2) * np.cos((lat1 + lat2) * 0.5)
    return 6371.0072 * np.sqrt(dlat**2 + dlon**2)


def _point_rate_within(model, lon, lat, dist):
    r = _npz(model, "ruptures.npz")
    xyz = np.load(POINT / model / "points.npy")
    rr = np.linalg.norm(xyz, axis=1)
    plat = np.degrees(np.arcsin(xyz[:, 2] / rr))
    plon = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    near = horz_distance_fast(lon, lat, plon, plat) <= dist
    return r["rate"][near[r["node_index"]]].astype(float).sum()


@pytest.mark.parametrize("model", sorted(MODELS))
def test_source_info(model):
    directory, component, cluster = MODELS[model]
    info = json.loads((directory / model / "source_info.json").read_text(encoding="utf-8"))
    assert info["name"] == model
    assert info["tectonic_region"] == "stable_crust"
    assert info["nshm_component"] == component
    assert info["cluster"] is cluster
    assert info["source"].startswith("nshm-conus 6.2.0")
    if directory == POINT:
        assert "point_source_type" in info["point_source"]
        assert "depth" in info["point_source"]


def test_fault_source_rates_match_hand_calculation():
    r = _npz("nshm23_ceus_fault", "ruptures.npz")

    def total(name):
        return r["rate"][r["source_id"] == SOURCE_IDS[name]].sum()

    # SINGLE magnitude trees (weights sum to 1) x rate trees, no epistemic or aleatory terms
    commerce = 0.75 * _mean(
        [(0.101, 2.5e-4), (0.244, 1.4e-4), (0.310, 8.0e-5), (0.244, 4.0e-5), (0.101, 1.4e-5)]
    ) + 0.25 * _mean(
        [(0.101, 3.3e-4), (0.244, 2.0e-4), (0.310, 1.3e-4), (0.244, 7.6e-5), (0.101, 3.4e-5)]
    )
    assert total("Commerce") == pytest.approx(commerce, rel=1e-9)
    erm_north = 0.9 * _mean(
        [(0.101, 2.9e-4), (0.244, 1.5e-4), (0.310, 8.0e-5), (0.244, 4.0e-5), (0.101, 1.4e-5)]
    ) + 0.1 * _mean(
        [(0.101, 3.9e-4), (0.244, 2.2e-4), (0.310, 1.3e-4), (0.244, 7.2e-5), (0.101, 3.2e-5)]
    )
    assert total("ERM north") == pytest.approx(erm_north, rel=1e-9)
    # Eastern Rift Margin (south): Crittenden Co 0.6 and Meeman-Shelby 0.4 share the 2/3/4-eq rate trees
    erm_south = (
        0.333 * _mean([(0.101, 3.5e-4), (0.244, 2.1e-4), (0.310, 1.4e-4), (0.244, 8.0e-5), (0.101, 3.6e-5)])
        + 0.334 * _mean([(0.101, 4.3e-4), (0.244, 2.8e-4), (0.310, 1.9e-4), (0.244, 1.2e-4), (0.101, 6.2e-5)])
        + 0.333 * _mean([(0.101, 5.0e-4), (0.244, 3.4e-4), (0.310, 2.4e-4), (0.244, 1.6e-4), (0.101, 9.0e-5)])
    )
    assert total("ERM south") == pytest.approx(erm_south, rel=1e-6)
    # New Madrid cluster-out: SSCn 0.5 x 0.05 x Reelfoot (short 0.7, extended 0.3); USGS 0.5 x 0.2 x
    # (west, mid-west, center, mid-east, east)
    nm = 0.5 * 0.05 * _mean(
        [(0.101, 1.3e-3), (0.244, 7.2e-4), (0.310, 4.2e-4), (0.244, 2.2e-4), (0.101, 8.05940594e-05)]
    ) + 0.5 * 0.2 * _mean([(0.45, 0.002), (0.45, 0.002), (0.05, 0.001), (0.05, 2.0e-5)])
    assert total("New Madrid") == pytest.approx(nm, rel=1e-6)
    # Moment balancing preserves the moment rate of the SINGLE branches; the epistemic and aleatory
    # magnitude spread changes the event rate. Meers: USGS 0.5 x 2.22e-4 (+-0.2, Gaussian) and SSCn
    # 0.5 x (cluster-in 0.8, cluster-out 0.2) rate trees (Gaussian).
    assert total("Meers") == pytest.approx(4.14644e-4, rel=1e-5)
    assert total("Cheraw") == pytest.approx(1.31434e-4, rel=1e-5)


def test_fault_ruptures_moment_balance_meers_usgs():
    # Meers USGS branch: M7.0 SINGLE at 2.22e-4/yr with epistemic +-0.2 (0.2, 0.6, 0.2) and Gaussian
    # aleatory variability, each moment balanced: total moment rate = 0.5 x 2.22e-4 x M0(7.0) plus
    # the SSCn branch (SINGLE rates x magnitude weights, each moment balanced).
    r = _npz("nshm23_ceus_fault", "ruptures.npz")
    s = r["source_id"] == SOURCE_IDS["Meers"]
    mo = np.sum(r["rate"][s] * 10.0 ** (1.5 * r["m"][s] + 9.05))
    mags = [(0.10, 6.6), (0.45, 6.7), (0.30, 6.9), (0.10, 7.3), (0.05, 7.4)]
    m0 = sum(w * 10.0 ** (1.5 * m + 9.05) for w, m in mags)
    cin = _mean([(0.101, 2.1e-3), (0.244, 1.2e-3), (0.310, 6.7e-4), (0.244, 3.4e-4), (0.101, 1.2e-4)])
    cout = _mean([(0.333, 5.0e-6), (0.334, 2.9e-6), (0.333, 4.94414e-6)])
    expected = 0.5 * 2.22e-4 * 10.0 ** (1.5 * 7.0 + 9.05) + 0.5 * (0.8 * cin + 0.2 * cout) * m0
    assert mo == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize(
    "lon, lat, dist, expected",
    [
        (-103.3, 38.3, 60.0, 1.31434e-04),  # Cheraw
        (-98.5, 34.8, 50.0, 4.14644e-04),  # Meers
        (-89.6, 36.6, 300.0, 6.17420e-04),  # Commerce, ERM, New Madrid cluster-out
    ],
)
def test_fault_rates_match_usgs_rate_service(lon, lat, dist, expected):
    p = geometry.point_to_xyz(np.array([lat, lon, 0.0]))
    out = plha.get_source_data("fault_source_models", "nshm23_ceus_fault", p, dist, None, ["bssa14"])
    assert out[2].sum() == pytest.approx(expected, rel=1e-5)


def test_cluster_model_structure():
    r = _npz("nshm23_ceus_fault_cluster", "ruptures.npz")
    c = _npz("nshm23_ceus_fault_cluster", "clusters.npz")
    assert len(c["cluster_id"]) == 18  # 8 SSCn + 5 USGS x (all, center-south)
    np.testing.assert_array_equal(c["cluster_id"], np.arange(18))
    assert set(np.unique(r["cluster_id"])) == set(range(18))
    # rupture "rate" is the magnitude branch weight within a section: sums to 1 per section
    sections = np.unique(r["cluster_section"])
    sums = np.array([r["rate"][r["cluster_section"] == s].sum() for s in sections])
    np.testing.assert_allclose(sums, 1.0, rtol=1e-9)
    # each section belongs to one cluster; clusters have 2 or 3 sections
    for s in sections:
        assert len(np.unique(r["cluster_id"][r["cluster_section"] == s])) == 1
    n = [len(np.unique(r["cluster_section"][r["cluster_id"] == i])) for i in range(18)]
    np.testing.assert_array_equal(n, c["n_sections"])
    assert set(n) == {2, 3}
    # SSCn clusters: weight 0.5 x 0.9 x cluster-set weight; rate = rate-tree mean
    sscn = c["source_id"] == 3099
    sscn_rate = _mean([(0.101, 0.006), (0.244, 0.0037), (0.310, 0.0024), (0.244, 0.0014), (0.101, 0.000619802)])
    usgs_all = _mean([(0.45, 0.002), (0.45, 0.0013333), (0.05, 0.001), (0.05, 2.0e-05)])
    usgs_cs = 0.45 * 0.0006667
    expected = 0.5 * 0.9 * sscn_rate + 0.5 * 0.8 * (usgs_all + usgs_cs)
    assert sscn.all()
    assert np.sum(c["weight"] * c["rate"]) == pytest.approx(expected, rel=1e-6)
    assert np.sort(c["weight"])[-1] == pytest.approx(0.5 * 0.7 * 0.8)  # USGS center
    assert c["weight"].sum() == pytest.approx(0.5 * 0.9 + 0.5 * 0.8 * 2)


def test_cluster_rates_match_usgs_rate_service():
    # FaultCluster component within 300 km of New Madrid: sum_c weight_c rate_c sum_r rate_r
    r = _npz("nshm23_ceus_fault_cluster", "ruptures.npz")
    c = _npz("nshm23_ceus_fault_cluster", "clusters.npz")
    total = np.sum(r["rate"] * (c["weight"] * c["rate"])[r["cluster_id"]])
    assert total == pytest.approx(5.68814e-03, rel=1e-5)
    # rate service magnitude bins: index int((m - 4.2) / 0.1), so m = 7.0 falls in the 6.95 bin
    m7 = ((r["m"] - 4.2) / 0.1).astype(int) >= 28
    assert np.sum((r["rate"] * (c["weight"] * c["rate"])[r["cluster_id"]])[m7]) == pytest.approx(
        2.4732e-04 + 8.0759e-04 + 1.4623e-03 + 4.9859e-04 + 1.2561e-03 + 3.4775e-04, rel=1e-4
    )


def test_cluster_source_data():
    p = geometry.point_to_xyz(np.array([35.15, -90.05, 0.0]))  # Memphis
    out = plha.get_source_data("fault_source_models", "nshm23_ceus_fault_cluster", p, None, None, ["bssa14", "cb14"])
    for a in out:
        assert np.all(np.isfinite(a))
    assert np.all(out[3] < 200.0) and np.all(out[4] >= out[3] * (1 - 1e-3))


def test_zone_rates_match_hand_calculation():
    r = _npz("nshm23_ceus_zone", "ruptures.npz")
    # sums of the zone rate files x source tree weights (magnitude weights sum to 1)
    expected = (
        0.5 * (8.88e-05 + 1.08e-04 + 4.42e-04 + 3.53300772e-04 + 1.1466e-03)  # AR zones (active 0.5)
        + 1.68832006e-04  # Wabash Valley
        + 1.372000004e-03  # Charlevoix
        + 0.2 * 1.8903612026e-03 + 0.5 * 1.890359248e-03 + 0.3 * 1.89035894e-03  # Charleston
        + 0.5 * (0.5 * 8.835e-04 + 0.5 * 8.84e-04)  # Central Virginia (active 0.5)
    )
    assert r["rate"].sum() == pytest.approx(expected, rel=1e-8)
    assert set(np.unique(r["style"])) == {3}
    strike = np.load(POINT / "nshm23_ceus_zone" / "strike.npy")
    assert strike.shape == np.load(POINT / "nshm23_ceus_zone" / "node_index.npy").shape
    # only Charlevoix (67 nodes) has no strike
    assert np.isnan(strike).sum() == 67
    # Charleston magnitudes 6.7-7.5 (weights 0.1, 0.25, 0.3, 0.25, 0.1)
    near = _point_rate_within("nshm23_ceus_zone", -80.0, 32.9, 300.0)
    assert near == pytest.approx(1.89036e-03, rel=1e-5)


@pytest.mark.parametrize(
    "model, lon, lat, dist, expected",
    [
        ("nshm23_ceus_zone", -90.5, 35.5, 150.0, 4.96050e-04),
        ("nshm23_ceus_zone", -70.2, 47.5, 100.0, 1.37200e-03),
        ("nshm23_ceus_grid", -80.0, 32.9, 300.0, 3.74511e-02),
        ("nshm23_ceus_grid", -90.5, 35.5, 150.0, 8.88089e-02),
        ("nshm23_ceus_grid", -70.2, 47.5, 100.0, 3.57707e-02),
        ("nshm23_ceus_grid", -97.0, 38.0, 300.0, 1.20315e-02),
        ("nshm23_ceus_grid", -83.0, 40.0, 1000.0, 5.73014e-01),
        # Montana: service Grid 2.79460e-03 = ceus-stable 9.5662e-05 + system-stable
        ("nshm23_ceus_grid_system", -106.0, 47.0, 200.0, 2.79460e-03 - 9.5662e-05),
    ],
)
def test_point_rates_match_usgs_rate_service(model, lon, lat, dist, expected):
    assert _point_rate_within(model, lon, lat, dist) == pytest.approx(expected, rel=2e-5)


def test_grid_model_totals_and_magnitudes():
    r = _npz("nshm23_ceus_grid", "ruptures.npz")
    rate = r["rate"].astype(float)
    assert rate.sum() == pytest.approx(0.8992453, rel=1e-6)
    assert rate[r["m"] >= 5.0].sum() == pytest.approx(0.4758039, rel=1e-6)
    m = np.unique(r["m"])
    np.testing.assert_allclose(m, np.round(4.75 + 0.1 * np.arange(33), 2), atol=1e-9)
    assert set(np.unique(r["style"])) == {3}
    rjb = _npz("nshm23_ceus_grid", "rjb_correction.npz")
    assert rjb["rjb"].shape == (26, 1001)
    # the corrected rjb is never larger than the epicentral distance
    assert np.all(rjb["rjb"] <= rjb["r"][None, :] + 1e-9)


def test_grid_system_model():
    r = _npz("nshm23_ceus_grid_system", "ruptures.npz")
    assert set(np.unique(r["style"])) <= {1, 2, 3}
    assert r["rate"].sum() == pytest.approx(0.1849000, rel=1e-5)
    info = json.loads((POINT / "nshm23_ceus_grid_system" / "source_info.json").read_text(encoding="utf-8"))
    assert info["gmm_max_distance_km"] == 300.0
    assert sum(b["weight"] for b in info["gmm_tree"]) == pytest.approx(1.0, abs=1e-3)


def _lon_lat(model):
    xyz = np.load(POINT / model / "points.npy")
    lat = np.degrees(np.arctan2(xyz[:, 2], np.hypot(xyz[:, 0], xyz[:, 1])))
    lon = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    return set(zip(np.round(lon, 3), np.round(lat, 3)))


def test_wus_grid_and_ceus_grid_system_partition_the_nodes():
    # nshmp-lib assigns each node of the WUS fault system grid (branch-avg-grid.csv, 68,883
    # nodes) to either the active (grid-system-active) or the stable (grid-system-stable)
    # region with java.awt.geom.Area.contains, so no node is in both models
    active = _lon_lat("nshm23_wus_grid")
    stable = _lon_lat("nshm23_ceus_grid_system")
    assert (len(active), len(stable)) == (54997, 13886)
    assert not active & stable
    assert len(active | stable) == 68883
    # nodes on the region boundary that are in the stable region in nshmp-lib
    assert (-105.5, 33.7) in stable and (-105.5, 33.7) not in active
    assert (-105.0, 38.9) in stable and (-105.0, 38.9) not in active

"""
Tests of tectonic regions (source_info.json), per-region ground motion model logic trees, the
pygmm ground motion model adapter, nshmp-lib point source distances, and cluster hazard, using
synthetic source models and hand calculations.
"""

import json
import shutil
from importlib.resources import files
from pathlib import Path

import numpy as np
import pytest
from scipy.special import ndtr

from ucla_plha import plha, pygmm_gmms, tectonic_regions
from ucla_plha.geometry import geometry, point_source

pygmm = pytest.importorskip("pygmm")
if not hasattr(pygmm, "NgaEastSeeds") or not hasattr(pygmm.KuehnEtAl2020, "PARAMS"):
    pytest.skip(
        "requires the vectorized pygmm with the nshmp-lib models on PYTHONPATH",
        allow_module_level=True,
    )

PGA = [0.05, 0.1, 0.2, 0.4]
SITE = {"latitude": 36.0, "longitude": -90.0, "elevation": 0.0, "vs30": 760.0}
SIGMA = 0.6


# ---------------------------------------------------------------------------------------
# Synthetic source models
# ---------------------------------------------------------------------------------------


def _xyz(lat, lon):
    return geometry.point_to_xyz(np.array([lat, lon, 0.0]))


def _write_fault_model(root, name, info, n_rup=4, extra=None, clusters=None):
    """A fault model with one vertical segment per rupture, 10 km north of the site."""
    d = root / "source_models" / "fault_source_models" / name
    d.mkdir(parents=True)
    lat0 = SITE["latitude"] + 10.0 / 111.2
    lon0, lon1 = SITE["longitude"] - 0.1, SITE["longitude"] + 0.1
    top0, top1 = _xyz(lat0, lon0), _xyz(lat0, lon1)
    up = top0 / np.linalg.norm(top0)
    bot0, bot1 = top0 - 10.0 * up, top1 - 10.0 * up
    tri = np.array([[top0, top1, bot0], [top1, bot1, bot0]] * n_rup)
    np.save(d / "tri_segment_id.npy", np.repeat(np.arange(n_rup), 2))
    np.save(d / "tri_rjb.npy", tri)
    np.save(d / "tri_rrup.npy", tri)
    np.save(d / "rect_segment_id.npy", np.arange(n_rup))
    np.save(d / "rect_rjb.npy", np.array([[top0, top1, bot1, bot0]] * n_rup))
    np.savez(d / "ruptures_segments.npz", segment_index=np.arange(n_rup), rupture_index=np.arange(n_rup))
    rup = dict(
        m=np.linspace(6.0, 7.5, n_rup),
        fault_type=np.full(n_rup, 3),
        rate=np.full(n_rup, 1e-3),
        dip=np.full(n_rup, 90.0),
        ztor=np.zeros(n_rup),
        zbor=np.full(n_rup, 10.0),
    )
    rup.update(extra or {})
    np.savez(d / "ruptures.npz", **rup)
    if clusters is not None:
        np.savez(d / "clusters.npz", **clusters)
    (d / "source_info.json").write_text(json.dumps(info))
    return d


def _write_point_model(root, name, info, m, style, rate, node_lat, node_lon, node_of_rup, depth=None):
    d = root / "source_models" / "point_source_models" / name
    d.mkdir(parents=True)
    np.save(d / "points.npy", np.array([_xyz(a, b) for a, b in zip(node_lat, node_lon)]))
    np.save(d / "node_index.npy", np.arange(len(node_lat)))
    rup = dict(node_index=np.asarray(node_of_rup), m=np.asarray(m, float), style=np.asarray(style), rate=np.asarray(rate, float))
    if depth is not None:
        rup["depth"] = np.asarray(depth, float)
    np.savez(d / "ruptures.npz", **rup)
    (d / "source_info.json").write_text(json.dumps(info))
    return d


@pytest.fixture
def package(tmp_path, monkeypatch):
    """A package root in tmp_path with the schema, used through plha.files."""
    real = Path(str(files("ucla_plha")))
    shutil.copy(real / "ucla_plha_schema.json", tmp_path / "ucla_plha_schema.json")
    shutil.copytree(real / "geometry", tmp_path / "geometry")
    monkeypatch.setattr(plha, "files", lambda _package: tmp_path)
    return tmp_path


def _hazard(tmp_path, config):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    return plha.get_hazard(str(path))


# ---------------------------------------------------------------------------------------
# Ground motion model logic trees
# ---------------------------------------------------------------------------------------


def test_flat_config_is_active_crust():
    config = {"ground_motion_models": {"bssa14": {"weight": 1.0}, "ask14": {"weight": 3.0}}}
    trees, notes = tectonic_regions.parse_ground_motion_models(config, {"active_crust"})
    assert list(trees) == ["active_crust"]
    assert [(b.key, b.weight) for b in trees["active_crust"]] == [("bssa14", 0.25), ("ask14", 0.75)]
    # normalized in place, as before
    assert config["ground_motion_models"]["ask14"]["weight"] == 0.75
    assert notes == []


def test_region_config_normalizes_within_regions():
    config = {
        "ground_motion_models": {
            "active_crust": {"cb14": {"weight": 2.0}, "cy14": {"weight": 2.0}},
            "stable_crust": {"nga_east_2026": {"weight": 1.0}, "nga_east_seeds_2026": {"weight": 3.0}},
        }
    }
    trees, _ = tectonic_regions.parse_ground_motion_models(config, {"active_crust", "stable_crust"})
    assert [b.weight for b in trees["active_crust"]] == [0.5, 0.5]
    assert [b.weight for b in trees["stable_crust"]] == [0.25, 0.75]
    assert trees["stable_crust"][0].spec.cls is pygmm.NgaEast


def test_mixed_config_is_rejected():
    with pytest.raises(ValueError):
        tectonic_regions.is_region_config({"active_crust": {}, "bssa14": {"weight": 1.0}})


def test_default_interface_tree_weights():
    config = {"ground_motion_models": {"active_crust": {"bssa14": {"weight": 1.0}}}}
    trees, notes = tectonic_regions.parse_ground_motion_models(config, {"subduction_interface"})
    tree = {b.spec.name: b.weight for b in trees["subduction_interface"]}
    assert tree == pytest.approx(tectonic_regions.DEFAULT_GMM_TREES["subduction_interface"])
    assert tree["AM_09_INTERFACE_BASIN"] == pytest.approx(0.125)
    assert any("default" in n for n in notes)


def test_default_slab_tree_weights():
    trees, _ = tectonic_regions.parse_ground_motion_models(
        {"ground_motion_models": {"subduction_slab": "default"}}, {"subduction_slab"}
    )
    tree = {b.spec.name: b.weight for b in trees["subduction_slab"]}
    np.testing.assert_allclose(tree["ZHAO_06_SLAB_BASIN"], 0.25)
    np.testing.assert_allclose(tree["AG_20_CASCADIA_SLAB_BASIN"], 0.0825)
    np.testing.assert_allclose(sum(tree.values()), 1.0)


@pytest.fixture
def pygmm_without_am09_zhao06(monkeypatch):
    """A registry without Atkinson and Macias (2009) and Zhao et al. (2006), as for a pygmm
    older than 1b30f5b."""
    reg = pygmm_gmms.registry()
    ids = {k: v for k, v in reg["ids"].items() if k not in pygmm_gmms.UNAVAILABLE}
    monkeypatch.setattr(pygmm_gmms, "_REGISTRY", {"ids": ids, "classes": reg["classes"]})


def test_default_tree_drops_unavailable_models_and_renormalizes(pygmm_without_am09_zhao06):
    config = {"ground_motion_models": {"active_crust": {"bssa14": {"weight": 1.0}}}}
    with pytest.warns(UserWarning, match="AM_09"):
        trees, notes = tectonic_regions.parse_ground_motion_models(config, {"subduction_interface"})
    tree = {b.spec.name: b.weight for b in trees["subduction_interface"]}
    # AM_09 and ZHAO_06 (0.125 each) removed, the three NGA-Subduction models 1/3 each
    assert set(tree) == {
        "AG_20_CASCADIA_INTERFACE_ADJUSTED_BASIN",
        "KBCG_20_CASCADIA_INTERFACE_BASIN",
        "PSBAH_20_CASCADIA_INTERFACE_BASIN",
    }
    np.testing.assert_allclose(list(tree.values()), 1.0 / 3.0)


# ---------------------------------------------------------------------------------------
# pygmm adapter
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "gmm_id, cls, options, scenario",
    [
        ("NGA_EAST_2026", "NgaEast", {"version": "2026"}, {}),
        ("NGA_EAST_SEEDS_2026_ADJUSTED", "NgaEastSeeds", {"version": "2026", "adjusted": True}, {}),
        ("KBCG_20_CASCADIA_INTERFACE_BASIN", "KuehnEtAl2020", {"basin": True}, {"event_type": "interface", "region": "cascadia"}),
        ("AG_20_CASCADIA_SLAB_ADJUSTED_BASIN", "AbrahamsonGulerce2020", {"adjusted": True, "basin": True}, {"event_type": "intraslab", "region": "cascadia"}),
        ("PSBAH_20_CASCADIA_SLAB_BASIN", "ParkerEtAl2020", {"basin": True}, {"event_type": "intraslab", "region": "cascadia"}),
        ("AM_09_INTERFACE_BASIN", "AtkinsonMacias2009", {"basin": True}, {}),
        ("ZHAO_06_INTERFACE_BASIN", "ZhaoEtAl2006", {"basin": True}, {"event_type": "interface"}),
        ("ZHAO_06_SLAB_BASIN", "ZhaoEtAl2006", {"basin": True}, {"event_type": "intraslab"}),
    ],
)
def test_resolve_nshmp_ids(gmm_id, cls, options, scenario):
    spec = pygmm_gmms.resolve(gmm_id.lower())
    assert spec.kind == "pygmm" and spec.cls is getattr(pygmm, cls)
    assert spec.options == options
    for k, v in scenario.items():
        assert spec.scenario[k] == v


def test_resolve_names():
    assert pygmm_gmms.resolve("bssa14").kind == "native"
    assert pygmm_gmms.resolve("parkeretal2020").cls is pygmm.ParkerEtAl2020
    assert pygmm_gmms.resolve("zhao_06_slab_basin").cls is pygmm.ZhaoEtAl2006
    with pytest.raises(ValueError):
        pygmm_gmms.resolve("not_a_model")


def test_resolve_unavailable(pygmm_without_am09_zhao06):
    with pytest.raises(pygmm_gmms.UnavailableGmmError):
        pygmm_gmms.resolve("zhao_06_slab_basin")


RUPTURES = dict(
    m=np.array([5.5, 6.5, 7.5]),
    fault_type=np.array([3, 1, 2]),
    rjb=np.array([5.0, 20.0, 80.0]),
    rrup=np.array([8.0, 25.0, 82.0]),
    rx=np.array([-5.0, 22.0, 10.0]),
    ry0=np.zeros(3),
    dip=np.array([90.0, 45.0, 60.0]),
    ztor=np.array([2.0, 4.0, 40.0]),
    zbor=np.array([12.0, 18.0, 55.0]),
)


def test_adapter_matches_scalar_pygmm_nga_east():
    spec = pygmm_gmms.resolve("NGA_EAST_SEEDS_2026")
    mu, sigma = pygmm_gmms.get_ground_motion(spec, RUPTURES, {"vs30": 500.0})
    for i in range(3):
        s = pygmm.Scenario(mag=RUPTURES["m"][i], dist_rup=RUPTURES["rrup"][i], dist_jb=RUPTURES["rjb"][i], v_s30=500.0)
        ref = pygmm.NgaEastSeeds(s, version="2026")
        np.testing.assert_allclose(mu[i], np.log(ref.pga), rtol=1e-10)
        np.testing.assert_allclose(sigma[i], ref.ln_std_pga, rtol=1e-10)


def test_adapter_subduction_scenario_mapping():
    spec = pygmm_gmms.resolve("PSBAH_20_CASCADIA_SLAB_BASIN")
    site = {"vs30": 400.0, "z2p5": 3.0}
    values = pygmm_gmms.build_scenario(spec, RUPTURES, site, "subduction_slab")
    assert values["event_type"] == "intraslab" and values["region"] == "cascadia"
    np.testing.assert_array_equal(values["depth_tor"], RUPTURES["ztor"])
    assert values["depth_2_5"] == 3.0
    assert "dist_jb" not in values
    mu, _ = pygmm_gmms.get_ground_motion(spec, RUPTURES, site, "subduction_slab")
    ref = pygmm.ParkerEtAl2020(
        pygmm.Scenario(mag=7.5, dist_rup=82.0, depth_tor=40.0, depth_2_5=3.0, v_s30=400.0, event_type="intraslab", region="cascadia"),
        basin=True,
    )
    np.testing.assert_allclose(mu[2], np.log(ref.pga), rtol=1e-10)


def test_adapter_crustal_mechanism_and_width():
    spec = pygmm_gmms.resolve("chiouyoungs2014")
    values = pygmm_gmms.build_scenario(spec, RUPTURES, {"vs30": 760.0, "measured_vs30": True})
    assert list(values["mechanism"]) == ["SS", "RS", "NS"]
    assert list(values["on_hanging_wall"]) == [False, True, True]
    assert values["vs_source"] == "measured"
    assert "depth_1_0" not in values  # None: model default


def test_distance_needs():
    assert plha._distance_needs(["bssa14"]) == {"rjb"}
    assert plha._distance_needs(["kbcg_20_cascadia_slab_basin"]) == {"rrup"}
    assert plha._distance_needs(["nga_east_2026"]) == {"rjb", "rrup"}
    assert plha._distance_needs(["zhao_06_slab_basin"]) == {"rrup"}
    assert plha._distance_needs(["am_09_interface_basin"]) == {"rrup"}


# ---------------------------------------------------------------------------------------
# nshmp-lib point source distances
# ---------------------------------------------------------------------------------------


def test_point_source_distance_table_lookup():
    table = np.load(str(files("ucla_plha").joinpath("geometry/nshmp_rjb_tables.npz")))["nshm_somerville"]
    m = np.array([5.95, 6.0, 6.1, 6.62, 9.0, 7.0])
    r = np.array([12.7, 12.7, 12.7, 55.2, 30.0, 1500.0])
    out = point_source.point_source_distance(m, r, "nshm_somerville")
    # M < 6: unchanged; Java Math.round((6.1 - 6.05) / 0.1) = 0 in floating point (0.49999)
    np.testing.assert_allclose(out, [12.7, table[0, 12], table[0, 12], table[6, 55], table[25, 30], table[10, 1000]])
    np.testing.assert_allclose(point_source.point_source_distance(m, r, "none"), r)


def test_horz_distance_fast():
    # 1 degree of latitude
    np.testing.assert_allclose(point_source.horz_distance_fast(35.0, -90.0, 36.0, -90.0), np.radians(1.0) * 6371.0072)


def test_finite_strike_slip_and_hanging_wall_by_hand():
    m = np.array([5.5, 6.5])
    style = np.array([3, 1])
    r = np.array([3.0, 3.0])
    ztor = np.array([5.0, 1.0])
    d = point_source.finite_point_source_distances(m, style, r, ztor, "nshm_point_wc94_length", "finite", max_depth=14.0)
    rjb = d["rjb"][1]  # corrected distance (table) for M 6.5
    assert d["rjb"][0] == 3.0 and rjb < 3.0
    # strike slip (footwall), reverse footwall and hanging wall with half rates
    np.testing.assert_array_equal(d["index"], [0, 1, 1])
    np.testing.assert_array_equal(d["rate_scale"], [1.0, 0.5, 0.5])
    np.testing.assert_allclose(d["rrup"][:2], [np.hypot(3.0, 5.0), np.hypot(rjb, 1.0)])
    np.testing.assert_allclose(d["rx"][:2], [-3.0, -rjb])
    # hanging wall: WC94 point scaling width = min((14 - 1) / sin 50, L / 1.5)
    dip = np.radians(50.0)
    w = min(13.0 / np.sin(dip), 10 ** (-3.22 + 0.69 * 6.5) / 1.5)
    w_h, z_bot = w * np.cos(dip), 1.0 + w * np.sin(dip)
    r_cut = z_bot * np.tan(dip)
    rrup0 = min(np.hypot(w_h, 1.0), z_bot * np.cos(dip))
    expected = (z_bot / np.cos(dip) - rrup0) * rjb / r_cut + rrup0
    np.testing.assert_allclose(d["rrup"][2], expected)
    np.testing.assert_allclose(d["rx"][2], rjb + w_h)
    np.testing.assert_allclose(d["zbor"][2], z_bot)


def test_slab_finite_width():
    d = point_source.finite_point_source_distances(
        np.array([7.0]), np.array([3]), np.array([40.0]), np.array([50.0]), "nshm_sub_geomat_length", "finite", max_width=8.0
    )
    table = np.load(str(files("ucla_plha").joinpath("geometry/nshmp_rjb_tables.npz")))["nshm_sub_geomat_length"]
    np.testing.assert_allclose(d["rjb"], table[10, 40])
    np.testing.assert_allclose(d["rrup"], np.hypot(table[10, 40], 50.0))
    np.testing.assert_allclose(d["zbor"], 58.0)


def test_fixed_strike_rupture():
    # E-W rupture centered on a node 20 km north of the site: rJB = 20 km, rRup = hypot(20, 5)
    lat_node = 36.0 + 20.0 / (np.radians(1.0) * 6371.0072)
    d = point_source.fixed_strike_distances(
        np.array([7.0]), np.array([3]), np.array([lat_node]), np.array([-90.0]), np.array([90.0]),
        36.0, -90.0, np.array([5.0]), "nshm_point_wc94_length", max_depth=22.0,
    )
    # the great-circle trace (33 km long) curves slightly toward the site
    np.testing.assert_allclose(d["rjb"], 20.0, rtol=5e-3)
    np.testing.assert_allclose(d["rrup"], np.hypot(d["rjb"], 5.0))
    np.testing.assert_allclose(abs(d["rx"]), 20.0, rtol=5e-3)


def test_smoothing_conserves_rate():
    idx, lat, lon, scale = point_source.smooth_nodes(36.0, -90.0, np.array([36.05, 40.0]), np.array([-90.0, -90.0]), 4, 40.0, 0.1)
    assert np.sum(idx == 0) == 16 and np.sum(idx == 1) == 1
    np.testing.assert_allclose(np.bincount(idx, weights=scale), [1.0, 1.0])
    np.testing.assert_allclose(sorted(set(np.round(lat[idx == 0] - 36.05, 6))), [-0.0375, -0.0125, 0.0125, 0.0375])


def test_nshmp_point_source_model(package):
    info = {
        "name": "synthetic_slab",
        "tectonic_region": "subduction_slab",
        "nshm_component": "Slab",
        "point_source": {"point_source_type": "FINITE", "rupture_scaling": "NSHM_SUB_GEOMAT_LENGTH", "max_width": 8.0, "depth": "rupture"},
    }
    node_lat = np.array([36.0 + 30.0 / 111.19, 37.0])
    _write_point_model(package, "synthetic_slab", info, [6.5, 7.0, 7.0], [3, 3, 3], [1e-3, 2e-3, 1e-3], node_lat, [-90.0, -90.0], [0, 0, 1], depth=[40.0, 50.0, 60.0])
    p_xyz = _xyz(36.0, -90.0)
    out = plha.get_source_data("point_source_models", "synthetic_slab", p_xyz, 100.0, None, ["kbcg_20_cascadia_slab_basin"])
    m, ft, rate, rjb, rrup, rx, rx1, ry0, dip, ztor, zbor = out
    assert len(m) == 2  # the node 111 km away is beyond the cutoff
    r = point_source.horz_distance_fast(36.0, -90.0, node_lat[0], -90.0)
    np.testing.assert_allclose(rjb, point_source.point_source_distance(m, np.full(2, r), "nshm_sub_geomat_length"))
    np.testing.assert_allclose(rrup, np.hypot(rjb, [40.0, 50.0]))
    np.testing.assert_allclose(ztor, [40.0, 50.0])
    np.testing.assert_allclose(zbor, [48.0, 58.0])


def test_depth_map_and_distance_bins(package):
    info = {
        "name": "synthetic_grid",
        "tectonic_region": "stable_crust",
        "nshm_component": "Grid",
        "point_source": {
            "point-source-type": "FINITE",
            "rupture-scaling": "NSHM_SOMERVILLE",
            "max-depth": 22.0,
            "grid-depth-map": {"all": {"mMin": 4.5, "mMax": 10.0, "depth-tree": [{"weight": 0.4, "value": 2.0}, {"weight": 0.6, "value": 6.0}]}},
            "opt-distance-bin": 5.0,
        },
    }
    lats = 36.0 + np.array([101.0, 102.0, 104.0]) / (np.radians(1.0) * 6371.0072)
    _write_point_model(package, "synthetic_grid", info, [5.0, 5.0, 5.0], [3, 3, 3], [1.0, 2.0, 4.0], lats, [-90.0] * 3, [0, 1, 2])
    out = plha.get_source_data("point_source_models", "synthetic_grid", _xyz(36.0, -90.0), 1000.0, None, ["nga_east_2026"])
    m, ft, rate, rjb = out[:4]
    ztor = out[9]
    # all three nodes are in the 100-105 km bin: one rupture per depth, rates summed and weighted
    np.testing.assert_allclose(sorted(ztor), [2.0, 6.0])
    np.testing.assert_allclose(rjb, 102.5)
    np.testing.assert_allclose(sorted(rate), [0.4 * 7.0, 0.6 * 7.0])


# ---------------------------------------------------------------------------------------
# Hazard: per-region models, backward compatibility, clusters
# ---------------------------------------------------------------------------------------


def _fake_gmm(mu_by_name):
    def ground_motion(spec, rupture, site, region=None):
        n = len(rupture["m"])
        return np.full(n, mu_by_name[spec.name]), np.full(n, SIGMA)

    return ground_motion


def test_each_source_uses_its_region_models(package, monkeypatch):
    _write_fault_model(package, "crust", {"name": "crust", "tectonic_region": "active_crust", "nshm_component": "Fault"}, n_rup=1)
    _write_fault_model(package, "stable", {"name": "stable", "tectonic_region": "stable_crust", "nshm_component": "Fault"}, n_rup=1)
    mu = {"NGA_EAST_2026": np.log(0.3), "NGA_EAST_SEEDS_2026": np.log(0.1)}
    monkeypatch.setattr(pygmm_gmms, "get_ground_motion", _fake_gmm(mu))
    monkeypatch.setattr(plha, "get_ground_motion_data", lambda gmm, *a: (np.full(len(a[5]), np.log(0.2)), np.full(len(a[5]), SIGMA)))
    config = {
        "site": SITE,
        "source_models": {"fault_source_models": {"crust": {"weight": 1.0}, "stable": {"weight": 1.0}}},
        "ground_motion_models": {
            "active_crust": {"bssa14": {"weight": 1.0}},
            "stable_crust": {"nga_east_2026": {"weight": 3.0}, "nga_east_seeds_2026": {"weight": 1.0}},
        },
        "output": {"psha": {"pga": PGA, "source_model_hazard": True}},
    }
    out = _hazard(package, config)["output"]["psha"]
    exc = lambda med: 1e-3 * (1.0 - ndtr((np.log(PGA) - np.log(med)) / SIGMA))
    # different regions are additive (each its own weight group)
    expected = exc(0.2) + 0.75 * exc(0.3) + 0.25 * exc(0.1)
    np.testing.assert_allclose(out["annual_rate_of_exceedance"], expected, rtol=1e-12)
    np.testing.assert_allclose(out["source_model_hazard"]["crust"]["annual_rate_of_exceedance"], exc(0.2), rtol=1e-12)
    assert out["source_model_hazard"]["stable"]["tectonic_region"] == "stable_crust"


def test_flat_config_with_old_models_unchanged(package, monkeypatch):
    # Old-format config and a model directory without source_info.json (active crust)
    d = _write_fault_model(package, "old", {}, n_rup=1)
    (d / "source_info.json").unlink()
    monkeypatch.setattr(plha, "get_ground_motion_data", lambda gmm, *a: (np.full(1, np.log(0.2 if gmm == "bssa14" else 0.4)), np.full(1, SIGMA)))
    config = {
        "site": SITE,
        "source_models": {"fault_source_models": {"old": {"weight": 2.0}}},
        "ground_motion_models": {"bssa14": {"weight": 1.0}, "cy14": {"weight": 1.0}},
        "output": {"psha": {"pga": PGA}},
    }
    out = _hazard(package, config)
    exc = lambda med: 1e-3 * (1.0 - ndtr((np.log(PGA) - np.log(med)) / SIGMA))
    np.testing.assert_allclose(out["output"]["psha"]["annual_rate_of_exceedance"], 0.5 * exc(0.2) + 0.5 * exc(0.4), rtol=1e-12)
    assert "notes" not in out and set(out["output"]) == {"psha"}


def test_truncation_level():
    eps, p = plha.get_exceedance(np.array([0.1, 1.0, 10.0]), np.log(np.array([0.1])), np.array([0.5]), 3.0)
    p_hi = ndtr(-3.0)
    np.testing.assert_allclose(p[:, 0], [(0.5 - p_hi) / (1 - p_hi), np.clip((1 - ndtr(np.log(10) / 0.5) - p_hi) / (1 - p_hi), 0, 1), 0.0])


def test_cluster_hazard_hand_calculation():
    # Two clusters: cluster 0 has faults A (two magnitude variants, weights 0.4/0.6) and B
    # (one variant); cluster 1 has fault C. One ground motion level.
    p = np.array([[0.1, 0.3, 0.2, 0.5]])
    weight = np.array([0.4, 0.6, 1.0, 1.0])
    cluster = np.array([0, 0, 0, 1])
    fault = np.array([0, 0, 1, 2])
    rate = np.array([0.01, 0.01, 0.01, 0.002])
    hazard, contrib = plha.get_cluster_hazard(p, weight, cluster, fault, rate)
    p_a = 0.4 * 0.1 + 0.6 * 0.3
    p_c0 = 1 - (1 - p_a) * (1 - 0.2)
    expected = 0.01 * p_c0 + 0.002 * 0.5
    np.testing.assert_allclose(hazard, [expected])
    np.testing.assert_allclose(contrib.sum(), expected)
    np.testing.assert_allclose(contrib[0, 2], 0.01 * p_c0 * 0.2 / (p_a + 0.2))


def test_cluster_model_hazard_and_disaggregation(package, monkeypatch):
    # Cluster model with an independent rupture (cluster_id -1) and one cluster of two sections
    info = {"name": "clus", "tectonic_region": "stable_crust", "nshm_component": "Fault", "cluster": True}
    extra = dict(
        cluster_id=np.array([-1, 0, 0, 0]),
        cluster_section=np.array([-1, 0, 0, 1]),
        rate=np.array([1e-3, 0.5, 0.5, 1.0]),
        m=np.array([6.0, 7.0, 7.2, 7.4]),
    )
    clusters = dict(cluster_id=np.array([0]), rate=np.array([0.002]), weight=np.array([0.5]))
    _write_fault_model(package, "clus", info, n_rup=4, extra=extra, clusters=clusters)
    medians = np.log(np.array([0.1, 0.2, 0.3, 0.4]))

    def ground_motion(spec, rupture, site, region=None):
        return medians[: len(rupture["m"])].copy(), np.full(len(rupture["m"]), SIGMA)

    monkeypatch.setattr(pygmm_gmms, "get_ground_motion", ground_motion)
    bins = {"magnitude_bin_edges": [5.0, 6.5, 8.0], "distance_bin_edges": [0.0, 100.0], "epsilon_bin_edges": [-20.0, 20.0]}
    config = {
        "site": SITE,
        "source_models": {"fault_source_models": {"clus": {"weight": 1.0}}},
        "ground_motion_models": {"stable_crust": {"nga_east_2026": {"weight": 1.0}}},
        "output": {"psha": {"pga": PGA, "disaggregation": bins, "source_model_hazard": True}},
    }
    out = _hazard(package, config)["output"]["psha"]
    p = 1.0 - ndtr((np.log(PGA)[:, None] - medians) / SIGMA)
    p_s0 = 0.5 * p[:, 1] + 0.5 * p[:, 2]
    p_c = 1 - (1 - p_s0) * (1 - p[:, 3])
    expected_cluster = 0.002 * 0.5 * p_c
    expected = 1e-3 * p[:, 0] + expected_cluster
    np.testing.assert_allclose(out["annual_rate_of_exceedance"], expected, rtol=1e-12)
    np.testing.assert_allclose(out["source_model_hazard"]["clus:cluster"]["annual_rate_of_exceedance"], expected_cluster, rtol=1e-12)
    assert out["source_model_hazard"]["clus:cluster"]["nshm_component"] == "FaultCluster"
    disagg = np.array(out["disaggregation"])
    # magnitude bin 6.5-8 holds the cluster share
    np.testing.assert_allclose(disagg[:, 1, 0, 0], 100.0 * expected_cluster / expected, rtol=1e-10)

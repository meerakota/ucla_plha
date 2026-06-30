import copy
import json
import shutil
from importlib.resources import files
from unittest import mock

import numpy as np
import pytest

from ucla_plha import plha

from constants import P_XYZ

FAULT_MODEL = "ucerf3_fm31"
POINT_MODEL = "ucerf3_fm31_grid_sub_seis"

ALL_GMMS = ["bssa14", "ask14", "cb14", "cy14"]


FAULT_CASES = [
    (ALL_GMMS, 60.0, 7.0),
    (["ask14"], 60.0, 7.0),
    (["bssa14"], 60.0, None),
    (["ask14"], 60.0, None),
    (["bssa14"], None, 7.0),
    (["bssa14"], None, None),
]


@pytest.mark.parametrize(
    "gmms,dist_cutoff,m_min",
    FAULT_CASES,
    ids=["all/d/m", "ask/d/m", "bssa/d", "ask/d", "bssa/m", "bssa/none"],
)
def test_get_source_data_fault(gmms, dist_cutoff, m_min):
    out = plha.get_source_data(
        "fault_source_models", FAULT_MODEL, P_XYZ, dist_cutoff, m_min, gmms
    )
    m, fault_type, rate, rjb, rrup, rx, rx1, ry0, dip, ztor, zbor = out

    lengths = {len(a) for a in out}
    assert len(lengths) == 1, "all returned arrays must share one length"
    if m_min is not None and len(m):
        assert m.min() >= m_min
    if dist_cutoff is not None and len(m):
        uses_rjb = any(g in ("bssa14", "cb14", "cy14") for g in gmms)
        distance = rjb if uses_rjb else rrup
        assert distance.max() < dist_cutoff
    assert set(np.unique(fault_type)).issubset({1, 2, 3})


POINT_FILTER_CASES = [
    (60.0, 7.0),
    (60.0, None),
    (None, 7.0),
    (None, None),
]


@pytest.mark.parametrize("dist_cutoff,m_min", POINT_FILTER_CASES)
def test_get_source_data_point_filters(dist_cutoff, m_min):
    out = plha.get_source_data(
        "point_source_models", POINT_MODEL, P_XYZ, dist_cutoff, m_min, ["bssa14"]
    )
    lengths = {len(a) for a in out}
    assert len(lengths) == 1
    m = out[0]
    rjb = out[3]
    if m_min is not None and len(m):
        assert m.min() >= m_min
    if dist_cutoff is not None and len(m):
        assert rjb.max() < dist_cutoff


def test_get_source_data_point_site_on_node(tmp_path, monkeypatch):
    grid = tmp_path / "source_models" / "point_source_models" / POINT_MODEL
    grid.mkdir(parents=True)
    points = np.array([[100.0, 0.0, 0.0], [200.0, 0.0, 0.0], [300.0, 0.0, 0.0]])
    node_ids = np.array([10, 20, 30])
    np.save(str(grid / "points.npy"), points)
    np.save(str(grid / "node_index.npy"), node_ids)
    np.savez(
        str(grid / "ruptures.npz"),
        rate=np.array([0.1, 0.2, 0.3]),
        m=np.array([5.0, 6.0, 7.0]),
        style=np.array([1, 2, 3]),
        node_index=node_ids,
    )
    monkeypatch.setattr(plha, "files", lambda _package: tmp_path)

    on_node = points[1]
    out = plha.get_source_data(
        "point_source_models", POINT_MODEL, on_node, None, None, ["bssa14"]
    )
    rjb = out[3]
    assert (rjb == 0).any()
    assert np.all(np.isfinite(out[4]))


def test_get_source_data_no_distance_models():
    out = plha.get_source_data("fault_source_models", FAULT_MODEL, P_XYZ, None, None, [])
    assert len({len(a) for a in out}) == 1
    assert len(out[0]) > 0


def test_get_source_data_unknown_type_returns_none():
    assert (
        plha.get_source_data(
            "bogus_source_models", FAULT_MODEL, P_XYZ, None, None, ["bssa14"]
        )
        is None
    )


def test_get_disagg_bins_and_conserves_hazard():
    m = np.array([5.2, 5.8, 6.4, 6.9, 7.2, 7.8])
    r = np.array([10.0, 30.0, 40.0, 60.0, 80.0, 95.0])
    eps = np.array([[-2.0, -1.0, 0.5, 1.0, 2.0, -0.5], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    hazards = np.array(
        [
            [2.3e-3, 1.1e-3, 8.4e-4, 3.5e-4, 1.7e-4, 6.2e-5],
            [4.1e-4, 1.9e-4, 1.4e-4, 5.8e-5, 2.6e-5, 9.3e-6],
        ]
    )
    m_edges = np.array([5.0, 6.0, 7.0, 8.0])
    r_edges = np.array([0.0, 50.0, 100.0])
    eps_edges = np.array([-3.0, 0.0, 3.0])

    disagg = plha.get_disagg(hazards, m, r, eps, m_edges, r_edges, eps_edges)

    assert disagg.shape == (2, 3, 2, 2)
    np.testing.assert_allclose(disagg.sum(axis=(1, 2, 3)), hazards.sum(axis=1))


def test_decompress_and_decompressed_read_branch(tmp_path, monkeypatch):
    real_pkg = files("ucla_plha")
    for branch in ("ucerf3_fm31", "ucerf3_fm32"):
        dst = tmp_path / "source_models" / "fault_source_models" / branch
        dst.mkdir(parents=True)
        src = real_pkg.joinpath("source_models/fault_source_models/" + branch)
        files_needed = ["ruptures.npz", "ruptures_segments.npz"]
        if branch == "ucerf3_fm31":
            files_needed += ["tri_segment_id.npy", "tri_rjb.npy"]
        for name in files_needed:
            shutil.copy(str(src.joinpath(name)), str(dst / name))

    monkeypatch.setattr(plha, "files", lambda _package: tmp_path)

    plha.decompress_ucerf3_source_data()

    made = sorted(p.name for p in tmp_path.rglob("*.npy.npz"))
    assert made == [
        "ruptures.npy.npz",
        "ruptures.npy.npz",
        "ruptures_segments.npy.npz",
        "ruptures_segments.npy.npz",
    ]

    out = plha.get_source_data(
        "fault_source_models", "ucerf3_fm31", P_XYZ, 200.0, 5.0, ["bssa14"]
    )
    assert len({len(a) for a in out}) == 1
    assert len(out[0]) > 0


SITE = {"latitude": 36.80547, "longitude": -121.786074, "elevation": 0, "vs30": 230.0}
TIGHT = {"dist_cutoff": 60, "m_min": 7.0}
SMALL_BINS = {
    "magnitude_bin_edges": [6.5, 7.0, 7.5, 8.5],
    "distance_bin_edges": [0.0, 30.0, 60.0],
    "epsilon_bin_edges": [-500.0, -1.0, 0.0, 1.0, 500.0],
}


FULL_STACK_CONFIG = {
    "site": SITE,
    "constraints": TIGHT,
    "source_models": {
        "fault_source_models": {
            "ucerf3_fm31": {"weight": 0.5},
            "ucerf3_fm32": {"weight": 0.5},
        },
        "point_source_models": {
            "ucerf3_fm31_grid_sub_seis": {"weight": 0.25},
            "ucerf3_fm31_grid_unassociated": {"weight": 0.25},
            "ucerf3_fm32_grid_sub_seis": {"weight": 0.25},
            "ucerf3_fm32_grid_unassociated": {"weight": 0.25},
        },
    },
    "ground_motion_models": {
        "bssa14": {"weight": 0.25},
        "ask14": {"weight": 0.25},
        "cb14": {"weight": 0.25},
        "cy14": {"weight": 0.0},
    },
    "liquefaction_models": {
        "moss_et_al_2006": {
            "depth": 2.85,
            "qc": 5.49,
            "fs": 0.02617,
            "sigmav": 46.55,
            "sigmavp": 43.22,
            "pa": 101.325,
            "weight": 0.0,
        },
        "boulanger_idriss_2016": {
            "depth": 2.2,
            "qc1ncs": 71.8,
            "sigmav": 40.0,
            "sigmavp": 28.0,
            "pa": 101.325,
            "weight": 0.5,
        },
        "ngl_smt_2024": {
            "ztop": [0.0, 1.0, 5.0],
            "zbot": [1.0, 5.0, 9.0],
            "qc1ncs": [50.0, 150.0, 250.0],
            "ic": [1.5, 2.0, 2.6],
            "sigmav": [20.0, 60.0, 120.0],
            "sigmavp": [15.0, 45.0, 90.0],
            "ksat": [1.0, 1.0, 1.0],
            "pa": 101.325,
            "weight": 0.5,
        },
    },
    "output": {
        "outputfile": "default",
        "psha": {"pga": [0.1, 0.5, 1.0], "disaggregation": SMALL_BINS},
        "plha": {"fsl": [0.1, 0.5, 1.0], "disaggregation": SMALL_BINS},
    },
}


def test_get_hazard_full_stack(tmp_path):
    cfg_path = tmp_path / "full_stack.json"
    cfg_path.write_text(json.dumps(FULL_STACK_CONFIG))

    out = plha.get_hazard(str(cfg_path))

    psha = out["output"]["psha"]
    plha_out = out["output"]["plha"]
    assert len(psha["annual_rate_of_exceedance"]) == 3
    assert len(plha_out["annual_rate_of_nonexceedance"]) == 3
    assert np.array(psha["disaggregation"]).shape == (3, 3, 2, 4)
    assert np.array(plha_out["disaggregation"]).shape == (3, 3, 2, 4)
    assert (tmp_path / "full_stack_output.json").exists()


def test_get_hazard_psha_only_no_disaggregation(tmp_path, monkeypatch):
    config = copy.deepcopy(FULL_STACK_CONFIG)
    del config["source_models"]["point_source_models"]
    del config["source_models"]["fault_source_models"]["ucerf3_fm32"]
    del config["ground_motion_models"]["ask14"]
    del config["ground_motion_models"]["cb14"]
    del config["ground_motion_models"]["cy14"]
    del config["liquefaction_models"]
    del config["output"]["plha"]
    del config["output"]["psha"]["disaggregation"]
    del config["output"]["outputfile"]
    cfg_path = tmp_path / "psha_only.json"
    cfg_path.write_text(json.dumps(config))

    liquefaction = mock.Mock(wraps=plha.get_liquefaction_cdfs)
    monkeypatch.setattr(plha, "get_liquefaction_cdfs", liquefaction)
    out = plha.get_hazard(str(cfg_path))

    assert set(out["output"].keys()) == {"psha"}
    assert "disaggregation" not in out["output"]["psha"]
    liquefaction.assert_not_called()


def test_get_hazard_plha_only_with_named_output(tmp_path, monkeypatch):
    config = copy.deepcopy(FULL_STACK_CONFIG)
    del config["source_models"]["fault_source_models"]
    del config["source_models"]["point_source_models"]["ucerf3_fm32_grid_unassociated"]
    del config["output"]["psha"]
    del config["output"]["plha"]["disaggregation"]
    config["output"]["outputfile"] = "named_out.json"
    cfg_path = tmp_path / "plha_only.json"
    cfg_path.write_text(json.dumps(config))

    seismic = mock.Mock(wraps=plha.ndtr)
    monkeypatch.setattr(plha, "ndtr", seismic)
    monkeypatch.chdir(tmp_path)
    out = plha.get_hazard(str(cfg_path))

    assert set(out["output"].keys()) == {"plha"}
    assert "disaggregation" not in out["output"]["plha"]
    assert (tmp_path / "named_out.json").exists()
    seismic.assert_not_called()


def test_get_hazard_liquefaction_defined_but_psha_only(tmp_path):
    config = copy.deepcopy(FULL_STACK_CONFIG)
    del config["output"]["plha"]
    cfg_path = tmp_path / "psha_unused_liq.json"
    cfg_path.write_text(json.dumps(config))

    out = plha.get_hazard(str(cfg_path))

    assert set(out["output"].keys()) == {"psha"}


def test_get_hazard_all_zero_weight_gmm(tmp_path):
    config = copy.deepcopy(FULL_STACK_CONFIG)
    for gmm in config["ground_motion_models"].values():
        gmm["weight"] = 0.0
    del config["liquefaction_models"]
    del config["output"]["plha"]
    del config["output"]["psha"]["disaggregation"]
    cfg_path = tmp_path / "zero_gmm.json"
    cfg_path.write_text(json.dumps(config))

    out = plha.get_hazard(str(cfg_path))

    assert out["output"]["psha"]["annual_rate_of_exceedance"] == [0.0, 0.0, 0.0]


def test_get_hazard_invalid_config_returns_none(tmp_path, capsys):
    cfg_path = tmp_path / "bad.json"
    cfg_path.write_text("{}")

    result = plha.get_hazard(str(cfg_path))

    assert result is None
    assert "Config File Error" in capsys.readouterr().out

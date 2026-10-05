"""
Tests of the hazard integral in plha.get_hazard against analytic solutions, and of how source
model, ground motion model, and liquefaction model weights are combined.
"""

import copy
import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import pytest
from scipy.special import ndtr

from ucla_plha import plha

SOURCE_MODELS = Path(str(files("ucla_plha").joinpath("source_models")))
PGA = [0.05, 0.1, 0.2, 0.4, 0.8]
FSL = [0.5, 1.0, 1.5]
SITE = {"latitude": 36.80547, "longitude": -121.786074, "elevation": 0, "vs30": 300.0}
BI16 = {"depth": 4.0, "qc1ncs": 90.0, "sigmav": 72.0, "sigmavp": 52.0, "pa": 101.325}

# One synthetic rupture per source model, and a fixed median and standard deviation per GMM
RUPTURE = {"m": 7.0, "rate": 0.01, "rjb": 10.0}
GMM_MU = {
    "bssa14": np.log(0.2),
    "ask14": np.log(0.3),
    "cb14": np.log(0.15),
    "cy14": np.log(0.25),
    "idriss14": np.log(0.22),
}
GMM_DIR = Path(str(files("ucla_plha").joinpath("ground_motion_models")))
SIGMA = 0.6


def _fake_source_data(source_type, source_model, p_xyz, dist_cutoff, m_min, gmms):
    one = lambda v: np.array([v], dtype=float)
    return (
        one(RUPTURE["m"]),
        np.array([3]),
        one(RUPTURE["rate"]),
        one(RUPTURE["rjb"]),
        one(RUPTURE["rjb"]),
        one(RUPTURE["rjb"]),
        one(RUPTURE["rjb"]),
        one(0.0),
        one(90.0),
        one(0.0),
        one(12.0),
    )


def _fake_ground_motion_data(gmm, *args):
    return np.array([GMM_MU[gmm]]), np.array([SIGMA])


@pytest.fixture
def synthetic(monkeypatch):
    monkeypatch.setattr(plha, "get_source_data", _fake_source_data)
    monkeypatch.setattr(plha, "get_ground_motion_data", _fake_ground_motion_data)


def _hazard(tmp_path, config):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    return plha.get_hazard(str(path))["output"]


def _config(fault_models=None, point_models=None, gmms=None, liquefaction=None):
    source_models = {}
    if fault_models:
        source_models["fault_source_models"] = {
            k: {"weight": w} for k, w in fault_models.items()
        }
    if point_models:
        source_models["point_source_models"] = {
            k: {"weight": w} for k, w in point_models.items()
        }
    config = {
        "site": SITE,
        "source_models": source_models,
        "ground_motion_models": {
            k: {"weight": w} for k, w in (gmms or {"bssa14": 1.0}).items()
        },
        "output": {"psha": {"pga": PGA}},
    }
    if liquefaction:
        config["liquefaction_models"] = liquefaction
        config["output"]["plha"] = {"fsl": FSL}
    return config


def _exceedance(gmm):
    return RUPTURE["rate"] * (1.0 - ndtr((np.log(PGA) - GMM_MU[gmm]) / SIGMA))


def test_single_rupture_hazard_matches_analytic(synthetic, tmp_path):
    out = _hazard(tmp_path, _config(fault_models={"ucerf3_fm31": 1.0}))
    np.testing.assert_allclose(
        out["psha"]["annual_rate_of_exceedance"], _exceedance("bssa14"), rtol=1e-12
    )


def test_ground_motion_model_weights_are_normalized(synthetic, tmp_path):
    config = _config(
        fault_models={"ucerf3_fm31": 1.0}, gmms={"bssa14": 1.0, "ask14": 3.0}
    )
    out = _hazard(tmp_path, config)
    expected = 0.25 * _exceedance("bssa14") + 0.75 * _exceedance("ask14")
    np.testing.assert_allclose(
        out["psha"]["annual_rate_of_exceedance"], expected, rtol=1e-12
    )


def test_fault_and_point_source_hazards_add(synthetic, tmp_path):
    # Alternative models within a source type are weighted, and the two source types are summed
    config = _config(
        fault_models={"ucerf3_fm31": 2.0, "ucerf3_fm32": 2.0},
        point_models={
            "ucerf3_fm31_grid_sub_seis": 1.0,
            "ucerf3_fm31_grid_unassociated": 1.0,
        },
    )
    out = _hazard(tmp_path, config)
    np.testing.assert_allclose(
        out["psha"]["annual_rate_of_exceedance"],
        2.0 * _exceedance("bssa14"),
        rtol=1e-12,
    )


def test_every_ground_motion_model_is_registered(synthetic, tmp_path):
    # Each module in ground_motion_models must be accepted by the config schema and have its
    # weight normalized in get_hazard
    gmms = sorted(p.stem for p in GMM_DIR.glob("*.py") if p.stem != "__init__")
    assert set(gmms) == set(GMM_MU)
    for gmm in gmms:
        config = _config(fault_models={"ucerf3_fm31": 1.0}, gmms={gmm: 4.0})
        out = _hazard(tmp_path, config)
        assert out is not None, f"{gmm} rejected by config schema"
        np.testing.assert_allclose(
            out["psha"]["annual_rate_of_exceedance"],
            _exceedance(gmm),
            rtol=1e-12,
            err_msg=gmm,
        )


@pytest.mark.parametrize("source_type", ["fault_source_models", "point_source_models"])
def test_every_source_model_is_registered(synthetic, tmp_path, source_type):
    # Each source model directory must be accepted by the config schema and have its weight
    # normalized in get_hazard. An unregistered model would be rejected by the schema, or its
    # weight would be applied without normalization.
    # Each source model uses the ground motion models of its tectonic region (source_info.json),
    # so bssa14 is given for the region of the model. Cluster models do not have independent
    # ruptures and are tested in test_tectonic_regions.py.
    models = sorted(
        p.name for p in (SOURCE_MODELS / source_type).iterdir() if p.is_dir()
    )
    assert models
    for model in models:
        info = plha.get_source_info(source_type, model)
        if info["cluster"]:
            continue
        key = "fault_models" if source_type == "fault_source_models" else "point_models"
        config = _config(**{key: {model: 4.0}})
        config["ground_motion_models"] = {
            info["tectonic_region"]: config["ground_motion_models"]
        }
        if isinstance(info.get("gmm_tree"), list):
            # source models with their own ground motion model tree (source_info.json)
            config["source_models"][source_type][model]["ground_motion_models"] = {
                "bssa14": {"weight": 1.0}
            }
        out = _hazard(tmp_path, config)
        assert out is not None, f"{model} rejected by config schema"
        np.testing.assert_allclose(
            out["psha"]["annual_rate_of_exceedance"],
            _exceedance("bssa14"),
            rtol=1e-12,
            err_msg=model,
        )


def test_liquefaction_hazard_matches_analytic(synthetic, tmp_path):
    liquefaction = {"boulanger_idriss_2016": dict(BI16, weight=1.0)}
    out = _hazard(
        tmp_path, _config(fault_models={"ucerf3_fm31": 1.0}, liquefaction=liquefaction)
    )

    # Boulanger and Idriss (2016) evaluated by hand for the synthetic rupture
    m, d, q = RUPTURE["m"], BI16["depth"], BI16["qc1ncs"]
    rd = np.exp(
        -1.012
        - 1.126 * np.sin(d / 11.73 + 5.133)
        + (0.106 + 0.118 * np.sin(d / 11.28 + 5.142)) * m
    )
    ln_csr = GMM_MU["bssa14"] + np.log(0.65 * BI16["sigmav"] / BI16["sigmavp"] * rd)
    msf = 1.0 + (min(1.09 + (q / 180.0) ** 3, 2.2) - 1.0) * (
        8.64 * np.exp(-m / 4.0) - 1.325
    )
    c_sigma = min(1.0 / (37.3 - 8.27 * q**0.264), 0.3)
    k_sigma = min(1.0 - c_sigma * np.log(BI16["sigmavp"] / BI16["pa"]), 1.1)
    ln_crr = (
        q / 113
        + (q / 1000) ** 2
        - (q / 140) ** 3
        + (q / 137) ** 4
        - 2.60
        + np.log(msf * k_sigma)
    )
    expected = RUPTURE["rate"] * ndtr(
        (np.log(FSL) - (ln_crr - ln_csr)) / np.sqrt(SIGMA**2 + 0.2**2)
    )

    np.testing.assert_allclose(
        out["plha"]["annual_rate_of_nonexceedance"], expected, rtol=1e-10
    )


REAL_CONFIG = {
    "site": SITE,
    "constraints": {"dist_cutoff": 100, "m_min": 6.0},
    "source_models": {
        "fault_source_models": {"ucerf3_fm31": {"weight": 1.0}},
        "point_source_models": {"ucerf3_fm31_grid_sub_seis": {"weight": 1.0}},
    },
    "ground_motion_models": {"bssa14": {"weight": 1.0}, "ask14": {"weight": 1.0}},
    "liquefaction_models": {"boulanger_idriss_2016": dict(BI16, weight=1.0)},
    "output": {
        "psha": {
            "pga": PGA,
            "disaggregation": {
                "magnitude_bin_edges": [6.0, 7.0, 8.0, 9.0],
                "distance_bin_edges": [0.0, 25.0, 50.0, 100.0],
                "epsilon_bin_edges": [-50.0, -1.0, 0.0, 1.0, 50.0],
            },
        },
        "plha": {"fsl": FSL},
    },
}


@pytest.fixture(scope="module")
def real_output(tmp_path_factory):
    path = tmp_path_factory.mktemp("real") / "config.json"
    path.write_text(json.dumps(REAL_CONFIG))
    return plha.get_hazard(str(path))["output"]


def test_real_hazard_curve_is_monotonic(real_output):
    psha = np.array(real_output["psha"]["annual_rate_of_exceedance"])
    plha_curve = np.array(real_output["plha"]["annual_rate_of_nonexceedance"])
    assert np.all(psha > 0.0)
    assert np.all(np.diff(psha) < 0.0)
    assert np.all(np.diff(plha_curve) > 0.0)


def test_real_disaggregation_sums_to_100_percent(real_output):
    disagg = np.array(real_output["psha"]["disaggregation"])
    np.testing.assert_allclose(disagg.sum(axis=(1, 2, 3)), 100.0, rtol=1e-9)
    assert np.all(disagg >= 0.0)


def test_magnitude_cutoff_reduces_hazard(tmp_path):
    config = copy.deepcopy(REAL_CONFIG)
    del config["liquefaction_models"]
    del config["output"]["plha"]
    del config["output"]["psha"]["disaggregation"]
    low = _hazard(tmp_path, config)["psha"]["annual_rate_of_exceedance"]
    config["constraints"]["m_min"] = 7.0
    high = _hazard(tmp_path, config)["psha"]["annual_rate_of_exceedance"]
    assert np.all(np.array(high) < np.array(low))

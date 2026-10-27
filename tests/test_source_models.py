"""Source models loaded once (plha.load_source_models) and several soil profiles per config."""

import copy
import json

import numpy as np
import pytest

from ucla_plha import plha

SITE = {"latitude": 36.80547, "longitude": -121.786074, "elevation": 0, "vs30": 230.0}
SECOND_SITE = {"latitude": 34.05, "longitude": -118.25, "elevation": 0, "vs30": 300.0}
BINS = {
    "magnitude_bin_edges": [6.5, 7.0, 7.5, 8.5],
    "distance_bin_edges": [0.0, 30.0, 60.0],
    "epsilon_bin_edges": [-500.0, -1.0, 0.0, 1.0, 500.0],
}

LIQUEFACTION = {
    "moss_et_al_2006": {
        "depth": 2.85, "qc": 5.49, "fs": 0.02617, "sigmav": 46.55, "sigmavp": 43.22,
        "pa": 101.325, "weight": 0.3,
    },
    "boulanger_idriss_2016": {
        "depth": 2.2, "qc1Ncs": 71.8, "sigmav": 40.0, "sigmavp": 28.0, "pa": 101.325,
        "weight": 0.3,
    },
    "ngl_smt_2024": {
        "ztop": [0.0, 1.0, 5.0], "zbot": [1.0, 5.0, 9.0], "qc1Ncs": [50.0, 150.0, 250.0],
        "Ic": [1.5, 2.0, 2.6], "sigmav": [20.0, 60.0, 120.0], "sigmavp": [15.0, 45.0, 90.0],
        "Ksat": [0.0, 1.0, 1.0], "pa": 101.325, "weight": 0.4,
    },
}

CONFIG = {
    "site": SITE,
    "constraints": {"dist_cutoff": 60, "m_min": 7.0},
    "source_models": {
        "fault_source_models": {"ucerf3_fm31": {"weight": 1.0}},
        "point_source_models": {"ucerf3_fm31_grid_unassociated": {"weight": 1.0}},
    },
    "ground_motion_models": {"bssa14": {"weight": 0.5}, "ask14": {"weight": 0.5}},
    "liquefaction_models": LIQUEFACTION,
    "output": {"psha": {"pga": [0.05, 0.2, 0.8]}, "plha": {"fsl": [0.5, 1.0, 2.0]}},
}


def _profile_models(scale):
    liq = copy.deepcopy(LIQUEFACTION)
    liq["moss_et_al_2006"]["qc"] *= scale
    liq["boulanger_idriss_2016"]["qc1Ncs"] *= scale
    liq["ngl_smt_2024"]["qc1Ncs"] = [q * scale for q in liq["ngl_smt_2024"]["qc1Ncs"]]
    return liq


PROFILES = [
    {"name": "Base", "liquefaction_models": _profile_models(1.0)},
    {"name": "loose", "liquefaction_models": _profile_models(0.8)},
    {"name": "dense", "liquefaction_models": _profile_models(1.3)},
]


def _with_profiles(config, disaggregation=False):
    config = copy.deepcopy(config)
    del config["liquefaction_models"]
    config["profiles"] = copy.deepcopy(PROFILES)
    if disaggregation:
        config["output"]["plha"]["disaggregation"] = BINS
    return config


def _single(config, profile, disaggregation=False):
    config = copy.deepcopy(config)
    config["liquefaction_models"] = copy.deepcopy(profile["liquefaction_models"])
    if disaggregation:
        config["output"]["plha"]["disaggregation"] = BINS
    return config


def test_config_dict_same_as_config_file(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(CONFIG))
    config = copy.deepcopy(CONFIG)
    from_file = plha.get_hazard(str(path))
    from_dict = plha.get_hazard(config)
    assert from_dict["output"] == from_file["output"]
    assert from_dict["input"] == from_file["input"]
    # the config dict is not changed (the weights are normalized in the output copy)
    assert config == CONFIG


def test_default_outputfile_needs_a_file_name():
    config = copy.deepcopy(CONFIG)
    config["output"]["outputfile"] = "default"
    with pytest.raises(ValueError, match="outputfile"):
        plha.get_hazard(config)


def test_loaded_models_give_identical_results():
    models = plha.load_source_models(CONFIG)
    for site in (SITE, SECOND_SITE):
        config = copy.deepcopy(CONFIG)
        config["site"] = site
        expected = plha.get_hazard(copy.deepcopy(config))
        for _ in range(2):
            assert plha.get_hazard(config, models=models)["output"] == expected["output"]


def test_loaded_models_are_read_only():
    models = plha.load_source_models(CONFIG)
    path = plha._source_model_path("fault_source_models", "ucerf3_fm31")
    ruptures = models.load(path, "ruptures.npz")
    assert "m" in ruptures.files
    with pytest.raises(ValueError):
        ruptures["m"][0] = 0.0


def test_source_data_with_models_identical():
    models = plha.SourceModels()
    for source_type, model in [
        ("fault_source_models", "ucerf3_fm31"),
        ("point_source_models", "ucerf3_fm31_grid_unassociated"),
    ]:
        p_xyz = plha.geometry.point_to_xyz(np.array([36.8, -121.8, 0.0]))
        args = (source_type, model, p_xyz, 200.0, 5.0, ["bssa14", "ask14"])
        expected = plha.get_source_data(*args)
        for _ in range(2):
            out = plha.get_source_data(*args, models=models)
            assert len(out) == len(expected)
            for a, b in zip(out, expected):
                assert np.array_equal(a, b)


@pytest.mark.parametrize("disaggregation", [False, True])
def test_profiles_identical_to_single_profile_runs(disaggregation):
    out = plha.get_hazard(_with_profiles(CONFIG, disaggregation))
    assert "plha" not in out["output"]
    assert [p["name"] for p in out["output"]["profiles"]] == ["base", "loose", "dense"]
    curves = []
    for profile, result in zip(PROFILES, out["output"]["profiles"]):
        single = plha.get_hazard(_single(CONFIG, profile, disaggregation))
        assert result["plha"] == single["output"]["plha"]
        assert out["output"]["psha"] == single["output"]["psha"]
        curves.append(result["plha"]["annual_rate_of_nonexceedance"])
    assert curves[0] != curves[1] != curves[2]


def test_profiles_with_loaded_models():
    models = plha.load_source_models(_with_profiles(CONFIG))
    config = _with_profiles(CONFIG)
    config["site"] = SECOND_SITE
    out = plha.get_hazard(config, models=models)
    for profile, result in zip(PROFILES, out["output"]["profiles"]):
        single = _single(CONFIG, profile)
        single["site"] = SECOND_SITE
        assert result["plha"] == plha.get_hazard(single)["output"]["plha"]


def test_profile_default_names():
    config = _with_profiles(CONFIG)
    for profile in config["profiles"]:
        del profile["name"]
    out = plha.get_hazard(config)
    assert [p["name"] for p in out["output"]["profiles"]] == ["0", "1", "2"]


def test_profiles_and_liquefaction_models_rejected(capsys):
    config = _with_profiles(CONFIG)
    config["liquefaction_models"] = LIQUEFACTION
    assert plha.get_hazard(config) is None
    assert "profiles" in capsys.readouterr().out


def test_profile_without_liquefaction_models_rejected(capsys):
    config = _with_profiles(CONFIG)
    config["profiles"].append({"name": "empty"})
    assert plha.get_hazard(config) is None
    assert "Config File Error" in capsys.readouterr().out

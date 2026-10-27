"""ucla_plha.batch: many sites with a pool of worker processes, results in chunk files."""

import copy
import csv
import json

import numpy as np
import pytest

from ucla_plha import batch, plha

LIQUEFACTION = {
    "boulanger_idriss_2016": {
        "depth": 2.2, "qc1ncs": 71.8, "sigmav": 40.0, "sigmavp": 28.0, "pa": 101.325,
        "weight": 0.5,
    },
    "ngl_smt_2024": {
        "ztop": [0.0, 1.0, 5.0], "zbot": [1.0, 5.0, 9.0], "qc1ncs": [50.0, 150.0, 250.0],
        "ic": [1.5, 2.0, 2.6], "sigmav": [20.0, 60.0, 120.0], "sigmavp": [15.0, 45.0, 90.0],
        "ksat": [0.0, 1.0, 1.0], "pa": 101.325, "weight": 0.5,
    },
}


def _profiles(scales):
    out = []
    for i, scale in enumerate(scales):
        liq = copy.deepcopy(LIQUEFACTION)
        liq["boulanger_idriss_2016"]["qc1ncs"] *= scale
        liq["ngl_smt_2024"]["qc1ncs"] = [q * scale for q in liq["ngl_smt_2024"]["qc1ncs"]]
        out.append({"name": f"p{i}", "liquefaction_models": liq})
    return out


BASE = {
    "site": {"latitude": 0.0, "longitude": 0.0, "elevation": 0.0, "vs30": 300.0},
    "constraints": {"dist_cutoff": 60, "m_min": 7.0},
    "source_models": {
        "fault_source_models": {"ucerf3_fm31": {"weight": 1.0}},
        "point_source_models": {"ucerf3_fm31_grid_unassociated": {"weight": 1.0}},
    },
    "ground_motion_models": {"bssa14": {"weight": 0.5}, "ask14": {"weight": 0.5}},
    "profiles": _profiles([1.0, 0.8]),
    "output": {"psha": {"pga": [0.05, 0.2, 0.8]}, "plha": {"fsl": [0.5, 1.0, 2.0]}},
}

SITES = [
    {"site_id": "monterey", "latitude": 36.80547, "longitude": -121.786074, "vs30": 230.0},
    {"site_id": "la", "latitude": 34.05, "longitude": -118.25, "vs30": 300.0,
     "profiles_file": "la_profiles.json"},
    {"site_id": "sf", "latitude": 37.77, "longitude": -122.42, "vs30": 250.0},
]
LA_PROFILES = _profiles([1.2, 0.9, 0.7])


def _write_inputs(tmp_path):
    with open(tmp_path / "sites.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, ["site_id", "latitude", "longitude", "vs30", "profiles_file"])
        writer.writeheader()
        for row in SITES:
            writer.writerow(row)
    (tmp_path / "la_profiles.json").write_text(json.dumps(LA_PROFILES))
    (tmp_path / "base.json").write_text(json.dumps(BASE))


def _expected(site):
    config = copy.deepcopy(BASE)
    config["site"] = {"latitude": site["latitude"], "longitude": site["longitude"],
                      "elevation": 0.0, "vs30": site["vs30"]}
    profiles = LA_PROFILES if "profiles_file" in site else BASE["profiles"]
    psha = None
    plha_curves = {}
    for profile in profiles:
        # each profile alone, as a single-profile config
        single = copy.deepcopy(config)
        del single["profiles"]
        single["liquefaction_models"] = copy.deepcopy(profile["liquefaction_models"])
        out = plha.get_hazard(single)
        psha = np.array(out["output"]["psha"]["annual_rate_of_exceedance"])
        plha_curves[profile["name"]] = np.array(out["output"]["plha"]["annual_rate_of_nonexceedance"])
    return psha, plha_curves


@pytest.mark.parametrize("processes", [1, 2])
def test_batch_same_as_single_runs(tmp_path, processes):
    _write_inputs(tmp_path)
    out_dir = tmp_path / "out"
    written = batch.run(str(tmp_path / "sites.csv"), str(tmp_path / "base.json"), str(out_dir),
                        processes=processes, chunk_size=2, verbose=False)
    assert len(written) == 2
    results = batch.read_results(str(out_dir))
    assert sorted(results) == sorted(s["site_id"] for s in SITES)
    for site in SITES:
        psha, plha_curves = _expected(site)
        got = results[site["site_id"]]
        assert np.array_equal(got["psha"], psha)
        assert sorted(got["plha"]) == sorted(plha_curves)
        for name, curve in plha_curves.items():
            assert np.array_equal(got["plha"][name], curve)
    # chunks that exist are not run again
    assert batch.run(str(tmp_path / "sites.csv"), str(tmp_path / "base.json"), str(out_dir),
                     processes=processes, chunk_size=2, verbose=False) == []


def test_batch_start_and_max_sites(tmp_path):
    _write_inputs(tmp_path)
    out_dir = tmp_path / "out"
    batch.run(str(tmp_path / "sites.csv"), str(tmp_path / "base.json"), str(out_dir),
              processes=1, start=1, max_sites=1, verbose=False)
    assert sorted(batch.read_results(str(out_dir))) == ["la"]
    assert (out_dir / "chunk_000000001.npz").exists()


def test_batch_hdf5(tmp_path):
    pytest.importorskip("h5py")
    _write_inputs(tmp_path)
    out_dir = tmp_path / "out"
    batch.run(str(tmp_path / "sites.csv"), str(tmp_path / "base.json"), str(out_dir),
              processes=1, max_sites=1, file_format="hdf5", verbose=False)
    results = batch.read_results(str(out_dir))
    psha, plha_curves = _expected(SITES[0])
    assert np.array_equal(results["monterey"]["psha"], psha)
    for name, curve in plha_curves.items():
        assert np.array_equal(results["monterey"]["plha"][name], curve)

"""Config and helpers shared by the tests of the liquefaction models conditioned on PGA
(test_ngl_shared_pga.py, test_conditional_pga_models.py)."""

import json
import os
import tempfile
from importlib.resources import files

import jsonschema
import numpy as np

from ucla_plha import plha

NGL = dict(
    ztop=np.array([0.0, 1.0, 5.0]),
    zbot=np.array([1.0, 5.0, 9.0]),
    qc1Ncs=np.array([50.0, 150.0, 250.0]),
    Ic=np.array([1.5, 2.0, 2.6]),
    sigmav=np.array([20.0, 60.0, 120.0]),
    sigmavp=np.array([15.0, 45.0, 90.0]),
    Ksat=np.array([0.0, 1.0, 1.0]),
)

CONFIG = {
    "site": {"latitude": 36.80547, "longitude": -121.786074, "elevation": 0, "vs30": 230.0},
    "constraints": {"dist_cutoff": 60, "m_min": 7.0},
    "source_models": {"fault_source_models": {"ucerf3_fm31": {"weight": 1.0}}},
    "ground_motion_models": {"bssa14": {"weight": 1.0}},
    "liquefaction_models": {
        "ngl_smt_2024": {
            **{k: v.tolist() for k, v in NGL.items()}, "pa": 101.325, "weight": 1.0,
        },
    },
    "output": {"plha": {"fsl": [0.5, 1.0, 2.0]}},
}
BINS = {
    "magnitude_bin_edges": [6.5, 7.0, 7.5, 8.5],
    "distance_bin_edges": [0.0, 30.0, 60.0],
    "epsilon_bin_edges": [-500.0, -1.0, 0.0, 1.0, 500.0],
}


def run_hazard(config):
    """plha.get_hazard of a config dict (written to a temporary config file)."""
    fd, path = tempfile.mkstemp(suffix=".json")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(config, f)
        return plha.get_hazard(path)
    finally:
        os.remove(path)


def validate(config):
    """Validates a config dict against the schema, as plha.get_hazard does (lower case)."""
    schema = json.loads(files("ucla_plha").joinpath("ucla_plha_schema.json").read_text().lower())
    jsonschema.validate(json.loads(json.dumps(config).lower()), schema)

"""
Generate gmm_reference.json, the reference PGA values used by tests/test_gmm_reference.py.

The reference values come from pygmm (https://github.com/arkottke/pygmm), an independent
implementation of the NGA-West2 ground motion models. Storing them in a JSON file means the
unit tests do not depend on pygmm. Regenerate the file only when the scenarios change:

    pip install pygmm
    python tests/data/make_gmm_reference.py

The main scenarios use the default basin depths (z1p0 and z2p5 not specified) and provide
reference values for all five models (ASK14, BSSA14, CB14, CY14, and Idriss 2014). A second set of scenarios specifies z1p0 and provides
reference values for ASK14 only. User-specified basin depths are not yet checked for CB14 and
CY14 (for PGA, the CY14 basin term is zero).
"""

import itertools
import json
import logging
import warnings
from pathlib import Path

import numpy as np
import pygmm

logging.disable(logging.WARNING)
# pygmm warns when scenarios are outside a model's recommended range (e.g., Vs30 < 450 m/s for
# Idriss 2014). The reference values are still valid tests of the model equations.
warnings.simplefilter("ignore", UserWarning)

MECHANISM = {1: "RS", 2: "NS", 3: "SS"}
DIP = {1: 40.0, 2: 55.0, 3: 90.0}
ZTOR = {5.0: 6.0, 6.0: 3.0, 7.0: 0.0, 7.8: 0.0}
SEISMOGENIC_DEPTH = 15.0
MODELS = {
    "bssa14": pygmm.BooreStewartSeyhanAtkinson2014,
    "ask14": pygmm.AbrahamsonSilvaKamai2014,
    "cb14": pygmm.CampbellBozorgnia2014,
    "cy14": pygmm.ChiouYoungs2014,
    "idriss14": pygmm.Idriss2014,
}


def rupture_width(fault_type, m, ztor, dip):
    """Down-dip width from Wells and Coppersmith (1994), limited by the seismogenic depth"""
    if fault_type == 3:
        w = 10 ** (-1.01 + 0.32 * m)
    elif fault_type == 1:
        w = 10 ** (-1.61 + 0.41 * m)
    else:
        w = 10 ** (-1.14 + 0.35 * m)
    return min(w, (SEISMOGENIC_DEPTH - ztor) / np.sin(np.radians(dip)))


def distances(rx, ztor, zbor, width, dip):
    """Rjb, Rrup, and Rx1 for a site on a line perpendicular to the strike of a planar rupture of infinite length"""
    d = np.radians(dip)
    wh = width * np.cos(d)
    if rx < 0:
        rjb = -rx
    elif rx <= wh:
        rjb = 0.0
    else:
        rjb = rx - wh
    if dip == 90.0 or rx < ztor * np.tan(d):
        rrup = np.hypot(rx, ztor)
    elif rx <= ztor * np.tan(d) + width / np.cos(d):
        rrup = rx * np.sin(d) + ztor * np.cos(d)
    else:
        rrup = np.hypot(rx - wh, zbor)
    return rjb, rrup, rx - wh


def scenario(fault_type, m, rx, vs30, measured_vs30, models, z1p0=None):
    """Distances, rupture geometry, and pygmm reference values for one scenario"""
    dip = DIP[fault_type]
    ztor = ZTOR[m]
    width = rupture_width(fault_type, m, ztor, dip)
    zbor = ztor + width * np.sin(np.radians(dip))
    rjb, rrup, rx1 = distances(rx, ztor, zbor, width, dip)
    s = pygmm.model.Scenario(
        mag=m,
        dist_jb=rjb,
        dist_rup=rrup,
        dist_x=rx,
        dist_y0=0.0,
        dip=dip,
        v_s30=vs30,
        mechanism=MECHANISM[fault_type],
        depth_tor=ztor,
        depth_bor=zbor,
        depth_1_0=z1p0,
        width=width,
        on_hanging_wall=bool(rx >= 0),
        region="california",
        vs_source="measured" if measured_vs30 else "inferred",
    )
    expected = {}
    for name in models:
        out = MODELS[name](s)
        expected[name] = {"mu": float(np.log(out.pga)), "sigma": float(out.ln_std_pga)}
    row = {
        "fault_type": fault_type,
        "m": m,
        "vs30": vs30,
        "measured_vs30": measured_vs30,
        "dip": dip,
        "ztor": ztor,
        "zbor": float(zbor),
        "rjb": float(rjb),
        "rrup": float(rrup),
        "rx": rx,
        "rx1": float(rx1),
        "ry0": 0.0,
        "expected": expected,
    }
    if z1p0 is not None:
        row["z1p0"] = z1p0
    return row


scenarios = [
    scenario(fault_type, m, rx, vs30, measured_vs30, MODELS)
    for fault_type, m, rx, vs30, measured_vs30 in itertools.product(
        [1, 2, 3],
        [5.0, 6.0, 7.0, 7.8],
        [-30.0, -5.0, 2.0, 10.0, 40.0, 150.0],
        [250.0, 450.0, 760.0],
        [False, True],
    )
]

# Vs30 values span all four soil depth bins of ASK14 Equation 17, including values between bin
# centers, and z1p0 values span shallower and deeper than the default Z1 for each Vs30.
ask14_basin_scenarios = [
    scenario(fault_type, m, rx, vs30, False, ["ask14"], z1p0=z1p0)
    for (fault_type, m, rx), vs30, z1p0 in itertools.product(
        [(3, 5.0, 5.0), (3, 7.0, 10.0), (1, 7.0, 10.0), (2, 7.8, 50.0)],
        [180.0, 220.0, 270.0, 350.0, 450.0, 550.0, 760.0, 1000.0],
        [0.03, 0.1, 0.3, 0.6, 1.2, 2.0],
    )
]

output = {
    "source": f"pygmm {pygmm.__version__}",
    "scenarios": scenarios,
    "ask14_basin_scenarios": ask14_basin_scenarios,
}
path = Path(__file__).with_name("gmm_reference.json")
path.write_text(json.dumps(output, indent=1))
print(
    f"Wrote {len(scenarios)} scenarios and {len(ask14_basin_scenarios)} ASK14 basin scenarios to {path}"
)

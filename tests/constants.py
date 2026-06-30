import numpy as np

from ucla_plha.geometry import geometry

SITE_LATLONELEV = np.asarray([36.80547, -121.786074, 0.0])
P_XYZ = geometry.point_to_xyz(SITE_LATLONELEV)

_M = np.array(
    [3.5, 4.0, 4.3, 4.5, 4.7, 4.97, 5.0, 5.3, 5.5,
     5.85, 6.0, 6.5, 6.6, 6.75, 7.0, 7.5, 8.0]
)
_N = len(_M)
GMM_INPUTS = {
    "m": _M,
    "fault_type": np.array([1, 2, 3] * 6, dtype=int)[:_N],
    "rjb": np.linspace(1.0, 150.0, _N),
    "rrup": np.linspace(1.0, 150.0, _N) + 3.0,
    "rx": np.linspace(-50.0, 100.0, _N),
    "rx1": np.linspace(-50.0, 100.0, _N) + 10.0,
    "ry0": np.linspace(0.0, 90.0, _N),
    "dip": np.array([30.0, 45.0, 60.0, 90.0] * 5, dtype=float)[:_N],
    "ztor": np.linspace(0.0, 18.0, _N),
    "zbor": np.linspace(0.0, 18.0, _N) + 5.0,
}

PGA_INPUTS = (
    _M,
    np.full(_N, np.log(0.2)),
    np.full(_N, 0.6),
    np.array([0.5, 1.0, 2.0]),
)

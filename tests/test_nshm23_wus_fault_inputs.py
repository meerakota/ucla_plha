"""Ground motion model inputs of NSHM23 WUS fault system ruptures, compared with nshmp-lib.

The expected values are the nshmp-lib (commit 44728a7d) inputs of the ruptures at a Los Angeles
site (-118.25, 34.05), from a Python port of SystemRuptureSet.InputGenerator,
DefaultGriddedSurface (1 km grid, aseismicity removes the top of the fault plane), and
Distance.compute. The rRup values agree with the USGS disaggregation service
(conus-2023/dynamic/disagg/-118.25/34.05/760/2475, PGA) to 1 m: Puente Hills (Los Angeles) (1)
4.672 km, Elysian Park (upper) (1) 5.395 km, Compton (2) 14.482 km, Raymond (2) 8.757 km,
Newport - Inglewood (8) 11.702 km.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

from ucla_plha import plha
from ucla_plha.geometry import geometry

LA = (34.05, -118.25)
GMMS = ["ask14", "bssa14", "cb14", "cy14"]

# rupture index: (closest section, rJB, rRup, rX) from the nshmp-lib port
NSHMP_INPUTS = {
    319060: ("Puente Hills (Los Angeles) (1)", 0.000, 4.672, 3.505),
    319192: ("Puente Hills (Los Angeles) (1)", 0.000, 4.672, 3.505),
    128059: ("Puente Hills (Los Angeles) (1)", 0.000, 4.672, 3.505),
    128063: ("Elysian Park (upper) (1)", 3.387, 5.395, -3.387),
    90014: ("Compton (2)", 0.425, 14.482, 24.848),
    320640: ("Raymond (2)", 8.617, 8.757, -8.246),
    461164: ("Newport - Inglewood (8)", 0.425, 11.702, 11.604),
}


def _la_inputs(rupture_rx="closest_section", **kwargs):
    p_xyz = geometry.point_to_xyz(np.array([LA[0], LA[1], 0.0]))
    arrays, extras = plha.get_source_data(
        "fault_source_models", "nshm23_wus", p_xyz, 300.0, None, GMMS,
        extras=True, rupture_rx=rupture_rx, **kwargs,
    )
    pos = {int(j): k for k, j in enumerate(extras["index"])}
    return arrays, pos


def test_closest_segment_selection_and_ties():
    # three ruptures: segments (0, 1), (2,), (3, 4, 5)
    segment_index = np.array([0, 1, 2, 5, 4, 3], dtype=np.int32)
    rrup_seg = np.array([3.0, 1.0, 7.0, 2.0, 2.0, 9.0])
    boundaries = np.array([0, 2, 3])
    rrup = np.minimum.reduceat(rrup_seg, boundaries)
    closest = plha._closest_segment(rrup_seg, rrup, segment_index, boundaries)
    # ties (segments 5 and 4) go to the lowest segment index, as in nshmp-lib
    np.testing.assert_array_equal(closest, [1, 2, 4])


@pytest.mark.parametrize("rupture", sorted(NSHMP_INPUTS))
def test_la_rupture_inputs_match_nshmp(rupture):
    (m, fault_type, rate, rjb, rrup, rx, rx1, ry0, dip, ztor, zbor), pos = _la_inputs()
    _, rjb_n, rrup_n, rx_n = NSHMP_INPUTS[rupture]
    k = pos[rupture]
    # the nshmp-lib gridded section surfaces (source_info "fault_distances": "nshmp_grid")
    assert rrup[k] == pytest.approx(rrup_n, abs=0.002)
    assert rx[k] == pytest.approx(rx_n, abs=0.002)
    assert rjb[k] == pytest.approx(rjb_n, abs=0.002)
    # rx1 is the rx1 of the same (closest) segment (equal to rx for a vertical segment)
    assert rx1[k] <= rx[k]


@pytest.mark.parametrize("rupture", sorted(NSHMP_INPUTS))
def test_la_rupture_inputs_with_triangles(rupture):
    # distances to the triangles of the planes: rrup within 0.06 km; nshmp-lib rJB is the
    # distance to the nearest 1 km grid point unless the site is inside the surface perimeter
    # (built from the trace, without the aseismic shift), so it can differ by a few hundred meters
    (m, fault_type, rate, rjb, rrup, rx, rx1, ry0, dip, ztor, zbor), pos = _la_inputs(
        fault_distances="triangles")
    _, rjb_n, rrup_n, rx_n = NSHMP_INPUTS[rupture]
    k = pos[rupture]
    assert rrup[k] == pytest.approx(rrup_n, abs=0.06)
    assert rx[k] == pytest.approx(rx_n, abs=0.06)
    assert rjb[k] <= rjb_n + 0.05
    assert rjb[k] >= rjb_n - 0.5


def test_minimum_rx_convention_is_available():
    # the ucla_plha 2.x conventions: minimum over the segments of the rupture and over the
    # planar pieces of each segment
    legacy = dict(fault_distances="triangles", section_rx="minimum")
    arrays, pos = _la_inputs("minimum", **legacy)
    rx_min = arrays[5]
    arrays_c, pos_c = _la_inputs(**legacy)
    rx_c = arrays_c[5]
    # Puente Hills rupture with sections of the Elysian Park and Raymond faults: the site is on
    # the hanging wall of the closest section, but the minimum rx is on a footwall section
    k = pos[319192]
    assert rx_min[k] < -5.0
    assert rx_c[pos_c[319192]] > 3.0
    assert np.all(rx_min <= rx_c + 1e-9)
    # the minimum over the segments with the nshmp-lib rx of each segment
    arrays_n, pos_n = _la_inputs("minimum")
    assert arrays_n[5][pos_n[319192]] < -4.0
    with pytest.raises(ValueError):
        _la_inputs("closest")


def test_la_fault_system_hazard_matches_usgs():
    pytest.importorskip("pygmm")
    from ucla_plha import pygmm_gmms

    try:
        pygmm_gmms.resolve("ask_14_basin")
    except Exception:
        pytest.skip("pygmm does not provide the nshmp-lib NGA-West2 models")
    validation = Path(__file__).resolve().parents[1] / "validation"
    cache = validation / "usgs_nshm_cache"
    if not (cache / "conus2023_-118.250_34.050_760.json").exists():
        pytest.skip("cached USGS response not available")
    sys.path.insert(0, str(validation))
    try:
        import compare_usgs_nshm as cu
    finally:
        sys.path.remove(str(validation))
    models = {k: v for k, v in cu.available_models().items() if k == "nshm23_wus"}
    result = cu.compare_site("los_angeles", -118.25, 34.05, 760.0, "nshmp", models, str(cache), "source_info")
    row = [r for r in result["components"] if r["component"] == "FaultSystem"][0]
    for rp in (475, 2475):
        ratio = row[f"pga_{rp}yr_ucla_plha"] / row[f"pga_{rp}yr_usgs"]
        assert ratio == pytest.approx(1.0, abs=0.03), rp
    xs = np.array(row["xs"])
    i = int(np.argmin(np.abs(np.log(xs / 1.0))))
    # rate ratio at 1 g: 0.64 with the minimum rx convention and the earlier geometry
    assert 0.92 < row["ratio"][i] < 1.08

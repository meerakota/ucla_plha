"""nshmp-lib distance conventions (model.Distance.compute) for the NSHM23 fault source models.

The expected values are the nshmp-lib (commit 44728a7d) rJB, rRup, and rX of ruptures from a
Python port of DefaultGriddedSurface (WUS sections with aseismicity, CEUS faults, 1 km),
ApproxGriddedSurface (Cascadia interface, 5 km), GriddedSubsetSurface (floating ruptures),
SystemRuptureSet.InputGenerator, and Distance.compute, which reproduces the rRup of the USGS
disaggregation service to 1 m. The grid points are stored as float32 (grid.npz), so the
distances agree to a few meters.
"""
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

from ucla_plha import plha
from ucla_plha.geometry import geometry

SITES = {
    "los_angeles": (-118.25, 34.05),
    "portland": (-122.68, 45.52),
    "seattle": (-122.30, 47.60),
    "coastal_oregon": (-124.05, 44.63),
    "new_madrid": (-89.60, 36.60),
}

# (site, model): {rupture index: (rJB, rRup, rX)} from the nshmp-lib port
NSHMP_INPUTS = {
    ("coastal_oregon", "nshm23_cascadia_interface"): {
        27: (0.0000, 25.3041, 114.0079),
        1064: (0.0000, 25.3984, 112.2366),
        3800: (41.5150, 44.3027, 114.0231),
        819: (432.2676, 432.2965, -21.5465),
    },
    ("seattle", "nshm23_cascadia_interface_cluster"): {
        41: (84.3154, 91.1123, 274.6578),
        14: (212.0308, 214.7240, 262.3074),
        224: (715.1248, 715.2308, 200.9467),
    },
    ("new_madrid", "nshm23_ceus_fault"): {
        1060: (2.5702, 10.3250, -2.5284),
        1057: (0.0000, 10.9157, 4.3568),
        1037: (10.2932, 14.3510, -10.2929),
    },
    ("new_madrid", "nshm23_ceus_fault_cluster"): {
        227: (0.7361, 10.0271, -0.6247),
        174: (0.0000, 10.9157, 4.3568),
        120: (64.0404, 64.8165, 63.6093),
    },
    # Portland Hills: the site is 0.8 km up dip of the top of the plane after the aseismic
    # reduction but inside the nshmp-lib perimeter (built from the trace), so rJB = 0
    ("portland", "nshm23_wus"): {
        13741: (0.0000, 1.7052, -0.7814),
        318572: (0.0000, 1.7052, -0.7814),
        483564: (206.2541, 206.2595, -205.4126),
    },
    ("los_angeles", "nshm23_wus"): {
        124706: (0.0000, 4.6715, 3.5048),
        232622: (0.0000, 6.3090, -0.7830),
        400832: (19.8573, 19.9080, -19.8573),
        357530: (55.5303, 55.5458, 55.4099),
        386024: (297.2852, 297.3676, -63.7167),
    },
}


def _inputs(site, model, **kwargs):
    lon, lat = SITES[site]
    p_xyz = geometry.point_to_xyz(np.array([lat, lon, 0.0]))
    # cb14 uses rjb, rrup, and rx
    arrays, extras = plha.get_source_data(
        "fault_source_models", model, p_xyz, 1000.0, None, ["cb14"], extras=True, **kwargs
    )
    pos = {int(j): k for k, j in enumerate(extras["index"])}
    return arrays, pos


@pytest.mark.parametrize("site,model", sorted(NSHMP_INPUTS))
def test_rupture_distances_match_nshmp(site, model):
    arrays, pos = _inputs(site, model)
    rjb, rrup, rx = arrays[3], arrays[4], arrays[5]
    for rupture, expected in NSHMP_INPUTS[(site, model)].items():
        k = pos[rupture]
        np.testing.assert_allclose([rjb[k], rrup[k], rx[k]], expected, atol=0.005, err_msg=str(rupture))


def test_source_info_settings():
    for model in ["nshm23_wus", "nshm23_cascadia_interface", "nshm23_cascadia_interface_cluster",
                  "nshm23_ceus_fault", "nshm23_ceus_fault_cluster"]:
        assert plha.get_source_info("fault_source_models", model)["fault_distances"] == "nshmp_grid"
    assert plha.get_source_info("fault_source_models", "nshm23_wus")["section_rx"] == "extended_trace"
    for model in ["ucerf3_fm31", "ucerf3_fm32"]:
        info = plha.get_source_info("fault_source_models", model)
        assert info.get("fault_distances", "triangles") == "triangles"
        assert info.get("section_rx", "minimum") == "minimum"


def test_triangles_option():
    # Cascadia interface triangles are flat strips 100-200 km wide: rJB of a few hundred
    # meters above the interface instead of 0, and rRup within ~0.5 km of nshmp-lib
    arrays, pos = _inputs("coastal_oregon", "nshm23_cascadia_interface", fault_distances="triangles")
    rjb, rrup = arrays[3], arrays[4]
    k = pos[27]
    assert 0.0 < rjb[k] < 0.6
    assert rrup[k] == pytest.approx(25.3041, abs=0.6)
    with pytest.raises(ValueError):
        _inputs("coastal_oregon", "nshm23_cascadia_interface", fault_distances="grid")
    with pytest.raises(ValueError):
        _inputs("portland", "ucerf3_fm31", fault_distances="nshmp_grid")
    with pytest.raises(ValueError):
        _inputs("portland", "nshm23_wus", section_rx="closest")


def test_ucerf3_defaults_unchanged():
    a, _ = _inputs("los_angeles", "ucerf3_fm31")
    b, _ = _inputs("los_angeles", "ucerf3_fm31", fault_distances="triangles", section_rx="minimum")
    for x, y in zip(a, b):
        np.testing.assert_array_equal(x, y)


def test_section_rx_extended_trace_wus():
    # Portland Hills: the minimum over the planar pieces of the section gives -1.75 km (a
    # piece whose line passes farther from the site); the extended trace gives -0.78 km
    a, pos = _inputs("portland", "nshm23_wus", fault_distances="triangles", section_rx="minimum")
    b, pos_b = _inputs("portland", "nshm23_wus", fault_distances="triangles")
    assert a[5][pos[318572]] == pytest.approx(-1.750, abs=0.01)
    assert b[5][pos_b[318572]] == pytest.approx(-0.7814, abs=0.01)


def test_extended_trace_rx():
    # trace from (0, 0) north to (0, 0.1) degrees: strike 0, dip direction east
    lon1, lat1 = np.array([0.0, 0.0]), np.array([0.0, 0.05])
    lon2, lat2 = np.array([0.0, 0.0]), np.array([0.05, 0.1])
    b = np.array([0])
    km = geometry.NSHMP_EARTH_RADIUS * np.pi / 180.0
    east = geometry.extended_trace_rx(lon1, lat1, lon2, lat2, b, 0.1, 0.05)
    west = geometry.extended_trace_rx(lon1, lat1, lon2, lat2, b, -0.1, 0.05)
    assert east[0] == pytest.approx(0.1 * km, rel=1e-3)
    assert west[0] == pytest.approx(-0.1 * km, rel=1e-3)
    # beyond the north end: the distance to the extension of the trace, not to its end point
    north = geometry.extended_trace_rx(lon1, lat1, lon2, lat2, b, -0.1, 0.5)
    assert north[0] == pytest.approx(-0.1 * km * np.cos(np.radians(0.25)), rel=2e-3)
    # on the trace: zero, counted as hanging wall
    on = geometry.extended_trace_rx(lon1, lat1, lon2, lat2, b, 0.0, 0.07)
    assert on[0] == 0.0
    # two polylines at once (the second reversed: dip direction west)
    two = geometry.extended_trace_rx(
        np.r_[lon1, lon2[::-1]], np.r_[lat1, lat2[::-1]], np.r_[lon2, lon1[::-1]],
        np.r_[lat2, lat1[::-1]], np.array([0, 2]), 0.1, 0.05)
    assert two[0] == pytest.approx(east[0]) and two[1] == pytest.approx(-east[0])


def test_gridded_rupture_distances_simple_grid():
    # one vertical surface along the equator (2 rows x 3 columns, 1 km apart) and a rupture on
    # its first two columns
    deg = 1.0 / (geometry.NSHMP_EARTH_RADIUS * np.pi / 180.0)
    lon = np.repeat(np.array([0.0, deg, 2 * deg]), 2)
    lat = np.zeros(6)
    depth = np.tile([0.0, 1.0], 3)
    grid = dict(lon=lon, lat=lat, depth=depth, surface_offset=np.array([0]),
                surface_rows=np.array([2]), surface_cols=np.array([3]),
                surface_dip=np.array([90.0]), surface_spacing=np.array([1.0]),
                rupture_surface=np.array([0, 0]), rupture_row0=np.array([0, 1]),
                rupture_rows=np.array([2, 1]), rupture_col0=np.array([0, 1]),
                rupture_cols=np.array([2, 2]))
    rjb, rrup, rx = geometry.gridded_rupture_distances(grid, 3 * deg, 0.0)
    # the site is 3 km north of the first column
    np.testing.assert_allclose(rjb, [3.0, np.hypot(3.0, 1.0)], rtol=1e-6)
    # dip > 89: top row only (the second rupture's top row is at 1 km depth)
    np.testing.assert_allclose(rrup, [3.0, np.hypot(np.hypot(3.0, 1.0), 1.0)], rtol=1e-6)
    # strike east (azimuth 90), dip direction south: the site (north) is on the footwall
    assert rx[0] == pytest.approx(-3.0, rel=1e-6)


def _compare(site, models):
    pytest.importorskip("pygmm")
    from ucla_plha import pygmm_gmms

    try:
        pygmm_gmms.resolve("ask_14_basin")
    except Exception:
        pytest.skip("pygmm does not provide the nshmp-lib models")
    validation = Path(__file__).resolve().parents[1] / "validation"
    cache = validation / "usgs_nshm_cache"
    sys.path.insert(0, str(validation))
    try:
        import compare_usgs_nshm as cu
    finally:
        sys.path.remove(str(validation))
    lon, lat = cu.SITES[site]
    if not (cache / f"conus2023_{lon:.3f}_{lat:.3f}_760.json").exists():
        pytest.skip("cached USGS response not available")
    available = cu.available_models()
    result = cu.compare_site(site, lon, lat, 760.0, "nshmp", {k: available[k] for k in models},
                             str(cache), "source_info")
    return {row["component"]: row for row in result["components"]}


def _ratio_at(row, g):
    xs = np.array(row["xs"])
    return row["ratio"][int(np.argmin(np.abs(np.log(xs / g))))]


@pytest.mark.parametrize("site,models,component,levels", [
    # rX of the Portland Hills sections: 0.94 at 1 g with the minimum over the planar pieces
    ("portland", ["nshm23_wus"], "FaultSystem", [0.4, 1.0]),
    # Cascadia interface: 1.014 at 1 g and +0.6% at 475 years with the triangles
    ("seattle", ["nshm23_cascadia_interface"], "Interface", [0.1, 0.4, 1.0]),
    ("seattle", ["nshm23_cascadia_interface_cluster"], "FaultCluster", [0.1, 0.4, 1.0]),
    ("coastal_oregon", ["nshm23_cascadia_interface"], "Interface", [0.1, 0.4, 1.0]),
    # New Madrid cluster 900 km away: +0.9% at 2475 years with the triangles
    ("charleston", ["nshm23_ceus_fault_cluster"], "FaultCluster", [0.01, 0.1]),
    ("new_madrid", ["nshm23_ceus_fault", "nshm23_ceus_fault_cluster"], "FaultCluster", [0.1, 0.4, 1.0]),
])
def test_usgs_comparison_with_nshmp_distances(site, models, component, levels):
    row = _compare(site, models)[component]
    for g in levels:
        assert _ratio_at(row, g) == pytest.approx(1.0, abs=0.006), g
    for rp in (475, 2475):
        u, c = row[f"pga_{rp}yr_usgs"], row[f"pga_{rp}yr_ucla_plha"]
        if np.isfinite(u) and np.isfinite(c):
            assert c / u == pytest.approx(1.0, abs=0.002), rp


def test_config_constraints_are_valid():
    import jsonschema
    from importlib.resources import files

    schema = json.loads(files("ucla_plha").joinpath("ucla_plha_schema.json").read_text())
    props = schema["properties"]["constraints"]["properties"]
    assert set(props["fault_distances"]["enum"]) == {"source_info", "triangles"}
    assert set(props["section_rx"]["enum"]) == {"source_info", "minimum", "extended_trace"}

"""Aseismicity in the UCERF3 fault source models (utilities/convert_ucerf3_fault_source_models.py).

UCERF3 (OpenSHA FaultSection.getReducedAveUpperDepth, StirlingGriddedSurface) moves the upper
edge of each subsection down-dip from its trace to the reduced upper depth
upper + aseismicity * (lower - upper) and keeps the lower edge. The expected rRup values are those
of the USGS disaggregation service for the 2018 NSHM, which uses the UCERF3 fault system solutions
(conus-2018.R1/dynamic/disagg/<lon>/<lat>/760/2475, PGA, nshmp-lib): the closest distance of the
listed subsections.
"""
import importlib.util
import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ucla_plha.geometry import geometry

UTILITIES = Path(__file__).parents[1] / "utilities"
MODELS = files("ucla_plha").joinpath("source_models/fault_source_models")
SECTIONS = {
    "ucerf3_fm31": UTILITIES / "FM3_1_branch_averaged/ruptures/fault_sections.geojson",
    "ucerf3_fm32": UTILITIES / "FM3_2_branch_averaged/ruptures/fault_sections.geojson",
}
pytestmark = pytest.mark.skipif(
    not all(p.exists() for p in SECTIONS.values()), reason="UCERF3 input files not available"
)


def _converter():
    spec = importlib.util.spec_from_file_location(
        "convert_ucerf3_fault_source_models", UTILITIES / "convert_ucerf3_fault_source_models.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sections(model):
    return {f["properties"]["FaultID"]: f for f in json.loads(SECTIONS[model].read_text())["features"]}


def _section_id(model, name):
    return next(i for i, f in _sections(model).items() if f["properties"]["FaultName"] == name)


def _load(model, name):
    return np.load(str(MODELS.joinpath(model, name)))


def _lat_lon_depth(xyz):
    r = np.linalg.norm(xyz, axis=-1)
    lat = np.degrees(np.arcsin(xyz[..., 2] / r))
    lon = np.degrees(np.arctan2(xyz[..., 1], xyz[..., 0]))
    a, b = 6378.1370, 6356.7523
    c, s = np.cos(np.radians(lat)), np.sin(np.radians(lat))
    radius = np.sqrt(((a**2 * c) ** 2 + (b**2 * s) ** 2) / ((a * c) ** 2 + (b * s) ** 2))
    return lat, lon, radius - r


def _horizontal(lat1, lon1, lat2, lon2):
    """Horizontal distance (km, spherical) and azimuth (degrees) from point 1 to point 2."""
    p1, p2, dl = np.radians(lat1), np.radians(lat2), np.radians(lon2 - lon1)
    h = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    az = np.degrees(np.arctan2(np.sin(dl) * np.cos(p2), np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)))
    return 2 * 6371.0 * np.arcsin(np.sqrt(h)), az % 360.0


@pytest.mark.parametrize("model", list(SECTIONS))
def test_reduced_section_areas_are_the_ucerf3_rupture_areas(model):
    # The UCERF3 rupture areas (and magnitudes) use the aseismicity-reduced subsection areas
    conv = _converter()
    _, _, _, area_reduced = conv.get_section_properties(str(SECTIONS[model]), True)
    _, _, _, area_full = conv.get_section_properties(str(SECTIONS[model]), False)
    props = pd.read_csv(SECTIONS[model].parent / "properties.csv")
    rs = _load(model, "ruptures_segments.npz")
    seg = rs["segment_index"]
    b = np.r_[0, np.flatnonzero(np.diff(rs["rupture_index"])) + 1]
    area_ucerf3 = props["Area (m^2)"].values / 1e6
    np.testing.assert_allclose(np.add.reduceat(area_reduced[seg], b), area_ucerf3, rtol=2e-3)
    assert np.median(np.add.reduceat(area_full[seg], b) / area_ucerf3) > 1.1


def test_airport_lake_reduced_top_edge():
    # Airport Lake, Subsection 0: dip 50 degrees, upper depth 0, lower depth 13 km, aseismic slip
    # factor 0.1, so the upper edge is at 1.3 km and 1.3 / tan(50) = 1.091 km down-dip of the trace
    # and the lower edge at 13 km and 13 / tan(50) = 10.908 km from the trace.
    model = "ucerf3_fm31"
    sec = _sections(model)[0]
    p = sec["properties"]
    assert (p["DipDeg"], p["UpDepth"], p["LowDepth"]) == (50.0, 0.0, 13.0)
    assert p["AseismicSlipFactor"] == pytest.approx(0.1)
    seg = _load(model, "rect_segment_id.npy")
    i = int(np.flatnonzero(seg == 0)[0])
    trace = np.asarray(sec["geometry"]["coordinates"])
    lat, lon, _ = _lat_lon_depth(_load(model, "rect_rjb.npy")[i])
    tlat, tlon, tdepth = _lat_lon_depth(_load(model, "tri_rrup.npy")[i])  # corners 1, 2, 4
    np.testing.assert_allclose(tdepth, [1.3, 1.3, 13.0], atol=1e-6)
    for corner, point, expected in ((0, 0, 1.3), (1, 1, 1.3), (2, 0, 13.0), (3, 1, 13.0)):
        dist, az = _horizontal(trace[point, 1], trace[point, 0], lat[corner], lon[corner])
        assert dist == pytest.approx(expected / np.tan(np.radians(50.0)), abs=0.005)
        assert az == pytest.approx(p["DipDir"], abs=0.5)
    r = _load(model, "ruptures.npz")
    # rupture 0 is Airport Lake subsections 0 and 1 (both with aseismic slip factor 0.1)
    assert r["ztor"][0] == pytest.approx(1.3)
    assert r["zbor"][0] == pytest.approx(13.0)


def test_multi_section_rupture_depths_are_area_weighted():
    # Rupture 2025 of FM3.1: Garlock (Central) subsections 18 and 19 (reduced upper depth 0, lower
    # depth 11.5 km) and Panamint Valley subsections 0 and 1 (reduced upper depth 1.3, lower depth
    # 13 km), all vertical. ztor and zbor are the averages weighted by the reduced subsection areas,
    # as OpenSHA (CompoundSurface) and the nshm-conus 5.3.1 UCERF3 ruptures.csv (depth 0.704 km,
    # depth + width = 12.312 km), not the shallowest and deepest depths (0 and 13 km).
    model = "ucerf3_fm31"
    sections = _sections(model)
    rs = _load(model, "ruptures_segments.npz")
    seg = rs["segment_index"][rs["rupture_index"] == 2025]
    assert sorted(seg) == [612, 613, 1523, 1524]
    upper, lower, area = [], [], []
    for s in seg:
        p = sections[s]["properties"]
        u = p["UpDepth"] + p["AseismicSlipFactor"] * (p["LowDepth"] - p["UpDepth"])
        c = np.asarray(sections[s]["geometry"]["coordinates"])
        length = sum(_horizontal(c[i, 1], c[i, 0], c[i + 1, 1], c[i + 1, 0])[0] for i in range(len(c) - 1))
        upper.append(u)
        lower.append(p["LowDepth"])
        area.append(length * (p["LowDepth"] - u) / np.sin(np.radians(p["DipDeg"])))
    r = _load(model, "ruptures.npz")
    assert r["ztor"][2025] == pytest.approx(np.average(upper, weights=area), abs=1e-3)
    assert r["zbor"][2025] == pytest.approx(np.average(lower, weights=area), abs=1e-3)
    assert r["ztor"][2025] == pytest.approx(0.704, abs=0.002)
    assert r["zbor"][2025] == pytest.approx(12.312, abs=0.002)
    # the earlier convention is still available
    conv = _converter()
    up, low, _, _ = conv.get_section_properties(str(SECTIONS[model]), True)
    assert min(up[seg]) == pytest.approx(0.0) and max(low[seg]) == pytest.approx(13.0)


def test_without_aseismicity_upper_edge_is_the_trace():
    conv = _converter()
    out = conv.get_lat_lon(str(SECTIONS["ucerf3_fm31"]), apply_aseismicity=False)
    segment_id, lat1, lon1, d1 = out[:4]
    sections = _sections("ucerf3_fm31")
    first = np.r_[True, np.diff(segment_id) != 0]
    trace = np.array([sections[s]["geometry"]["coordinates"][0][:2] for s in segment_id[first]])
    np.testing.assert_array_equal(lon1[first], trace[:, 0])
    np.testing.assert_array_equal(lat1[first], trace[:, 1])
    np.testing.assert_array_equal(d1, [sections[s]["properties"]["UpDepth"] for s in segment_id])


def test_vertical_sections_do_not_move():
    conv = _converter()
    filename = str(SECTIONS["ucerf3_fm31"])
    new = conv.get_lat_lon(filename, apply_aseismicity=True)
    old = conv.get_lat_lon(filename, apply_aseismicity=False)
    vertical = new[13] == 90.0
    assert vertical.sum() > 1000
    for k in (1, 2, 4, 5, 7, 8, 10, 11):  # latitudes and longitudes of the corners
        np.testing.assert_array_equal(new[k][vertical], old[k][vertical])


# site (lon, lat): [(model, subsection, USGS 2018 rRup km)]
USGS_2018_RRUP = {
    (-119.29, 34.28): [
        ("ucerf3_fm31", "Ventura-Pitas Point, Subsection 4", 2.645),
        ("ucerf3_fm31", "Ventura-Pitas Point, Subsection 3", 2.917),
        ("ucerf3_fm31", "Red Mountain, Subsection 0", 7.703),
        ("ucerf3_fm31", "Oak Ridge (Onshore), Subsection 0", 8.514),
        ("ucerf3_fm32", "Oak Ridge (Offshore), Subsection 5", 4.766),
    ],
    (-118.25, 34.05): [
        ("ucerf3_fm31", "Elysian Park (Upper), Subsection 1", 5.400),
        ("ucerf3_fm31", "Puente Hills, Subsection 4", 5.825),
        ("ucerf3_fm31", "Hollywood, Subsection 0", 8.341),
        ("ucerf3_fm32", "San Vicente, Subsection 0", 6.307),
    ],
}


@pytest.mark.parametrize("site", list(USGS_2018_RRUP))
def test_section_rrup_matches_usgs_2018(site):
    p_xyz = geometry.point_to_xyz(np.array([site[1], site[0], 0.0]))
    for model, name, expected in USGS_2018_RRUP[site]:
        rrup = geometry.point_triangle_distance(
            _load(model, "tri_rrup.npy"), p_xyz, _load(model, "tri_segment_id.npy")
        )
        # nshmp-lib uses a 1 km grid; the full (unreduced) planes were 0.4 to 1.6 km closer
        assert rrup[_section_id(model, name)] == pytest.approx(expected, abs=0.05), name

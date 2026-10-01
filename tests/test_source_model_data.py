"""
Integrity tests for the source model files shipped in src/ucla_plha/source_models. These check
that every fault and point source model is internally consistent and physically reasonable,
and compare against the source model input files in utilities/ when they are available.
"""

import json
from importlib.resources import files
from pathlib import Path

import numpy as np
import pytest

from ucla_plha import plha
from ucla_plha.geometry import geometry

SOURCE_MODELS = files("ucla_plha").joinpath("source_models")
FAULT_MODELS = sorted(
    p.name
    for p in Path(str(SOURCE_MODELS), "fault_source_models").iterdir()
    if p.is_dir()
)
POINT_MODELS = sorted(
    p.name
    for p in Path(str(SOURCE_MODELS), "point_source_models").iterdir()
    if p.is_dir()
)
UTILITIES = Path(__file__).parents[1] / "utilities"
SECTION_FILES = {
    "ucerf3_fm31": (
        UTILITIES / "FM3_1_branch_averaged/ruptures/fault_sections.geojson",
        "FaultID",
        "DipDeg",
        "DipDir",
    ),
    "ucerf3_fm32": (
        UTILITIES / "FM3_2_branch_averaged/ruptures/fault_sections.geojson",
        "FaultID",
        "DipDeg",
        "DipDir",
    ),
    "nshm23_wus": (
        UTILITIES / "nshm23_wus_branch_avg/sections.geojson",
        "index",
        "dip",
        "dip-direction",
    ),
}
GRID_INPUT_FILES = {
    "ucerf3_fm31_grid_sub_seis": (
        "FM3_1_branch_averaged/solution/grid_sub_seis_mfds.csv",
        "FM3_1_branch_averaged/solution/grid_mech_weights.csv",
    ),
    "ucerf3_fm31_grid_unassociated": (
        "FM3_1_branch_averaged/solution/grid_unassociated_mfds.csv",
        "FM3_1_branch_averaged/solution/grid_mech_weights.csv",
    ),
    "ucerf3_fm32_grid_sub_seis": (
        "FM3_2_branch_averaged/solution/grid_sub_seis_mfds.csv",
        "FM3_2_branch_averaged/solution/grid_mech_weights.csv",
    ),
    "ucerf3_fm32_grid_unassociated": (
        "FM3_2_branch_averaged/solution/grid_unassociated_mfds.csv",
        "FM3_2_branch_averaged/solution/grid_mech_weights.csv",
    ),
}
ALL_GMMS = ["bssa14", "ask14", "cb14", "cy14"]


def _load(source_type, model, name):
    return np.load(str(SOURCE_MODELS.joinpath(source_type, model, name)))


def _xyz_to_latlon(xyz):
    r = np.linalg.norm(xyz, axis=-1)
    return np.degrees(np.arcsin(xyz[..., 2] / r)), np.degrees(
        np.arctan2(xyz[..., 1], xyz[..., 0])
    )


@pytest.fixture(scope="module", params=FAULT_MODELS)
def fault_model(request):
    model = request.param
    data = {
        "model": model,
        "ruptures": dict(_load("fault_source_models", model, "ruptures.npz")),
        "ruptures_segments": dict(
            _load("fault_source_models", model, "ruptures_segments.npz")
        ),
    }
    for name in [
        "tri_segment_id",
        "tri_rrup",
        "tri_rjb",
        "rect_segment_id",
        "rect_rjb",
    ]:
        data[name] = _load("fault_source_models", model, name + ".npy")
    return data


def test_fault_geometry_arrays_are_consistent(fault_model):
    d = fault_model
    n_rect = len(d["rect_segment_id"])
    assert d["rect_rjb"].shape == (n_rect, 4, 3)
    assert d["tri_rrup"].shape == d["tri_rjb"].shape == (2 * n_rect, 3, 3)
    np.testing.assert_array_equal(
        d["tri_segment_id"], np.concatenate((d["rect_segment_id"],) * 2)
    )
    # Distance arrays are indexed by segment_id, so ids must be grouped and run 0, 1, ..., N-1
    assert np.all(np.diff(d["rect_segment_id"]) >= 0)
    np.testing.assert_array_equal(
        np.unique(d["rect_segment_id"]), np.arange(d["rect_segment_id"].max() + 1)
    )
    for name in ["tri_rrup", "tri_rjb", "rect_rjb"]:
        assert np.all(np.isfinite(d[name])), name


def test_fault_surface_projections_lie_on_surface(fault_model):
    d = fault_model
    # Rjb triangles and Rx rectangles are surface projections, so all points have zero depth.
    # Rrup triangles have the same latitude and longitude as Rjb triangles but are at depth.
    lat, lon = _xyz_to_latlon(d["tri_rjb"])
    surface = geometry.point_to_xyz
    sample = np.unique(np.linspace(0, len(lat) - 1, 200).astype(int))
    for i in sample:
        for j in range(3):
            np.testing.assert_allclose(
                d["tri_rjb"][i, j],
                surface(np.array([lat[i, j], lon[i, j], 0.0])),
                atol=1e-6,
            )
    lat_rrup, lon_rrup = _xyz_to_latlon(d["tri_rrup"])
    np.testing.assert_allclose(lat_rrup, lat, atol=1e-9)
    np.testing.assert_allclose(lon_rrup, lon, atol=1e-9)


def test_fault_rupture_properties_are_physical(fault_model):
    r = fault_model["ruptures"]
    n = len(r["m"])
    for name in ["rate", "fault_type", "dip", "ztor", "zbor"]:
        assert len(r[name]) == n, name
    # UCERF3 assigns a zero rate to a small number of ruptures
    assert np.all(r["rate"] >= 0.0)
    assert np.mean(r["rate"] > 0.0) > 0.99
    assert np.all((r["m"] > 4.0) & (r["m"] < 9.5))
    assert set(np.unique(r["fault_type"])) <= {1, 2, 3}
    assert np.all((r["dip"] > 0.0) & (r["dip"] <= 90.0))
    assert np.all(r["ztor"] >= 0.0)
    assert np.all(r["zbor"] > r["ztor"])
    assert np.all(r["zbor"] < 40.0)
    # Depths and dips should vary across ruptures rather than taking a handful of values
    assert len(np.unique(r["ztor"])) > 5
    assert len(np.unique(r["zbor"])) > 5
    assert len(np.unique(r["dip"])) > 5


def test_fault_rupture_segment_mapping(fault_model):
    rs = fault_model["ruptures_segments"]
    n_rup = len(fault_model["ruptures"]["m"])
    n_seg = fault_model["rect_segment_id"].max() + 1
    assert len(rs["rupture_index"]) == len(rs["segment_index"])
    # get_source_data uses np.minimum.reduceat, which requires ruptures to be grouped in order
    assert np.all(np.diff(rs["rupture_index"]) >= 0)
    np.testing.assert_array_equal(np.unique(rs["rupture_index"]), np.arange(n_rup))
    assert rs["segment_index"].min() >= 0
    assert rs["segment_index"].max() < n_seg


def test_fault_down_dip_edge_follows_dip_direction(fault_model):
    filename, id_key, dip_key, dipdir_key = SECTION_FILES[fault_model["model"]]
    if not filename.exists():
        pytest.skip(f"{filename} not available")
    features = json.loads(filename.read_text())["features"]
    props = {f["properties"][id_key]: f["properties"] for f in features}
    seg = fault_model["rect_segment_id"]
    dip = np.array([props[s][dip_key] for s in seg])
    dipdir = np.array([props[s][dipdir_key] for s in seg])
    lat, lon = _xyz_to_latlon(fault_model["rect_rjb"])
    # Azimuth from top corner p1 to the down-dip corner p3
    dlat = np.radians(lat[:, 2] - lat[:, 0])
    dlon = np.radians(lon[:, 2] - lon[:, 0]) * np.cos(np.radians(lat[:, 0]))
    azimuth = np.degrees(np.arctan2(dlon, dlat)) % 360.0
    check = (dip < 85.0) & (np.hypot(dlat, dlon) * 6371.0 > 0.5)
    assert check.sum() > 100
    error = (azimuth[check] - dipdir[check] + 180.0) % 360.0 - 180.0
    assert np.max(np.abs(error)) < 2.0


@pytest.mark.parametrize("model", FAULT_MODELS)
@pytest.mark.parametrize(
    "site", [(36.80547, -121.786074), (34.05, -118.25), (40.76, -111.89)]
)
def test_fault_source_data_is_finite(model, site):
    p_xyz = geometry.point_to_xyz(np.array([site[0], site[1], 0.0]))
    out = plha.get_source_data(
        "fault_source_models", model, p_xyz, 200.0, 5.0, ALL_GMMS
    )
    names = [
        "m",
        "fault_type",
        "rate",
        "rjb",
        "rrup",
        "rx",
        "rx1",
        "ry0",
        "dip",
        "ztor",
        "zbor",
    ]
    for name, a in zip(names, out):
        assert np.all(np.isfinite(a)), f"{model} {name} has non-finite values at {site}"
    rjb, rrup = out[3], out[4]
    # Rrup is at least Rjb, apart from a small effect of Earth's curvature (about 30 m at 200 km)
    assert np.all(rrup >= rjb * (1.0 - 1e-3))


@pytest.mark.parametrize("model", POINT_MODELS)
def test_point_source_model_is_consistent(model):
    points = _load("point_source_models", model, "points.npy")
    node_index = _load("point_source_models", model, "node_index.npy")
    r = dict(_load("point_source_models", model, "ruptures.npz"))
    assert points.shape == (len(node_index), 3)
    assert len(np.unique(node_index)) == len(node_index)
    lat, lon = _xyz_to_latlon(points)
    assert np.all((lat > 20.0) & (lat < 55.0) & (lon > -130.0) & (lon < -100.0))
    n = len(r["m"])
    for name in ["rate", "style", "node_index"]:
        assert len(r[name]) == n, name
    assert np.all(r["rate"] > 0.0)
    assert np.all((r["m"] > 4.0) & (r["m"] < 9.5))
    assert set(np.unique(r["style"])) <= {1, 2, 3}
    assert np.all(np.isin(r["node_index"], node_index))


@pytest.mark.parametrize("model", sorted(GRID_INPUT_FILES))
def test_point_source_style_rates_match_input(model):
    pd = pytest.importorskip("pandas")
    mfd_file, weights_file = (UTILITIES / f for f in GRID_INPUT_FILES[model])
    if not (mfd_file.exists() and weights_file.exists()):
        pytest.skip("UCERF3 grid input files not available")
    mfds = pd.read_csv(mfd_file)
    weights = pd.read_csv(weights_file).set_index("Node Index").loc[mfds["Node Index"]]
    m_columns = [
        c
        for c in mfds.columns
        if c not in ("Node Index", "Latitude", "Longitude") and float(c) > 4.0
    ]
    total = mfds[m_columns].sum(axis=1).values
    r = dict(_load("point_source_models", model, "ruptures.npz"))
    # fault_type: 1 = reverse, 2 = normal, 3 = strike slip
    for style, column in [
        (1, "Fraction Reverse"),
        (2, "Fraction Normal"),
        (3, "Fraction Strike-Slip"),
    ]:
        expected = np.sum(weights[column].values * total)
        assert r["rate"][r["style"] == style].sum() == pytest.approx(
            expected, rel=1e-6
        ), column

"""
Known-value tests for the distance calculations in ucla_plha.geometry. Expected values are
computed analytically or by brute force, independently of the code under test.
"""

import numpy as np
import pytest

from ucla_plha.geometry import geometry


def _brute_force_triangle_distance(tri, p, n=400):
    """Minimum distance from p to a dense barycentric sampling of triangle tri"""
    u, v = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    keep = u + v <= 1.0
    u = u[keep]
    v = v[keep]
    pts = tri[0] + u[:, None] * (tri[1] - tri[0]) + v[:, None] * (tri[2] - tri[0])
    return np.min(np.linalg.norm(pts - p, axis=1))


@pytest.mark.parametrize(
    "p,expected",
    [
        ((0.2, 0.2, 3.0), 3.0),  # above the interior
        ((0.5, -2.0, 0.0), 2.0),  # beside edge p0-p1
        ((-3.0, -4.0, 0.0), 5.0),  # beyond vertex p0
        ((2.0, 2.0, 0.0), np.sqrt(4.5)),  # beyond hypotenuse
        ((0.25, 0.25, 0.0), 0.0),  # in the plane of the triangle
    ],
)
def test_point_triangle_distance_analytic(p, expected):
    tri = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]])
    dist = geometry.point_triangle_distance(tri, np.array(p), np.array([0]))
    np.testing.assert_allclose(dist, [expected], atol=1e-9)


def test_point_triangle_distance_matches_brute_force():
    rng = np.random.default_rng(42)
    tri = rng.uniform(-10.0, 10.0, size=(50, 3, 3))
    p = rng.uniform(-15.0, 15.0, size=3)
    dist = geometry.point_triangle_distance(tri, p, np.arange(50))
    expected = np.array([_brute_force_triangle_distance(t, p) for t in tri])
    # Sampling is coarser than the exact solution, so it can only overestimate the distance
    assert np.all(dist <= expected + 1e-9)
    np.testing.assert_allclose(dist, expected, atol=0.05)


def test_point_triangle_distance_takes_minimum_over_segment():
    tri = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[0.0, 0.0, 5.0], [1.0, 0.0, 5.0], [0.0, 1.0, 5.0]],
            [[0.0, 0.0, 9.0], [1.0, 0.0, 9.0], [0.0, 1.0, 9.0]],
        ]
    )
    # Triangles 0 and 2 belong to segment 0, and triangle 1 belongs to segment 1
    dist = geometry.point_triangle_distance(
        tri, np.array([0.2, 0.2, 8.0]), np.array([0, 1, 0])
    )
    np.testing.assert_allclose(dist, [1.0, 3.0], atol=1e-9)


# Surface projection of a rupture whose top edge runs from p1 = (0, 0) to p2 = (10, 0) and dips
# toward +y with a horizontal width of 5. Points p3 and p4 are the projections of the bottom
# edge directly down-dip of p1 and p2.
RECT = np.array([[[0, 0, 0], [10, 0, 0], [0, 5, 0], [10, 5, 0]]], dtype=float)


@pytest.mark.parametrize(
    "site,rx,rx1",
    [
        ((5.0, 8.0), 8.0, 3.0),  # hanging wall, beyond the bottom edge
        ((5.0, 2.0), 2.0, -3.0),  # hanging wall, above the rupture
        ((5.0, -4.0), -4.0, -9.0),  # footwall
    ],
    ids=["beyond-bottom-edge", "above-rupture", "footwall"],
)
def test_rx_rx1_sign_convention(site, rx, rx1):
    # Rx is positive on the hanging wall side of the top edge, and Rx1 = Rx - W cos(dip)
    out_rx, out_rx1, out_ry0 = geometry.get_Rx_Rx1_Ry0(
        RECT, np.array([site[0], site[1], 0.0]), np.array([0])
    )
    np.testing.assert_allclose(out_rx, [rx], atol=1e-9)
    np.testing.assert_allclose(out_rx1, [rx1], atol=1e-9)
    np.testing.assert_allclose(out_ry0, [0.0], atol=1e-9)


@pytest.mark.parametrize(
    "site,ry0",
    [((15.0, 2.0), 5.0), ((-3.0, 8.0), 3.0), ((4.0, 2.0), 0.0)],
)
def test_ry0_distance_off_end_of_rupture(site, ry0):
    _, _, out_ry0 = geometry.get_Rx_Rx1_Ry0(
        RECT, np.array([site[0], site[1], 0.0]), np.array([0])
    )
    np.testing.assert_allclose(out_ry0, [ry0], atol=1e-9)


def test_point_to_xyz_elevation_and_poles():
    np.testing.assert_allclose(
        geometry.point_to_xyz(np.array([90.0, 0.0, 0.0])),
        [0.0, 0.0, 6356.7523],
        atol=1e-6,
    )
    p0 = geometry.point_to_xyz(np.array([34.0, -118.0, 0.0]))
    p1 = geometry.point_to_xyz(np.array([34.0, -118.0, 1.5]))
    np.testing.assert_allclose(np.linalg.norm(p1) - np.linalg.norm(p0), 1.5, atol=1e-9)


def test_point_to_xyz_great_circle_distance():
    # One degree of latitude along a meridian is 110.95 km on the WGS84 spheroid at 34.5 degrees.
    # point_to_xyz treats latitude as geocentric, which is within 0.3 km of that.
    p0 = geometry.point_to_xyz(np.array([34.0, -118.0, 0.0]))
    p1 = geometry.point_to_xyz(np.array([35.0, -118.0, 0.0]))
    assert np.linalg.norm(p1 - p0) == pytest.approx(110.95, abs=0.3)

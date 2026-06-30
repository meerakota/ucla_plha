import numpy as np

from ucla_plha.geometry import geometry


def test_point_to_xyz_equator_reference():
    xyz = geometry.point_to_xyz(np.array([0.0, 0.0, 0.0]))
    assert xyz.shape == (3,)
    np.testing.assert_allclose(xyz, [6378.1370, 0.0, 0.0], atol=1e-6)
    xyz2 = geometry.point_to_xyz(np.array([36.8, -121.8, 0.5]))
    assert np.linalg.norm(xyz2) > 6356.0


def test_point_triangle_distance_regions_and_degenerate():
    tri = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [[5.0, 5.0, 5.0], [5.0, 5.0, 5.0], [5.0, 5.0, 5.0]],
            [[2.0, 2.0, 0.0], [3.0, 2.0, 0.0], [2.0, 3.0, 0.0]],
        ]
    )
    p_xyz = np.array([0.0, 0.0, 1.0])
    tri_segment_id = np.array([0, 1, 1])

    dist = geometry.point_triangle_distance(tri, p_xyz, tri_segment_id)

    assert dist.shape == (2,)
    assert np.all(dist >= 0.0)
    assert np.all(np.isfinite(dist))
    np.testing.assert_allclose(dist[0], 1.0, atol=1e-6)


def test_get_rx_rx1_ry0_shapes_and_grouping():
    rect = np.array(
        [
            [[0, 0, 0], [10, 0, 0], [0, 0, -5], [10, 0, -5]],
            [[10, 0, 0], [20, 0, 0], [10, 0, -5], [20, 0, -5]],
            [[0, 30, 0], [10, 30, 0], [0, 30, -5], [10, 30, -5]],
        ],
        dtype=float,
    )
    point = np.array([5.0, 8.0, 0.0])
    rect_segment_id = np.array([0, 0, 1])

    rx, rx1, ry0 = geometry.get_Rx_Rx1_Ry0(rect, point, rect_segment_id)

    assert rx.shape == rx1.shape == ry0.shape == (2,)
    assert np.all(np.isfinite(rx))
    assert np.all(ry0 >= 0.0)

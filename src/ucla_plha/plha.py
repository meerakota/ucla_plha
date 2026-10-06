import json
import os
from importlib.resources import files

import jsonschema
import numpy as np
from scipy.special import ndtr

from ucla_plha.liquefaction_models import (
    boulanger_idriss_2012,
    cetin_et_al_2018,
    ngl_smt_2024,
    moss_et_al_2006,
    boulanger_idriss_2016,
)
from ucla_plha.ground_motion_models import ask14, bssa14, cb14, cy14, idriss14
from ucla_plha.geometry import geometry, point_source
from ucla_plha import nshm_site_data, pygmm_gmms, tectonic_regions


def decompress_ucerf3_source_data():
    """Decompresses ruptures.gz and ruptures_segments.gz in source_models/fault_source_models/ucerf3_fm31
    and source_models/fault_source_models/ucerf3_fm32. These files are compressed to reduce the
    size of the ucla_plha package, but decompressing the files each time the code is run is inefficient.
    This function decompresses them within the installation directory. The get_source_data() function
    checks for .pkl versions of the files, and uses them if they exist. Otherwise it uses the
    .gz versions of the files and decompresses them at runtime. This function needs to be run only
    once.
    """
    branches = ["ucerf3_fm31", "ucerf3_fm32"]
    for branch in branches:
        path = files("ucla_plha").joinpath(
            "source_models/fault_source_models/" + branch
        )
        ruptures = np.load(str(path.joinpath("ruptures.npz")))
        np.savez(str(path.joinpath("ruptures.npy")), **ruptures)
        ruptures_segments = np.load(str(path.joinpath("ruptures_segments.npz")))
        np.savez(str(path.joinpath("ruptures_segments.npy")), **ruptures_segments)


def _source_model_path(source_type, source_model):
    """Directory of a source model in the package (source_models/<source_type>/<source_model>)."""
    return files("ucla_plha").joinpath(
        "source_models/" + source_type + "/" + source_model
    )


def get_source_info(source_type, source_model):
    """Returns the source_info.json contents of a source model, with defaults.

    Args:
        source_type (string): Either "fault_source_models" or "point_source_models"
        source_model (string): Directory for source_model within the source_type directory

    Returns:
        info (dict): see tectonic_regions.read_source_info. Source models without a
            source_info.json file are active crust models.
    """
    path = _source_model_path(source_type, source_model)
    if not os.path.isdir(str(path)):
        raise ValueError(f'source model "{source_model}" not found in {source_type}')
    return tectonic_regions.read_source_info(str(path), source_type, source_model)


def _distance_needs(gmms):
    """Set of distances ("rjb", "rrup", "rx", "ry0") used by the ground motion models."""
    needs = set()
    for gmm in gmms:
        try:
            needs |= pygmm_gmms.resolve(gmm).distances
        except (pygmm_gmms.UnavailableGmmError, ValueError, ImportError):
            continue
    return needs


def _filter_mask(n, distance, m, dist_cutoff, m_min):
    if (dist_cutoff is not None) and (m_min is not None):
        return (distance < dist_cutoff) & (m >= m_min)
    elif (dist_cutoff is not None) and (m_min is None):
        return distance < dist_cutoff
    elif (dist_cutoff is None) and (m_min is not None):
        return m >= m_min
    return np.full(n, True)


def _closest_segment(rrup_seg, rrup, segment_index, boundaries):
    """Segment of each rupture with the smallest rrup (the lowest segment index if tied).

    rrup_seg: rrup of each (rupture, segment) entry; rrup: minimum rrup of each rupture;
    segment_index: segment of each entry; boundaries: first entry of each rupture.
    """
    counts = np.diff(np.r_[boundaries, len(rrup_seg)])
    is_min = rrup_seg == np.repeat(rrup, counts)
    big = np.iinfo(segment_index.dtype).max
    return np.minimum.reduceat(np.where(is_min, segment_index, big), boundaries)


_SECTION_GRIDS = {}


def _section_grids(path):
    """nshmp-lib gridded surfaces of the fault sections in sections.npz (built once per model)."""
    key = str(path)
    if key not in _SECTION_GRIDS:
        with np.load(str(path.joinpath("sections.npz"))) as sections:
            _SECTION_GRIDS[key] = geometry.nshmp_section_grids(sections)
    return _SECTION_GRIDS[key]


def get_source_data(
    source_type,
    source_model,
    p_xyz,
    dist_cutoff,
    m_min,
    gmms,
    extras=False,
    source_info=None,
    rupture_rx="closest_section",
    fault_distances=None,
    section_rx=None,
):
    """Returns magnitude, fault type, rate, distance, and fault geometry terms.

    Args:
        source_type (string): Either "fault_source_models" or "point_source_models"
        source_model (string): Directory for source_model within the source_type directory.
            For example "ucerf3_fm31", "ucerf3_fm32", or "nshm23_wus" for fault_source_models, and
            "ucerf3_fm31_grid_sub_seis", "ucerf3_fm31_grid_unassociated", "ucerf3_fm32_grid_sub_seis",
            "ucerf3_fm32_grid_unassociated", or "nshm23_wus_grid" for point_source_models
        p_xyz (numpy array, dtype=float): Array containing x, y, z coordinates for point of interest, length = 3
        dist_cutoff (float): maximum distance to consider in seismic hazard analysis
        m_min (float): minimum magnitude to consider in seismic hazard analysis
        gmms (array, dtype=string): An array of strings defining ground motion models to use in seismic
            hazard analysis: "ask14", "bssa14", "cb14", "cy14", "idriss14", nshmp-lib Gmm ids, or
            pygmm class names. They determine which distances are computed and whether the distance
            cutoff is applied to rjb (if any model uses rjb) or rrup.
        extras (bool): if True, also return a dict of additional rupture data (see Returns)
        source_info (dict): source_info of the source model (read from source_info.json if None)
        rupture_rx (str): how rx and rx1 of a rupture with several fault segments are found.
            "closest_section" (default): rx and rx1 of the segment with the smallest rrup (ties go
            to the lowest segment index), as nshmp-lib does for fault system ruptures
            (SystemRuptureSet.InputGenerator). "minimum": the minimum rx and rx1 over the
            segments of the rupture (ucla_plha 1.x and earlier); because rx is signed, this is
            the most footwall-side segment of the rupture, which may be far from the site.
            rjb, rrup, and ry0 are always the minima over the segments.
        fault_distances (str): how the distances of fault source ruptures are computed. None
            (default): the "fault_distances" of source_info ("triangles" if not given).
            "triangles": distances to the triangles (tri_rjb.npy, tri_rrup.npy) and rectangles
            (rect_rjb.npy) of the fault segments. "nshmp_grid": rjb, rrup, and rx as nshmp-lib
            computes them (model.Distance.compute) from its gridded surfaces: per rupture from
            grid.npz (grid points and the grid subset of each rupture;
            geometry.gridded_rupture_distances; NSHM23 Cascadia interface and CEUS fault
            models), with rx1 = rx - W cos(dip), or per segment from sections.npz (traces and
            properties of the fault sections, from which the 1 km nshmp-lib grids are built;
            geometry.gridded_section_distances; NSHM23 WUS), combined over the segments of a
            rupture as for "triangles" (rx1 from the horizontal width of the segment). ry0 is
            always from the rectangles.
        section_rx (str): how rx and rx1 of a fault segment with several planar pieces
            (rectangles) are found, for "triangles" distances. None (default): the
            "section_rx" of source_info ("minimum" if not given). "minimum": the minimum over
            the pieces (ucla_plha 2.x and earlier). "extended_trace": as nshmp-lib does
            (model.Distance.getDistanceX), the distance to the upper edge of the whole segment
            extended 1000 km along strike at both ends, positive on the dip side
            (geometry.section_rx_extended_trace). The NSHM23 WUS model uses "extended_trace".

    Returns: A tuple containing the following arrays
        m (array, dtype=float): Numpy array of magnitudes, length = N
        fault_type (array, dtype=int): Numpy array of fault types. 1 = reverse, 2 = normal, 3 = strike slip, length = N
        rate (array, dtype=float): Numpy array of the rate of occurrence of each event, length = N
        rjb (array, dtype=float): Numpy array of Joyner-Boore distances between the site and each rupture, length = N
        rrup (array, dtype=float): Numpy array of rupture distance between the site and each rupture, length = N
        rx (array, dtype=float): Numpy array of distance from site to the surface projection of the top of each
            rupture, measured perpendicular to the strike, length = N
        rx1 (array, dtype=float): Numpy array of distance from site to the surface projection of the bottom of each
            rupture, measured perpendicular to the strike, length = N
        ry0 (array, dtype=float): Numpy array of the distance from the site to the surface projection of each
            rupture, measured parallel to the strike, length = N
        dip (array, dtype=float): Numpy array of dip angle for each rupture in degrees, length = N
        ztor (array, dtype=float): Numpy array of depth to the top of each rupture in km, length = N
        zbor (array, dtype=float): Numpy array of depth to the bottom of each rupture in km, length = N

        If extras is True, a tuple (arrays, extras) is returned, where extras is a dict with
        "index" (index of each returned rupture in ruptures.npz) and, for cluster models,
        "cluster_id" (-1 for independent ruptures), "cluster_fault_id", "cluster_weight"
        (weight of the rupture within its fault of the cluster), and "cluster_rate" (rate times
        logic tree weight of the cluster of the rupture), each of length N.

    Notes:
        N = number of events
        Point source models use the distance treatment given by "point_source" in their
        source_info.json: the ucla_plha crustal approximations (default, used by the UCERF3 and
        NSHM23 WUS grids) or the nshmp-lib point source conventions (geometry.point_source).
    """
    if source_type not in ("fault_source_models", "point_source_models"):
        return None
    path = _source_model_path(source_type, source_model)
    if source_info is None:
        source_info = get_source_info(source_type, source_model)
    needs = _distance_needs(gmms)
    need_rjb = "rjb" in needs
    need_rrup = bool(needs & {"rrup", "rx", "ry0"})

    if source_type == "fault_source_models":
        # Read files required by all ground motion models
        # Read decompressed version of files if they exist. Otherwise read zipped version.
        if os.path.exists(str(path.joinpath("ruptures.npy.npz"))):
            ruptures = np.load(str(path.joinpath("ruptures.npy.npz")))
        else:
            ruptures = np.load(str(path.joinpath("ruptures.npz")))
        m = ruptures["m"]
        fault_type = ruptures["fault_type"]
        if os.path.exists(str(path.joinpath("ruptures_segments.npy.npz"))):
            ruptures_segments = np.load(str(path.joinpath("ruptures_segments.npy.npz")))
        else:
            ruptures_segments = np.load(str(path.joinpath("ruptures_segments.npz")))
        segment_index = ruptures_segments["segment_index"]
        ruptures_index = ruptures_segments["rupture_index"]
        rate = ruptures["rate"]
        dip = ruptures["dip"]
        ztor = ruptures["ztor"]
        zbor = ruptures["zbor"]

        # Now read files required by the ground motion models
        # bssa14: rjb
        # ask14: rrup,rx,rx1,ry0,
        # cb14: rjb,rrup,rx
        # cy14: rjb,rrup,rx
        # idriss14: rrup
        # pygmm models: the distances in their parameters
        empty_array = np.empty(len(m))
        if fault_distances is None:
            fault_distances = source_info.get("fault_distances", "triangles")
        if fault_distances not in ("triangles", "nshmp_grid"):
            raise ValueError(
                f'fault_distances must be "triangles" or "nshmp_grid", not "{fault_distances}"'
            )
        if section_rx is None:
            section_rx = source_info.get("section_rx", "minimum")
        if section_rx not in ("minimum", "extended_trace"):
            raise ValueError(
                f'section_rx must be "minimum" or "extended_trace", not "{section_rx}"'
            )
        # nshmp_grid: rupture grids (grid.npz, distances per rupture) or section grids
        # (sections.npz, distances per segment, combined over the segments of a rupture)
        rupture_grid = section_grid = False
        if fault_distances == "nshmp_grid":
            rupture_grid = os.path.exists(str(path.joinpath("grid.npz")))
            section_grid = not rupture_grid and os.path.exists(
                str(path.joinpath("sections.npz"))
            )
            if not (rupture_grid or section_grid):
                raise ValueError(
                    f'source model "{source_model}" has no grid.npz or sections.npz for '
                    '"nshmp_grid" distances'
                )
        if rupture_grid or section_grid:
            site_lat = np.degrees(np.arctan2(p_xyz[2], np.hypot(p_xyz[0], p_xyz[1])))
            site_lon = np.degrees(np.arctan2(p_xyz[1], p_xyz[0]))
        if (need_rjb or need_rrup) and not (rupture_grid or section_grid):
            tri_segment_id = np.load(str(path.joinpath("tri_segment_id.npy")))
        if need_rjb and not (rupture_grid or section_grid):
            tri_rjb = np.load(str(path.joinpath("tri_rjb.npy")))
            rjb_all = geometry.point_triangle_distance(tri_rjb, p_xyz, tri_segment_id)
        if need_rrup:
            rect_segment_id = np.load(str(path.joinpath("rect_segment_id.npy")))
            rect = np.load(str(path.joinpath("rect_rjb.npy")))
            rx_all, rx1_all, ry0_all = geometry.get_Rx_Rx1_Ry0(
                rect, p_xyz, rect_segment_id
            )
            if section_rx == "extended_trace" and not rupture_grid:
                rx_all, rx1_all = geometry.section_rx_extended_trace(
                    rect, p_xyz, rect_segment_id
                )
            if not (rupture_grid or section_grid):
                tri_rrup = np.load(str(path.joinpath("tri_rrup.npy")))
                rrup_all = geometry.point_triangle_distance(
                    tri_rrup, p_xyz, tri_segment_id
                )
        if section_grid and (need_rjb or need_rrup):
            grids = _section_grids(path)
            rjb_all, rrup_all, rx_grid = geometry.gridded_section_distances(
                grids, site_lat, site_lon, need_rx=need_rrup
            )
            if need_rrup:
                # rx of the top row of the gridded section; rx1 with the horizontal width of
                # the section (from its rectangles)
                rx1_all = rx_grid - (rx_all - rx1_all)
                rx_all = rx_grid
        split_indices = np.where(np.diff(ruptures_index) != 0)[0] + 1
        boundaries = np.r_[0, split_indices]
        if rupture_grid:
            # nshmp-lib distances to the gridded rupture surfaces
            with np.load(str(path.joinpath("grid.npz"))) as grid:
                rjb_g, rrup_g, rx_g = geometry.gridded_rupture_distances(
                    grid, site_lat, site_lon, need_rx=need_rrup
                )
        if need_rrup:
            if rupture_grid:
                rrup = rrup_g
                rx = rx_g
                # rx1 = rx - W cos(dip), the horizontal width of the rupture
                rx1 = rx - (zbor - ztor) / np.tan(np.radians(np.minimum(dip, 89.999)))
            else:
                rrup_seg = rrup_all[segment_index]
                rrup = np.minimum.reduceat(rrup_seg, boundaries)
                if rupture_rx == "closest_section":
                    closest = _closest_segment(rrup_seg, rrup, segment_index, boundaries)
                    rx = rx_all[closest]
                    rx1 = rx1_all[closest]
                elif rupture_rx == "minimum":
                    rx = np.minimum.reduceat(rx_all[segment_index], boundaries)
                    rx1 = np.minimum.reduceat(rx1_all[segment_index], boundaries)
                else:
                    raise ValueError(
                        f'rupture_rx must be "closest_section" or "minimum", not "{rupture_rx}"'
                    )
            ry0 = np.minimum.reduceat(ry0_all[segment_index], boundaries)
        else:
            rrup = empty_array
            rx = empty_array
            rx1 = empty_array
            ry0 = empty_array

        if need_rjb:
            if rupture_grid:
                rjb = rjb_g
            else:
                rjb = np.minimum.reduceat(rjb_all[segment_index], boundaries)
        else:
            rjb = empty_array

        filter = _filter_mask(len(m), rjb if need_rjb else rrup, m, dist_cutoff, m_min)
        if source_info.get("cluster", False) and dist_cutoff is not None and "cluster_id" in ruptures.files:
            # nshmp-lib (ClusterRuptureSet.checkDistance) uses all ruptures of a cluster if
            # any of its faults is within the maximum distance
            cluster_id = np.asarray(ruptures["cluster_id"])
            near = (rjb if need_rjb else rrup) < dist_cutoff
            near_clusters = np.unique(cluster_id[near & (cluster_id >= 0)])
            whole = (cluster_id >= 0) & np.isin(cluster_id, near_clusters)
            filter = filter | (whole & _filter_mask(len(m), rrup, m, None, m_min))
        arrays = (
            m[filter],
            fault_type[filter],
            rate[filter],
            rjb[filter],
            rrup[filter],
            rx[filter],
            rx1[filter],
            ry0[filter],
            dip[filter],
            ztor[filter],
            zbor[filter],
        )
        if not extras:
            return arrays
        index = np.flatnonzero(filter)
        out_extras = {"index": index}
        if source_info.get("cluster", False):
            out_extras.update(_cluster_extras(path, ruptures, index, source_info))
        return arrays, out_extras

    elif source_type == "point_source_models":
        point_source_info = source_info.get("point_source") or {
            "distance_method": "ucla_plha_crustal"
        }
        if point_source_info["distance_method"] == "nshmp":
            arrays, index = _nshmp_point_source_data(
                path, p_xyz, dist_cutoff, m_min, point_source_info
            )
            return (arrays, {"index": index}) if extras else arrays

        ruptures = np.load(str(path.joinpath("ruptures.npz")))
        rate = ruptures["rate"]
        m = ruptures["m"]
        fault_type = ruptures["style"]
        node_index = ruptures["node_index"]
        nodes_index = np.load(str(path.joinpath("node_index.npy")))
        points = np.load(str(path.joinpath("points.npy")))

        repi_all = np.sqrt(
            (points[:, 0] - p_xyz[0]) ** 2
            + (points[:, 1] - p_xyz[1]) ** 2
            + (points[:, 2] - p_xyz[2]) ** 2
        )
        repi = repi_all[
            np.searchsorted(nodes_index, node_index, sorter=np.argsort(nodes_index))
        ]
        rjb = repi / (1.0 + np.exp(-1.05 * (np.log(repi) - 1.037 * m + 4.2776)))

        filter = _filter_mask(len(rjb), rjb, m, dist_cutoff, m_min)
        dip = np.empty(len(m), dtype=float)

        # using Kaklamanos et al. 2011 guidance for unknown dip, ztor, and zbor
        # Note fault_type = 1 reverse, 2 normal, 3 strike slip
        dip[fault_type == 1] = 40.0
        dip[fault_type == 2] = 50.0
        dip[fault_type == 3] = 90.0
        w = 10.0 ** (-0.76 + 0.27 * m)
        w[fault_type == 1] = 10.0 ** (-1.61 + 0.41 * m[fault_type == 1])
        w[fault_type == 2] = 10.0 ** (-1.14 + 0.35 * m[fault_type == 2])
        zhyp = 5.63 + 0.68 * m
        zhyp[fault_type == 1] = 11.24 - 0.2 * m[fault_type == 1]
        zhyp[fault_type == 2] = 7.08 + 0.61 * m[fault_type == 2]
        ztor = zhyp - 0.6 * w * np.sin(dip * np.pi / 180.0)
        ztor[ztor < 0.0] = 0.0
        # Expression for Rrup is fit to default parameters (first row in table 1) in Thompson and Worden (2017)
        d = np.radians(80)
        zbor = ztor + w * np.sin(dip * np.pi / 180.0)
        rrup1 = np.sqrt((repi - 0.6 * w * np.cos(d)) ** 2 + ztor**2)
        filt = repi < 0.6 * w * np.cos(d) - ztor * np.tan(d)
        rrup1[filt] = (zhyp[filt] - repi[filt] * np.tan(d)) / (
            np.cos(d) + np.sin(d) * np.tan(d)
        )
        rrup2 = zhyp
        filt = repi > 1.7 * w
        rrup2[filt] = np.sqrt((repi[filt] - 1.7 * w[filt] / 2.0) ** 2 + zhyp[filt] ** 2)
        rrup = 0.5 * rrup1 + 0.5 * rrup2
        # Some assumptions here for rx. The rrup model is dominated by the site being updip from fault, so make rx positive.
        rx = -rjb / np.sqrt(2.0)
        rx[rjb == 0] = 0.5 * w[rjb == 0] * np.cos(d)
        # Calculate rx1 after computing rx
        rx1 = rx + w * np.cos(d)
        # Compute ry0 the same way as rx, but use fault length instead of width. Assume aspect ratio of 1.7
        ry0 = rjb / np.sqrt(2.0)
        ry0[rjb == 0] = 0.5 * 1.7 * w[rjb == 0] * np.cos(d)

        arrays = (
            m[filter],
            fault_type[filter],
            rate[filter],
            rjb[filter],
            rrup[filter],
            rx[filter],
            rx1[filter],
            ry0[filter],
            dip[filter],
            ztor[filter],
            zbor[filter],
        )
        return (arrays, {"index": np.flatnonzero(filter)}) if extras else arrays


def _cluster_extras(path, ruptures, index, source_info):
    """Cluster data of the ruptures of a cluster fault source model (see get_source_data).

    ruptures.npz has "cluster_id" (cluster of each rupture; -1 for independent ruptures) and
    "cluster_section" or "cluster_fault_id" (fault section, i.e. nshmp-lib rupture set, of
    each rupture within its cluster). clusters.npz has "cluster_id", "rate" (annual rate of the cluster), and "weight"
    (source logic tree weight of the cluster). With source_info "cluster_rupture_rate" =
    "conditional" (default), the "rate" of a cluster rupture is its weight within its fault
    (nshmp-lib stores the magnitude-variant weight in the rate field); with "absolute" it is
    the cluster rate times the cluster weight times that weight.
    """
    keys = ruptures.files
    if "cluster_id" not in keys:
        raise ValueError(f"cluster source model {path} has no cluster_id in ruptures.npz")
    cluster_id = np.asarray(ruptures["cluster_id"])[index].astype(np.int64)
    fault_key = next((k for k in ("cluster_fault_id", "cluster_section") if k in keys), None)
    if fault_key is not None:
        cluster_fault_id = np.asarray(ruptures[fault_key])[index].astype(np.int64)
    else:
        # every rupture is its own fault
        cluster_fault_id = index.astype(np.int64)
    clusters = np.load(str(path.joinpath("clusters.npz")))
    ids = np.asarray(clusters["cluster_id"]).astype(np.int64)
    c_rate = np.asarray(clusters["rate"], dtype=float)
    c_weight = (
        np.asarray(clusters["weight"], dtype=float)
        if "weight" in clusters.files
        else np.ones(len(ids))
    )
    order = np.argsort(ids)
    in_cluster = cluster_id >= 0
    pos = order[np.searchsorted(ids, cluster_id[in_cluster], sorter=order)]
    if not np.all(ids[pos] == cluster_id[in_cluster]):
        raise ValueError(f"ruptures.npz of {path} has cluster ids missing in clusters.npz")
    cluster_rate = np.zeros(len(index))
    cluster_rate[in_cluster] = c_rate[pos] * c_weight[pos]
    rate = np.asarray(ruptures["rate"], dtype=float)[index]
    if source_info.get("cluster_rupture_rate", "conditional") == "absolute":
        with np.errstate(divide="ignore", invalid="ignore"):
            cluster_weight = np.where(in_cluster, rate / cluster_rate, 0.0)
    else:
        cluster_weight = np.where(in_cluster, rate, 0.0)
    return {
        "cluster_id": cluster_id,
        "cluster_fault_id": cluster_fault_id,
        "cluster_weight": cluster_weight,
        "cluster_rate": cluster_rate,
    }


def _nshmp_point_source_data(path, p_xyz, dist_cutoff, m_min, ps):
    """Point source ruptures and distances with the nshmp-lib conventions.

    See geometry.point_source and tectonic_regions.normalize_point_source. The distance
    cutoff is applied to the horizontal distance from the site to the node (as nshmp-lib
    selects grid nodes within the maximum distance).

    Returns:
        (arrays, index): the get_source_data arrays, and the index of the rupture in
        ruptures.npz of each returned rupture
    """
    ruptures = np.load(str(path.joinpath("ruptures.npz")))
    rate = np.asarray(ruptures["rate"], dtype=float)
    m = np.asarray(ruptures["m"], dtype=float)
    style = np.asarray(ruptures["style"])
    node_index = ruptures["node_index"]
    nodes_index = np.load(str(path.joinpath("node_index.npy")))
    points = np.load(str(path.joinpath("points.npy")))
    sorter = np.argsort(nodes_index)
    pos = sorter[np.searchsorted(nodes_index, node_index, sorter=sorter)]
    node_lat, node_lon = point_source.xyz_to_latlon(points)
    site_lat, site_lon = point_source.xyz_to_latlon(p_xyz)
    site_lat, site_lon = site_lat[0], site_lon[0]

    # Depth to top of rupture
    index = np.arange(len(m))
    depth_source = ps["depth"]
    if depth_source == "auto":
        if "depth" in ruptures.files:
            depth_source = "rupture"
        elif os.path.exists(str(path.joinpath("depth.npy"))):
            depth_source = "node"
        elif ps.get("grid_depth_map"):
            depth_source = "depth_map"
        else:
            raise ValueError(f"{path}: no rupture depths (ruptures.npz depth, depth.npy, or grid_depth_map)")
    if depth_source == "rupture":
        if "depth" not in ruptures.files:
            raise ValueError(f"{path}: ruptures.npz has no depth")
        ztor = np.asarray(ruptures["depth"], dtype=float)
    elif depth_source == "node":
        ztor = np.load(str(path.joinpath("depth.npy")))[pos].astype(float)
    elif depth_source == "depth_map":
        index, ztor, depth_weight = point_source.expand_depths(m, ps["grid_depth_map"])
        rate = rate[index] * depth_weight
        m, style, pos = m[index], style[index], pos[index]
    else:
        raise ValueError(f'unknown point source depth "{ps["depth"]}"')

    # Magnitude filter first (cheap)
    if m_min is not None:
        keep = m >= m_min
        index, ztor, rate, m, style, pos = (
            a[keep] for a in (index, ztor, rate, m, style, pos)
        )

    r_node = point_source.horz_distance_fast(site_lat, site_lon, node_lat, node_lon)
    if dist_cutoff is not None:
        keep = r_node[pos] < dist_cutoff
        index, ztor, rate, m, style, pos = (
            a[keep] for a in (index, ztor, rate, m, style, pos)
        )

    bin_width = ps.get("distance_bin")
    smoothing = ps.get("smoothing") if bin_width else None
    if smoothing and len(pos):
        # Distribute nodes near the site (nshmp-lib grid optimization smoothing)
        used = np.unique(pos)
        sub_node, sub_lat, sub_lon, sub_scale = point_source.smooth_nodes(
            site_lat,
            site_lon,
            node_lat[used],
            node_lon[used],
            smoothing["density"],
            smoothing["limit"],
            smoothing["grid_spacing"],
        )
        sub_node = used[sub_node]
        order = np.argsort(sub_node, kind="stable")
        sub_node, sub_lat, sub_lon, sub_scale = (
            a[order] for a in (sub_node, sub_lat, sub_lon, sub_scale)
        )
        counts = np.bincount(sub_node, minlength=len(points))
        starts = np.searchsorted(sub_node, np.arange(len(points)))
        n_sub = counts[pos]
        rup = np.repeat(np.arange(len(pos)), n_sub)
        offset = np.arange(len(rup)) - np.repeat(np.cumsum(n_sub) - n_sub, n_sub)
        sub = starts[pos][rup] + offset
        index, ztor, rate, m, style = (a[rup] for a in (index, ztor, rate, m, style))
        rate = rate * sub_scale[sub]
        r_h = point_source.horz_distance_fast(
            site_lat, site_lon, sub_lat[sub], sub_lon[sub]
        )
    else:
        r_h = r_node[pos]

    if bin_width and len(r_h):
        # nshmp-lib grid optimization: horizontal distances at the centers of distance bins,
        # and ruptures with the same magnitude, mechanism, depth, and distance bin combined
        i_bin = np.floor(r_h / bin_width).astype(np.int64)
        r_h = (i_bin + 0.5) * bin_width
        key = np.stack(
            [
                np.round(m * 1000.0).astype(np.int64),
                style.astype(np.int64),
                np.round(ztor * 1000.0).astype(np.int64),
                i_bin,
            ],
            axis=1,
        )
        _, first, inverse = np.unique(
            key, axis=0, return_index=True, return_inverse=True
        )
        rate = np.bincount(inverse.ravel(), weights=rate)
        index, ztor, m, style, r_h = (a[first] for a in (index, ztor, m, style, r_h))

    max_depth = ps.get("max_depth")
    max_width = None if max_depth is not None else ps.get("max_width")
    if ps["type"] == "fixed_strike":
        # nshmp-lib zones: FIXED_STRIKE sources at nodes with a strike (strike.npy), FINITE
        # sources (with the rJB correction) at nodes without one (NaN). Not optimized.
        strike = np.load(str(path.joinpath("strike.npy")))[pos].astype(float)
        fixed = np.isfinite(strike)
        parts = []
        if np.any(fixed):
            d = point_source.fixed_strike_distances(
                m[fixed], style[fixed], node_lat[pos][fixed], node_lon[pos][fixed],
                strike[fixed], site_lat, site_lon, ztor[fixed], ps["rupture_scaling"],
                max_depth=max_depth, max_width=max_width,
            )
            d["index"] = np.flatnonzero(fixed)[d["index"]]
            parts.append(d)
        if np.any(~fixed):
            sel = np.flatnonzero(~fixed)
            d = point_source.finite_point_source_distances(
                m[sel], style[sel], r_h[sel], ztor[sel], ps["rupture_scaling"], "finite",
                max_depth=max_depth, max_width=max_width,
            )
            d["index"] = sel[d["index"]]
            parts.append(d)
        if parts:
            d = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
        else:
            d = point_source.finite_point_source_distances(
                m, style, r_h, ztor, ps["rupture_scaling"], "point"
            )
    else:
        d = point_source.finite_point_source_distances(
            m,
            style,
            r_h,
            ztor,
            ps["rupture_scaling"],
            ps["type"],
            max_depth=max_depth,
            max_width=max_width,
        )
    k = d["index"]
    arrays = (
        m[k],
        np.asarray(style[k]),
        rate[k] * d["rate_scale"],
        d["rjb"],
        d["rrup"],
        d["rx"],
        d["rx1"],
        d["ry0"],
        d["dip"],
        d["ztor"],
        d["zbor"],
    )
    return arrays, index[k]


def get_ground_motion_data(
    gmm,
    vs30,
    measured_vs30,
    z1p0,
    z2p5,
    fault_type,
    rjb,
    rrup,
    rx,
    rx1,
    ry0,
    m,
    ztor,
    zbor,
    dip,
):
    """Computes arrays containing mean and standard deviation of the natural log of a ground
    motion intensity measure.

    Inputs:
        gmm (string): Ground motion model. One of "ask14", "bssa14", "cb14", "cy14", "idriss14"
        vs30 (float): Time-averaged shear wave velocity in the upper 30m in m/s
        measured_vs30 (bool): boolean field indicating whether vs30 is measured (True) or inferred (False)
        z1p0 (float): Isosurface depth to a shear wave velocity of 1.0 km/s in km, length = N
        z2p5 (float): Isosurface depth to a shear wave velocity of 2.5 km/s in km, length = N
        fault_type (Numpy array, dtype=int): Numpy array of fault type. 1 = reverse, 2 = normal, 3 = strike slip, length = N
        rjb (Numpy array, dtype=float): Numpy array of Joyner-Boore distances between the site and each rupture, length = N
        rrup (Numpy array, dtype=float): Numpy array of rupture distance between the site and each rupture, length = N
        rx (Numpy array, dtype=float): Numpy array of distance from site to the surface projection of the top of each
            rupture, measured perpendicular to the strike, length = N
        rx1 (Numpy array, dtype=float): Numpy array of distance from site to the surface projection of the bottom of each
            rupture, measured perpendicular to the strike, length = N
        ry0 (Numpy array, dtype=float): Numpy array of the distance from the site to the surface projection of each
            rupture, measured parallel to the strike, length = N
        m (Numpy array, dtype=float): Numpy array of magnitudes, length = N
        ztor (Numpy array, dtype=float): Numpy array of depth to the top of each rupture in km, length = N
        zbor (Numpy array, dtype=float): Numpy array of depth to the bottom of each rupture in km, length = N
        dip (Numpy array, dtype=float): Numpy array of dip angle for each rupture in degrees, length = N

    Returns:
        mu_ln_pga (Numpy array, dtype=float): array of the mean of the natural logs of the ground motion intensity measure, IM
        sigma_ln_pga (Numpy array, dtype=float): array of the standard deviation of the natural logs of the IM

    Note:
        N = number of events
    """
    if gmm == "bssa14":
        mu_ln_pga, sigma_ln_pga = bssa14.get_im(vs30, rjb, m, fault_type)
    elif gmm == "cb14":
        mu_ln_pga, sigma_ln_pga = cb14.get_im(
            vs30, rjb, rrup, rx, rx1, m, fault_type, ztor, zbor, dip, z2p5=z2p5
        )
    elif gmm == "cy14":
        mu_ln_pga, sigma_ln_pga = cy14.get_im(
            vs30, rjb, rrup, rx, m, fault_type, measured_vs30, dip, ztor, z1p0=z1p0
        )
    elif gmm == "ask14":
        mu_ln_pga, sigma_ln_pga = ask14.get_im(
            vs30, rrup, rx, rx1, ry0, m, fault_type, measured_vs30, dip, ztor, z1p0=z1p0
        )
    elif gmm == "idriss14":
        mu_ln_pga, sigma_ln_pga = idriss14.get_im(vs30, rrup, m, fault_type)
    else:
        raise ValueError(
            f'incorrect ground motion model "{gmm}", '
            'expected one of "ask14", "bssa14", "cb14", "cy14", "idriss14"'
        )
    return [mu_ln_pga, sigma_ln_pga]


def get_liquefaction_cdfs(m, mu_ln_pga, sigma_ln_pga, fsl, liquefaction_model, config):
    """Computes log-normal cumulative distribution functions for factor of safety of liquefaction

    Inputs:
        m (array, dtype=float): Numpy array of magnitudes, length = N
        mu_ln_pga (array, dtype=float): Numpy array of the mean of the natural logs of the earthquake
            ground motion intensity measure, length = N
        sigma_ln_pga (array, dtype=float): Numpy array of the standard deviation of the natural logs
            of the ground motion intensity measure, length = N
        fsl (array, dtype=float): Numpy array of factor of safety values at which to compute the
            liquefaction hazard curve, length = L
        liquefaction_model (string): Liquefaction model to use. One of "cetin_et_al_2018", "moss_et_al_2006",
            "boulanger_idriss_2016", "boulanger_idriss_2012", "ngl_smt_2024".
        config (dict): A Python dictionary read from the config file

    Returns:
        fsl_cdfs (Numpy ndarray, dtype=float): Numpy array of cumulative distribution functions representing
            down-crossing rate of fsl for each event, length = N x L
        eps (Numpy ndarray, dtype=float): Numpy array of epsilon values, representing number of standard deviations
            of the natural log of fsl relative to the mean of the natural log of fsl, length = N x L

    Notes:
        N = number of events
        L = number of fsl values at which to compute hazard
    """
    if liquefaction_model == "cetin_et_al_2018":
        c = config["liquefaction_models"]["cetin_et_al_2018"]
        return cetin_et_al_2018.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            c["sigmav"],
            c["sigmavp"],
            c["vs12"],
            c["depth"],
            c["n160"],
            c["fc"],
            fsl,
            c["pa"],
        )
    elif liquefaction_model == "moss_et_al_2006":
        c = config["liquefaction_models"]["moss_et_al_2006"]
        return moss_et_al_2006.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            c["sigmav"],
            c["sigmavp"],
            c["depth"],
            c["qc"],
            c["fs"],
            fsl,
            c["pa"],
        )
    elif liquefaction_model == "boulanger_idriss_2016":
        c = config["liquefaction_models"]["boulanger_idriss_2016"]
        return boulanger_idriss_2016.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            c["sigmav"],
            c["sigmavp"],
            c["depth"],
            c["qc1ncs"],
            fsl,
            c["pa"],
        )
    elif liquefaction_model == "boulanger_idriss_2012":
        c = config["liquefaction_models"]["boulanger_idriss_2012"]
        return boulanger_idriss_2012.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            c["sigmav"],
            c["sigmavp"],
            c["depth"],
            c["n160"],
            c["fc"],
            fsl,
            c["pa"],
        )
    elif liquefaction_model == "ngl_smt_2024":
        c = config["liquefaction_models"]["ngl_smt_2024"]
        return ngl_smt_2024.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            np.asarray(c["ztop"], dtype=float),
            np.asarray(c["zbot"], dtype=float),
            np.asarray(c["qc1ncs"], dtype=float),
            np.asarray(c["ic"], dtype=float),
            np.asarray(c["sigmav"], dtype=float),
            np.asarray(c["sigmavp"], dtype=float),
            np.asarray(c["ksat"], dtype=float),
            fsl,
            float(c.get("pa", 101.325)),
        )


def get_disagg(hazards, m, r, eps, m_bin_edges, r_bin_edges, eps_bin_edges):
    """Computes disaggregation for seismic hazard and/or liquefaction hazard

    Inputs:
        hazards (Numpy ndarray, dtype = float) = array of hazard values, shape = L x N.
        m (Numpy array, dtype = float) = array of magnitude values, length = N
        r (Numpy array, dtype = float) = array of distance values, length = N
        eps (Numpy ndarray, dtype = float) = array of epsilon values, shape = L x N
        m_bin_edges (Numpy array, dtype = float) = array defining magnitude bin edges, length = M + 1
        r_bin_edges (Numpy array, dtype = float) = array defining distance bin edges, length = R + 1
        eps_bin_edges (Numpy array, dtype = float) = array defining epsilon bin edges, length = E + 1

    Returns:
        disagg (Numpy ndarray, dtype=float) = Numpy array of contribution to hazard within each
            magnitude, distance, and epsilon bin for each pga (PSHA) or fsl (PLHA) value, shape = L x M x R x E

    Notes:
        N = number of events
        L = number of intensity measure values
        M = number of magnitude bins (note that m_bin_edges has a length of M + 1)
        R = number of distance bins (note that r_bin_edges has a length of R + 1)
        E = number of epsilon bins (note that eps_bin_edges has a length of E + 1)
    """
    # use Numpy digitize function to assign bin numbers
    m_hazard = np.digitize(m, m_bin_edges)
    r_hazard = np.digitize(r, r_bin_edges)
    eps_hazard = np.digitize(eps, eps_bin_edges)

    # compute number of bins per intensity measure value
    Nbins = (len(m_bin_edges) - 1) * (len(r_bin_edges) - 1) * (len(eps_bin_edges) - 1)

    # use tensor rank reduction to create L x Nbins length array of indices
    bin_indices = (
        eps_hazard
        - 1
        + (r_hazard - 1) * (len(eps_bin_edges) - 1)
        + (m_hazard - 1) * (len(eps_bin_edges) - 1) * (len(r_bin_edges) - 1)
    )

    # create empty array to store hazard sums, and use Numpy bincount to efficiently sum hazards
    # within each bin. The np.bincount function must be performed on a 1D array, so we need to loop
    # over the pga values. The cost of that loop is minimal since the number of pga values
    # is generally small
    sum_hazard = np.empty((len(bin_indices), Nbins))
    for i in range(len(bin_indices)):
        sum_hazard[i] = np.bincount(bin_indices[i], weights=hazards[i], minlength=Nbins)

    # reshape the sum_hazard array
    disagg = sum_hazard.reshape(
        (
            len(bin_indices),
            len(m_bin_edges) - 1,
            len(r_bin_edges) - 1,
            len(eps_bin_edges) - 1,
        )
    )

    return disagg


def get_exceedance(pga, mu_ln_pga, sigma_ln_pga, truncation_level=None):
    """Probability that each ground motion value is exceeded for each event.

    Inputs:
        pga (Numpy array, dtype=float): ground motion values (g), length = L
        mu_ln_pga (Numpy array, dtype=float): mean of ln(PGA) of each event, length = N
        sigma_ln_pga (Numpy array, dtype=float): standard deviation of ln(PGA), length = N
        truncation_level (float): if given, the log-normal distribution is truncated at
            mu + truncation_level * sigma and renormalized, as in the nshmp-lib
            TRUNCATION_UPPER_ONLY exceedance model (USGS NSHM: 3). None (default) is untruncated.

    Returns:
        eps (Numpy ndarray): epsilon of each ground motion value and event, shape = L x N
        p (Numpy ndarray): probability of exceedance, shape = L x N
    """
    eps = (np.log(pga[:, np.newaxis]) - mu_ln_pga) / sigma_ln_pga
    p = 1 - ndtr(eps)
    if truncation_level is not None:
        p_hi = ndtr(-float(truncation_level))
        p = np.clip((p - p_hi) / (1.0 - p_hi), 0.0, 1.0)
    return eps, p


def get_cluster_hazard(p, weight, cluster_index, fault_index, cluster_rate):
    """Hazard of cluster sources (nshmp-lib ClusterRuptureSet hazard).

    In nshmp-lib (calc/Transforms.ClusterGroundMotionsToCurves and
    ExceedanceModel.clusterExceedance), a cluster is a set of faults that rupture together
    with the rate of the cluster. For each ground motion model, the probability that a fault
    of the cluster produces an exceedance is the weighted sum over its magnitude variants,
    P_f = sum_i w_i P_i; the probability that the cluster event produces an exceedance is
    P_c = 1 - prod_f (1 - P_f); and the hazard is sum_c rate_c P_c, where rate_c includes the
    logic tree weight of the cluster. This must be evaluated separately for each ground
    motion model, before the ground motion model weights are applied.

    Inputs:
        p (Numpy ndarray): conditional probability of exceedance (PSHA) or of non-exceedance of
            the factor of safety (PLHA) for each cluster rupture, shape = L x N
        weight (Numpy array): weight w_i of each rupture within its fault, length = N
        cluster_index (Numpy array, dtype=int): cluster of each rupture, length = N
        fault_index (Numpy array, dtype=int): fault of each rupture within its cluster, length = N
        cluster_rate (Numpy array): rate (times logic tree weight) of the cluster of each
            rupture, length = N

    Returns:
        hazard (Numpy array): annual rate, length = L
        contributions (Numpy ndarray): hazard attributed to each rupture for disaggregation,
            shape = L x N. The hazard of a cluster is shared among its ruptures in proportion
            to w_i P_i, so the contributions sum to the hazard.
    """
    p = np.atleast_2d(p)
    n_values = p.shape[0]
    if p.shape[1] == 0:
        return np.zeros(n_values), np.zeros_like(p)
    _, c = np.unique(cluster_index, return_inverse=True)
    c = c.ravel()
    n_clusters = c.max() + 1
    _, g = np.unique(np.stack([c, fault_index], axis=1), axis=0, return_inverse=True)
    g = g.ravel()
    n_groups = g.max() + 1
    group_cluster = np.zeros(n_groups, dtype=int)
    group_cluster[g] = c
    rate_c = np.zeros(n_clusters)
    rate_c[c] = cluster_rate
    wp = weight * p
    p_fault = np.empty((n_values, n_groups))
    p_sum = np.empty((n_values, n_clusters))
    log_survival = np.empty((n_values, n_clusters))
    for i in range(n_values):
        p_fault[i] = np.bincount(g, weights=wp[i], minlength=n_groups)
        p_sum[i] = np.bincount(c, weights=wp[i], minlength=n_clusters)
    p_fault = np.clip(p_fault, 0.0, 1.0)
    with np.errstate(divide="ignore"):
        log_fault = np.log1p(-p_fault)
    for i in range(n_values):
        log_survival[i] = np.bincount(group_cluster, weights=log_fault[i], minlength=n_clusters)
    cluster_hazard = -np.expm1(log_survival) * rate_c
    with np.errstate(divide="ignore", invalid="ignore"):
        share = np.where(p_sum[:, c] > 0.0, wp / p_sum[:, c], 0.0)
    contributions = share * cluster_hazard[:, c]
    return cluster_hazard.sum(axis=1), contributions


def nshm_site_parameters(site_config, regions_used=()):
    """zSed and location-dependent ground motion model trees of a site from the NSHM site data.

    The USGS NSHM (nshmp-lib SiteData, used by the USGS hazard web service for every site) gives
    a site on the Gulf and Atlantic coastal plain its sediment thickness zSed (Boyd, 2023) and,
    inside the "Coastal Plain CPA region", a stable crust ground motion model logic tree with
    the Chapman and Guo (2021) coastal plain amplification models (see
    :mod:`ucla_plha.nshm_site_data`).

    Args:
        site_config (dict): config["site"]. "zsed" (km, or null for a site off the coastal
            plain) overrides the NSHM value. "nshm_site_data": false turns off both the zSed
            look-up and the gmm-region trees (default true).
        regions_used (iterable): tectonic regions of the source models that are used

    Returns:
        (zsed, site_trees, notes): zsed in km or None; site_trees is {region: (description,
        {Gmm id: weight})} for :func:`tectonic_regions.parse_ground_motion_models`
    """
    use = site_config.get("nshm_site_data", True)
    lon, lat = site_config["longitude"], site_config["latitude"]
    notes = []
    if "zsed" in site_config:
        zsed = site_config["zsed"]
    elif use:
        zsed = nshm_site_data.coastal_plain_zsed(lon, lat)
        if zsed is not None and "stable_crust" in regions_used:
            notes.append(
                f"site: zsed = {zsed:g} km (coastal plain sediment thickness of the USGS NSHM, "
                "Boyd 2023)"
            )
    else:
        zsed = None
    site_trees = {}
    if use and "stable_crust" in regions_used and nshm_site_data.in_coastal_plain_region(lon, lat):
        site_trees["stable_crust"] = (
            f'NSHM gmm-region "{nshm_site_data.region_name()}"',
            nshm_site_data.coastal_plain_stable_crust_tree(),
        )
    return zsed, site_trees, notes


def get_hazard(config_file):
    """Reads config file and runs PSHA and PLHA

    Inputs:
        config_file (string): Filename, including path, of config file. Must follow the schema
            defined by ucla_plha_schema.json. See documentation for more thorough documentation
            of the config file.

    Returns:
        output (dict): Python dictionary containing output of analysis. The output contains all
            of the inputs for preservation, along with the hazard curve(s) and any requested
            disaggregation data. See documentation for more thorough description of output.

    Notes:
        Every source model has a tectonic region (source_info.json in its directory; active
        crust if there is none). Ground motion models are given per tectonic region in the
        config file ("ground_motion_models": {"active_crust": {...}, "stable_crust": {...},
        ...}), or, in the earlier format, as a single set of models that is used for the active
        crust. Each source model uses the ground motion models of its region (the default USGS
        NSHM models if the region is not in the config file). Weights are normalized within
        each tectonic region, and within each group of alternative source models (source type,
        tectonic region, and NSHM component).
    """
    # Validate config_file against schema. If ngl_smt_2024 liquefaction model is used, the cpt_data file is
    # validated in the get_liquefaction_hazards function.

    schema = json.loads(
        open(files("ucla_plha").joinpath("ucla_plha_schema.json")).read().lower()
    )
    config = json.loads(open(config_file).read().lower())

    # validate config file and return messageg if errors are encountered
    try:
        jsonschema.validate(config, schema)
    except jsonschema.ValidationError as e:
        print("Config File Error:", e.message)
        return

    # Source model information, and normalization of source model weights within each group
    # of alternative source models
    source_infos = {}
    group_sums = {}
    for source_type, models in config["source_models"].items():
        for source_model, entry in models.items():
            info = get_source_info(source_type, source_model)
            source_infos[(source_type, source_model)] = info
            group = tectonic_regions.weight_group(info, source_type)
            group_sums[group] = group_sums.get(group, 0.0) + entry.get("weight", 0.0)
    for (source_type, source_model), info in source_infos.items():
        group = tectonic_regions.weight_group(info, source_type)
        if group_sums[group] > 0:
            config["source_models"][source_type][source_model]["weight"] /= group_sums[
                group
            ]

    # Ground motion model logic trees of the tectonic regions, weights normalized within
    # each region
    regions_used = {
        info["tectonic_region"]
        for key, info in source_infos.items()
        if config["source_models"][key[0]][key[1]]["weight"] > 0
    }
    # Location-dependent NSHM site data (nshmp-lib SiteData, as used by the USGS hazard
    # service): the coastal plain sediment thickness zSed and the stable crust ground motion
    # models of the Coastal Plain CPA region
    site_zsed, site_trees, site_notes = nshm_site_parameters(config["site"], regions_used)
    used_default = set()
    gmm_trees, notes = tectonic_regions.parse_ground_motion_models(
        config, regions_used, site_trees, used_default
    )
    notes = site_notes + notes
    # Ground motion models of individual source models: the "ground_motion_models" of the
    # source model in the config, or a "gmm_tree" list in its source_info.json (e.g. the NSHM
    # system grid in the stable crust, which mixes NGA-East and NGA-West2 models). As in
    # nshmp-lib, the tree of a site's NSHM gmm-region replaces the source_info.json trees of
    # its tectonic region (when the region uses the default tree).
    model_trees = {}
    for (source_type, source_model), info in source_infos.items():
        entry = config["source_models"][source_type][source_model]
        if entry["weight"] <= 0:
            continue
        if (
            not isinstance(entry.get("ground_motion_models"), dict)
            and info["tectonic_region"] in site_trees
            and info["tectonic_region"] in used_default
        ):
            continue
        tree = tectonic_regions.source_model_gmm_entries(entry, info)
        if tree is not None:
            model_trees[(source_type, source_model)] = tectonic_regions.build_region_tree(
                source_model, tree, notes
            )

    liquefaction_model_weight_sum = 0.0
    liquefaction_models = [
        "boulanger_idriss_2012",
        "boulanger_idriss_2016",
        "cetin_et_al_2018",
        "moss_et_al_2006",
        "ngl_smt_2024",
    ]
    for liquefaction_model in liquefaction_models:
        liquefaction_model_weight_sum += (
            config.get("liquefaction_models", {})
            .get(liquefaction_model, {})
            .get("weight", 0.0)
        )
    if liquefaction_model_weight_sum > 0:
        for liquefaction_model in liquefaction_models:
            if config.get("liquefaction_models", {}).get(liquefaction_model, {}):
                config["liquefaction_models"][liquefaction_model][
                    "weight"
                ] /= liquefaction_model_weight_sum

    # Read site properties
    latitude = config["site"]["latitude"]
    longitude = config["site"]["longitude"]
    elevation = config["site"]["elevation"]
    point = np.asarray([latitude, longitude, elevation])
    p_xyz = geometry.point_to_xyz(point)
    vs30 = config["site"]["vs30"]
    measured_vs30 = config["site"].get("measured_vs30", False)
    z1p0 = config["site"].get("z1p0", None)
    z2p5 = config["site"].get("z2p5", None)
    site = {
        "vs30": vs30,
        "measured_vs30": measured_vs30,
        "z1p0": z1p0,
        "z2p5": z2p5,
        "zsed": site_zsed,
    }
    dist_cutoff = config.get("constraints", {}).get("dist_cutoff", None)
    m_min = config.get("constraints", {}).get("m_min", None)
    truncation_level = config.get("constraints", {}).get("truncation_level", None)
    rupture_rx = config.get("constraints", {}).get("rupture_rx", "closest_section")
    # "source_info" (default): the settings of each source model's source_info.json
    fault_distances = config.get("constraints", {}).get("fault_distances", "source_info")
    fault_distances = None if fault_distances == "source_info" else fault_distances
    section_rx = config.get("constraints", {}).get("section_rx", "source_info")
    section_rx = None if section_rx == "source_info" else section_rx

    # Read output properties
    if "psha" in config["output"].keys():
        pga = np.asarray(config["output"]["psha"]["pga"], dtype=float)
        output_psha = True
        output_source_hazard = config["output"]["psha"].get(
            "source_model_hazard", False
        )
        source_hazard = {}
        if "disaggregation" in config["output"]["psha"].keys():
            output_psha_disaggregation = True
            psha_magnitude_bin_edges = np.asarray(
                config["output"]["psha"]["disaggregation"]["magnitude_bin_edges"],
                dtype=float,
            )
            psha_distance_bin_edges = np.asarray(
                config["output"]["psha"]["disaggregation"]["distance_bin_edges"],
                dtype=float,
            )
            psha_epsilon_bin_edges = np.asarray(
                config["output"]["psha"]["disaggregation"]["epsilon_bin_edges"],
                dtype=float,
            )
            psha_magnitude_bin_center = 0.5 * (
                psha_magnitude_bin_edges[0:-1] + psha_magnitude_bin_edges[1:]
            )
            psha_distance_bin_center = 0.5 * (
                psha_distance_bin_edges[0:-1] + psha_distance_bin_edges[1:]
            )
            psha_epsilon_bin_center = 0.5 * (
                psha_epsilon_bin_edges[0:-1] + psha_epsilon_bin_edges[1:]
            )
            psha_disagg = np.zeros(
                (
                    len(pga),
                    len(psha_magnitude_bin_center),
                    len(psha_distance_bin_center),
                    len(psha_epsilon_bin_center),
                ),
                dtype=float,
            )
        else:
            output_psha_disaggregation = False
        seismic_hazard = np.zeros(len(pga))
    else:
        output_psha = False
        output_psha_disaggregation = False
        output_source_hazard = False

    if "plha" in config["output"].keys():
        fsl = np.asarray(config["output"]["plha"]["fsl"], dtype=float)
        output_plha = True
        if "disaggregation" in config["output"]["plha"].keys():
            output_plha_disaggregation = True
            plha_magnitude_bin_edges = np.asarray(
                config["output"]["plha"]["disaggregation"]["magnitude_bin_edges"],
                dtype=float,
            )
            plha_distance_bin_edges = np.asarray(
                config["output"]["plha"]["disaggregation"]["distance_bin_edges"],
                dtype=float,
            )
            plha_epsilon_bin_edges = np.asarray(
                config["output"]["plha"]["disaggregation"]["epsilon_bin_edges"],
                dtype=float,
            )
            plha_magnitude_bin_center = 0.5 * (
                plha_magnitude_bin_edges[0:-1] + plha_magnitude_bin_edges[1:]
            )
            plha_distance_bin_center = 0.5 * (
                plha_distance_bin_edges[0:-1] + plha_distance_bin_edges[1:]
            )
            plha_epsilon_bin_center = 0.5 * (
                plha_epsilon_bin_edges[0:-1] + plha_epsilon_bin_edges[1:]
            )
            plha_disagg = np.zeros(
                (
                    len(fsl),
                    len(plha_magnitude_bin_center),
                    len(plha_distance_bin_center),
                    len(plha_epsilon_bin_center),
                ),
                dtype=float,
            )
        else:
            output_plha_disaggregation = False
        liquefaction_hazard = np.zeros(len(fsl))
    else:
        output_plha = False
        output_plha_disaggregation = False

    # Loop over source models. We have fault_source_models and point_source_models, so there are two loops
    for source_model in config["source_models"].keys():
        for fault_source_model in config["source_models"][source_model].keys():
            source_model_weight = config["source_models"][source_model][
                fault_source_model
            ]["weight"]
            # move on to next source model if weight is less than or equal to zero
            if source_model_weight <= 0:
                continue
            info = source_infos[(source_model, fault_source_model)]
            region = info["tectonic_region"]
            branches = model_trees.get(
                (source_model, fault_source_model), gmm_trees.get(region)
            )
            # all ground motion models of the region determine the distance types
            gmms = [branch.key for branch in branches]
            region_cutoff = config["source_models"][source_model][
                fault_source_model
            ].get("dist_cutoff", tectonic_regions.region_value(dist_cutoff, region))
            cluster = info["cluster"]
            if cluster:
                data, extras = get_source_data(
                    source_model,
                    fault_source_model,
                    p_xyz,
                    region_cutoff,
                    m_min,
                    gmms,
                    extras=True,
                    source_info=info,
                    rupture_rx=rupture_rx,
                    fault_distances=fault_distances,
                    section_rx=section_rx,
                )
                in_cluster = extras["cluster_id"] >= 0
            else:
                data = get_source_data(
                    source_model,
                    fault_source_model,
                    p_xyz,
                    region_cutoff,
                    m_min,
                    gmms,
                    rupture_rx=rupture_rx,
                    fault_distances=fault_distances,
                    section_rx=section_rx,
                )
            m, fault_type, rate, rjb, rrup, rx, rx1, ry0, dip, ztor, zbor = data
            # distance used for disaggregation: rjb, or rrup if no model uses rjb
            r_disagg = rjb if "rjb" in _distance_needs(gmms) or not len(gmms) else rrup
            rupture = {
                "m": m,
                "fault_type": fault_type,
                "rjb": rjb,
                "rrup": rrup,
                "rx": rx,
                "ry0": ry0,
                "dip": dip,
                "ztor": ztor,
                "zbor": zbor,
            }
            source_key = fault_source_model
            # Loop over ground motion models.
            for branch in branches:
                # move on to next ground motion model if weight is less than or equal to zero
                if branch.weight <= 0:
                    continue
                if branch.spec.kind == "native":
                    gm_branches = [(1.0, *get_ground_motion_data(
                        branch.spec.name,
                        vs30,
                        measured_vs30,
                        z1p0,
                        z2p5,
                        fault_type,
                        rjb,
                        rrup,
                        rx,
                        rx1,
                        ry0,
                        m,
                        ztor,
                        zbor,
                        dip,
                    ))]
                else:
                    gm_branches = pygmm_gmms.get_ground_motion_branches(
                        branch.spec, rupture, site, region
                    )
                # Loop over the branches of the ground motion distribution of the model (e.g. the
                # USGS epistemic branches of the nshmp-lib models). As in nshmp-lib
                # (ExceedanceModel.treeExceedanceCombined), the exceedance probabilities of the
                # branches are weighted and summed (the collapsed median is not the mixture).
                for gm_branch_weight, mu_ln_pga, sigma_ln_pga in gm_branches:
                    ground_motion_model_weight = branch.weight * gm_branch_weight
                    # Compute seismic hazard if requested in config file
                    if output_psha:
                        if truncation_level is None and not cluster:
                            eps = (np.log(pga[:, np.newaxis]) - mu_ln_pga) / sigma_ln_pga
                            seismic_hazards = (1 - ndtr(eps)) * rate
                        else:
                            eps, p = get_exceedance(
                                pga, mu_ln_pga, sigma_ln_pga, truncation_level
                            )
                            seismic_hazards = p * rate
                        if cluster:
                            seismic_hazards[:, in_cluster] = 0.0
                            cluster_curve, contributions = get_cluster_hazard(
                                p[:, in_cluster],
                                extras["cluster_weight"][in_cluster],
                                extras["cluster_id"][in_cluster],
                                extras["cluster_fault_id"][in_cluster],
                                extras["cluster_rate"][in_cluster],
                            )
                            if output_source_hazard:
                                key = source_key + ":cluster"
                                source_hazard[key] = source_hazard.get(
                                    key, 0.0
                                ) + source_model_weight * ground_motion_model_weight * (
                                    cluster_curve
                                )
                        curve = np.sum(seismic_hazards, axis=1)
                        if output_source_hazard:
                            source_hazard[source_key] = (
                                source_hazard.get(source_key, 0.0)
                                + source_model_weight * ground_motion_model_weight * curve
                            )
                        if cluster:
                            seismic_hazards[:, in_cluster] = contributions
                            curve = curve + cluster_curve
                        seismic_hazard += (
                            source_model_weight * ground_motion_model_weight * curve
                        )
                        # Compute seismic hazard disaggregation if requested in config file
                        if output_psha_disaggregation:
                            psha_disagg += (
                                source_model_weight
                                * ground_motion_model_weight
                                * get_disagg(
                                    seismic_hazards,
                                    m,
                                    r_disagg,
                                    eps,
                                    psha_magnitude_bin_edges,
                                    psha_distance_bin_edges,
                                    psha_epsilon_bin_edges,
                                )
                            )
                    # Compute liquefaction hazard if requested in config file
                    if "liquefaction_models" in config.keys():
                        for liquefaction_model in config["liquefaction_models"].keys():
                            liquefaction_model_weight = config["liquefaction_models"][
                                liquefaction_model
                            ]["weight"]
                            if liquefaction_model_weight <= 0:
                                continue
                            if output_plha:
                                liquefaction_hazards, eps = get_liquefaction_cdfs(
                                    m,
                                    mu_ln_pga,
                                    sigma_ln_pga,
                                    fsl,
                                    liquefaction_model,
                                    config,
                                )
                                if cluster:
                                    p_liq = liquefaction_hazards[in_cluster].T
                                liquefaction_hazards *= rate[:, np.newaxis]
                                if cluster:
                                    liquefaction_hazards[in_cluster] = 0.0
                                    cluster_curve, contributions = get_cluster_hazard(
                                        p_liq,
                                        extras["cluster_weight"][in_cluster],
                                        extras["cluster_id"][in_cluster],
                                        extras["cluster_fault_id"][in_cluster],
                                        extras["cluster_rate"][in_cluster],
                                    )
                                    liquefaction_hazards[in_cluster] = contributions.T
                                liquefaction_hazard += (
                                    source_model_weight
                                    * ground_motion_model_weight
                                    * liquefaction_model_weight
                                    * np.sum(liquefaction_hazards, axis=0)
                                )
                                eps = eps.T
                                liquefaction_hazards = liquefaction_hazards.T
                                # Compute liquefaction hazard disaggregation if requested in config file
                                if output_plha_disaggregation:
                                    plha_disagg += (
                                        source_model_weight
                                        * ground_motion_model_weight
                                        * liquefaction_model_weight
                                        * get_disagg(
                                            liquefaction_hazards,
                                            m,
                                            r_disagg,
                                            eps,
                                            plha_magnitude_bin_edges,
                                            plha_distance_bin_edges,
                                            plha_epsilon_bin_edges,
                                        )
                                    )
    # Now prepare output
    output = {}
    output["input"] = config
    output["output"] = {}
    if output_psha:
        if output_psha_disaggregation:
            for i in range(len(pga)):
                psha_disagg[i] = psha_disagg[i] / seismic_hazard[i] * 100.0
            output["output"]["psha"] = {
                "PGA": pga.tolist(),
                "annual_rate_of_exceedance": seismic_hazard.tolist(),
                "disaggregation": psha_disagg.tolist(),
            }
        else:
            output["output"]["psha"] = {
                "PGA": pga.tolist(),
                "annual_rate_of_exceedance": seismic_hazard.tolist(),
            }
        if output_source_hazard:
            output["output"]["psha"]["source_model_hazard"] = {}
            for key, curve in source_hazard.items():
                name, _, part = key.partition(":")
                info = next(i for k, i in source_infos.items() if k[1] == name)
                component = (
                    info["cluster_nshm_component"] if part else info["nshm_component"]
                )
                output["output"]["psha"]["source_model_hazard"][key] = {
                    "tectonic_region": info["tectonic_region"],
                    "nshm_component": component,
                    "annual_rate_of_exceedance": np.asarray(curve).tolist(),
                }
    if output_plha:
        if output_plha_disaggregation:
            for i in range(len(fsl)):
                plha_disagg[i] = plha_disagg[i] / liquefaction_hazard[i] * 100.0
            output["output"]["plha"] = {
                "FSL": fsl.tolist(),
                "annual_rate_of_nonexceedance": liquefaction_hazard.tolist(),
                "disaggregation": plha_disagg.tolist(),
            }
        else:
            output["output"]["plha"] = {
                "FSL": fsl.tolist(),
                "annual_rate_of_nonexceedance": liquefaction_hazard.tolist(),
            }
    if notes:
        output["notes"] = notes

    if "outputfile" in config["output"].keys():
        if config["output"]["outputfile"] == "default":
            outputfilename = config_file.split(".json")[0] + "_output.json"
        else:
            outputfilename = config["output"]["outputfile"]
        with open(outputfilename, "w") as outputfile:
            json.dump(output, outputfile, indent=4)

    return output

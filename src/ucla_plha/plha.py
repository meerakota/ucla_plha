import json
import os
import warnings
from importlib.resources import files

import jsonschema
import numpy as np
from scipy.special import ndtr, ndtri

from ucla_plha import gc2_cache, pga_grid
from ucla_plha.liquefaction_models import (
    boulanger_idriss_2012,
    cetin_et_al_2018,
    ngl_smt_2024,
    moss_et_al_2006,
    boulanger_idriss_2016,
)
from ucla_plha.ground_motion_models import ask14, bssa14, cb14, cy14
from ucla_plha.geometry import geometry


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


class _NpzArrays(dict):
    """Contents of an npz file (all members read), with the member names in .files as NpzFile."""

    @property
    def files(self):
        return list(self.keys())


def _read_only(value):
    """Marks the numpy arrays of value (an array, or a dict or tuple of them) read-only."""
    if isinstance(value, np.ndarray):
        value.flags.writeable = False
    elif isinstance(value, dict):
        for v in value.values():
            _read_only(v)
    elif isinstance(value, (tuple, list)):
        for v in value:
            _read_only(v)
    return value


class SourceModels:
    """Source model data that is read and prepared once and reused for many sites.

    Holds the contents of the source model files (npy files and all members of npz files) and
    their site-independent preparation (the rupture boundaries of fault models, the node of each
    rupture of point source models, and the geometry of the point source ruptures). The arrays are read-only, so
    they cannot be changed by accident and, after load_source_models in the parent process,
    are shared by worker processes created with fork. Pass the object to get_hazard (or
    get_source_data) with models=...; data that has not been loaded yet is loaded on first use.
    """

    def __init__(self):
        self._cache = {}

    def cached(self, key, make):
        """The value stored under key, made with make() (and made read-only) on first use."""
        if key not in self._cache:
            self._cache[key] = _read_only(make())
        return self._cache[key]

    def load(self, path, filename):
        """Contents of the file filename (npy: array; npz: dict-like with .files) in path."""

        def read():
            name = str(path.joinpath(filename))
            if filename.endswith(".npz"):
                with np.load(name) as data:
                    return _NpzArrays((k, data[k]) for k in data.files)
            return np.load(name)

        return self.cached(("file", str(path), filename), read)


def load_source_models(config=None):
    """Reads and prepares the source models of a config once, for many get_hazard calls.

    Args:
        config (dict or string): a config (or config file) as for get_hazard. The source
            models with a weight > 0 are read and prepared by computing the PSHA of the config
            site with them (without liquefaction, output, or disaggregation), so that the
            files, the site-independent preparation, and the ground motion models are loaded.
            None returns an empty SourceModels that loads the models on first use.

    Returns:
        models (SourceModels): pass it to get_hazard(config, models=models). The results are
        identical to those of get_hazard without models. Call this before creating worker
        processes with fork, so that the workers share the data.
    """
    models = SourceModels()
    if config is None:
        return models
    config = _read_config(config)
    warmup = {
        key: value
        for key, value in config.items()
        if key not in ("liquefaction_models", "profiles", "output")
    }
    warmup["output"] = {"psha": {"pga": [0.1]}}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        get_hazard(warmup, models=models)
    return models


_SCHEMA = {}


def _schema_validator():
    """Validator of ucla_plha_schema.json (lower case, as the config), made once per schema file."""
    path = str(files("ucla_plha").joinpath("ucla_plha_schema.json"))
    if path not in _SCHEMA:
        with open(path) as f:
            schema = json.loads(f.read().lower())
        cls = jsonschema.validators.validator_for(schema)
        cls.check_schema(schema)
        _SCHEMA[path] = cls(schema)
    return _SCHEMA[path]


def _read_config(config):
    """Config dict, in lower case, from a config file name or a config dict."""
    if isinstance(config, dict):
        return json.loads(json.dumps(config).lower())
    with open(config) as f:
        return json.loads(f.read().lower())


# Distances used by the ground motion models
GMM_DISTANCES = {
    "bssa14": {"rjb"},
    "ask14": {"rrup", "rx"},
    "cb14": {"rjb", "rrup", "rx"},
    "cy14": {"rjb", "rrup", "rx"},
}


def _distance_needs(gmms):
    """Set of distances ("rjb", "rrup", "rx") used by the ground motion models."""
    needs = set()
    for gmm in gmms:
        needs |= GMM_DISTANCES.get(gmm, set())
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


_LOCAL_GEOMETRY = {}


class _LocalGeometry:
    """Site-independent geometry of a fault source model for the local projection and GC2,
    loaded on first use and kept for the session: "tri_rjb" and "tri_rrup" (triangle vertices
    as latitude, longitude, depth), and, from the rectangles, "top_lat" and "top_lon" (upper
    edge of each piece, oriented so that the fault dips to its right), "piece_first", and
    "piece_count" (first piece and number of pieces of each segment)."""

    def __init__(self, path):
        self.path = path
        self.data = {}

    def __getitem__(self, name):
        if name not in self.data:
            if name in ("tri_rjb", "tri_rrup"):
                xyz = np.load(str(self.path.joinpath(f"{name}.npy")))
                lat, lon, depth = geometry.xyz_to_lat_lon_depth(xyz)
                self.data[name] = np.stack([lat, lon, np.maximum(depth, 0.0)], axis=-1)
            else:
                self.data.update(_piece_geometry(self.path))
        return self.data[name]


def _local_geometry(path):
    key = str(path)
    if key not in _LOCAL_GEOMETRY:
        _LOCAL_GEOMETRY[key] = _LocalGeometry(path)
    return _LOCAL_GEOMETRY[key]


def _piece_geometry(path):
    """Upper edges of the rectangles (pieces) of a fault source model and the pieces of each
    segment, with each segment oriented so that the fault dips to the right of its trace: a
    segment whose bottom edge is on the left (sum over its pieces of the flat-earth cross
    product of the upper edge, corner 0 -> 1, and the down-dip edge, corner 0 -> 2) is
    reversed as a whole (its pieces in reverse order, each from corner 1 to corner 0).
    Vertical segments keep their order."""
    rect = np.load(str(path.joinpath("rect_rjb.npy")))
    rid = np.load(str(path.joinpath("rect_segment_id.npy")))
    lat, lon, _ = geometry.xyz_to_lat_lon_depth(rect)
    c = np.cos(np.radians(lat[:, 0]))
    sx, sy = (lon[:, 1] - lon[:, 0]) * c, lat[:, 1] - lat[:, 0]
    bx, by = (lon[:, 2] - lon[:, 0]) * c, lat[:, 2] - lat[:, 0]
    n_seg = int(rid.max()) + 1
    first = np.searchsorted(rid, np.arange(n_seg))
    count = np.diff(np.r_[first, len(rid)])
    left = np.add.reduceat(sx * by - sy * bx, first) > 0.0
    # permutation of the pieces: reversed within the left-dipping segments
    seg = np.repeat(np.arange(n_seg), count)
    k = np.arange(len(rid)) - first[seg]
    order = np.where(left[seg], first[seg] + count[seg] - 1 - k, np.arange(len(rid)))
    lat, lon, rev = lat[order], lon[order], left[seg]
    return {
        "top_lat": np.where(rev[:, None], lat[:, [1, 0]], lat[:, [0, 1]]),
        "top_lon": np.where(rev[:, None], lon[:, [1, 0]], lon[:, [0, 1]]),
        "piece_first": first.astype(np.int64),
        "piece_count": count.astype(np.int64),
    }


def _gc2_setup(local, segment_index, boundaries):
    """Site-independent GC2 set-up of all ruptures of a model (geometry.gc2_rupture_setup),
    kept with the model's local geometry, and read from (or written to) the on-disk cache of
    ucla_plha.gc2_cache."""
    if "gc2_setup" not in local.data:
        model = local.path.name
        model_key = gc2_cache.key(local.path)
        expected = {
            "seg_length": len(local["piece_first"]),
            "piece_before": len(local["top_lat"]),
            "flip": len(segment_index),
            "offset": len(segment_index),
            "u_max": len(boundaries),
        }
        setup = gc2_cache.load(model, model_key, expected)
        if setup is None:
            setup = _compute_gc2_setup(local, segment_index, boundaries)
            gc2_cache.save(model, model_key, setup)
        local.data["gc2_setup"] = setup
    return local.data["gc2_setup"]


def _compute_gc2_setup(local, segment_index, boundaries):
    """The GC2 set-up of _gc2_setup, computed."""
    top_lat, top_lon = local["top_lat"], local["top_lon"]
    first = local["piece_first"]
    last = first + local["piece_count"] - 1
    piece_length = geometry.haversine_distance(
        top_lat[:, 0], top_lon[:, 0], top_lat[:, 1], top_lon[:, 1]
    )
    seg_length = np.add.reduceat(piece_length, first)
    csum = np.cumsum(piece_length)
    seg_of_piece = np.repeat(np.arange(len(first)), local["piece_count"])
    piece_before = csum - piece_length - (csum[first] - piece_length[first])[seg_of_piece]
    seg_lat = np.c_[top_lat[first, 0], top_lat[last, 1]]
    seg_lon = np.c_[top_lon[first, 0], top_lon[last, 1]]
    flip, offset, u_max = geometry.gc2_rupture_setup(
        seg_lat, seg_lon, seg_length, segment_index, boundaries
    )
    return {
        "seg_length": seg_length, "piece_before": piece_before,
        "flip": flip, "offset": offset, "u_max": u_max,
    }


def _gc2_rx_ry0(local, site_lat, site_lon, rupture, boundaries, segment_index):
    """GC2 Rx (= T) and Ry0 of the ruptures `rupture` (indices into the rupture boundaries),
    with each rupture as one GC2 coordinate system."""
    setup = _gc2_setup(local, segment_index, boundaries)
    x1, y1 = geometry.local_projection(local["top_lat"][:, 0], local["top_lon"][:, 0],
                                       site_lat, site_lon)
    x2, y2 = geometry.local_projection(local["top_lat"][:, 1], local["top_lon"][:, 1],
                                       site_lat, site_lon)
    W, A, B, on_seg, on_u = geometry.gc2_segment_sums(
        x1, y1, x2, y2, local["piece_first"], setup["piece_before"]
    )
    T, U = geometry.gc2_rupture_coordinates(
        W, A, B, on_seg, on_u, setup["seg_length"], segment_index, boundaries,
        setup["flip"], setup["offset"],
    )
    return T[rupture], geometry.gc2_ry0(U[rupture], setup["u_max"][rupture])


def get_source_data(
    source_type,
    source_model,
    p_xyz,
    dist_cutoff,
    m_min,
    gmms,
    extras=False,
    rupture_rx="gc2",
    triangle_frame="local",
    extra_distances=(),
    models=None,
):
    """Returns magnitude, fault type, rate, distance, and fault geometry terms.

    Args:
        source_type (string): Either "fault_source_models" or "point_source_models"
        source_model (string): Directory for source_model within the source_type directory.
            "ucerf3_fm31" or "ucerf3_fm32" for fault_source_models, and
            "ucerf3_fm31_grid_sub_seis", "ucerf3_fm31_grid_unassociated", "ucerf3_fm32_grid_sub_seis",
            or "ucerf3_fm32_grid_unassociated" for point_source_models
        p_xyz (numpy array, dtype=float): Array containing x, y, z coordinates for point of interest, length = 3
        dist_cutoff (float): maximum distance to consider in seismic hazard analysis
        m_min (float): minimum magnitude to consider in seismic hazard analysis
        gmms (array, dtype=string): An array of strings defining ground motion models to use in seismic
            hazard analysis: "ask14", "bssa14", "cb14", or "cy14". They determine which distances are
            computed and whether the distance cutoff is applied to rjb (if any model uses rjb) or rrup.
        extras (bool): if True, also return a dict of additional rupture data (see Returns)
        rupture_rx (str): how rx, rx1, and ry0 of fault ruptures are found. "gc2" (default):
            rx = T and ry0 from U of the GC2 coordinates of Spudich and Chiou (2015), with the
            upper edges of all pieces of the rupture as its trace (geometry.gc2_coordinates;
            jumps of more than 1 km between consecutive pieces start a new trace), as in the
            NGA-West2 and NGA-West3 flatfiles; rx1 = rx - the horizontal width of the segment
            with the smallest rrup. "closest_section": rx and rx1 of the segment with the
            smallest rrup (ties go to the lowest segment index), as nshmp-lib does for fault
            system ruptures (SystemRuptureSet.InputGenerator), and ry0 the minimum over the
            segments. "minimum": the minimum rx, rx1, and ry0 over the segments of the rupture
            (ucla_plha 1.x and earlier); because rx is signed, this is the most footwall-side
            segment of the rupture, which may be far from the site. rjb and rrup are always the
            minima over the segments.
        triangle_frame (str): coordinates of the rjb and rrup of the fault triangles. "local"
            (default): the triangles are projected on a sphere centred on the site
            (geometry.local_triangle_distance), so the distances are measured along the Earth's
            surface and rrup >= rjb. "geocentric": straight-line distances in Earth-centred
            Cartesian coordinates (ucla_plha 2.x and earlier), which are chords (about 0.1 km
            short at 400 km) and, because the vertical at the fault tilts toward a distant site,
            make rrup smaller than rjb beyond about 100 km.
        extra_distances (iterable): distances to compute in addition to those of the ground
            motion models, e.g. {"rrup"} for a disaggregation by rrup. They do not change the
            distance the cutoff is applied to.
        models (SourceModels): source model data read and prepared once
            (load_source_models). None (default): the files are read for this call.

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
        "index" (index of each returned rupture in ruptures.npz), of length N.

    Notes:
        N = number of events
        Point source distances use the Kaklamanos et al. (2011) and Thompson and Worden (2017)
        approximations.
    """
    if source_type not in ("fault_source_models", "point_source_models"):
        return None
    path = _source_model_path(source_type, source_model)
    if not os.path.isdir(str(path)):
        raise ValueError(f'source model "{source_model}" not found in {source_type}')
    if models is None:
        models = SourceModels()
    needs = _distance_needs(gmms)
    need_rjb = "rjb" in needs
    need_rrup = bool((needs | set(extra_distances)) & {"rrup", "rx", "ry0"})

    if source_type == "fault_source_models":
        # Read files required by all ground motion models
        # Read decompressed version of files if they exist. Otherwise read zipped version.
        if os.path.exists(str(path.joinpath("ruptures.npy.npz"))):
            ruptures = models.load(path, "ruptures.npy.npz")
        else:
            ruptures = models.load(path, "ruptures.npz")
        m = ruptures["m"]
        fault_type = ruptures["fault_type"]
        if os.path.exists(str(path.joinpath("ruptures_segments.npy.npz"))):
            ruptures_segments = models.load(path, "ruptures_segments.npy.npz")
        else:
            ruptures_segments = models.load(path, "ruptures_segments.npz")
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
        empty_array = np.empty(len(m))
        if rupture_rx not in ("gc2", "closest_section", "minimum"):
            raise ValueError(
                f'rupture_rx must be "gc2", "closest_section", or "minimum", not "{rupture_rx}"'
            )
        if triangle_frame not in ("local", "geocentric"):
            raise ValueError(
                f'triangle_frame must be "local" or "geocentric", not "{triangle_frame}"'
            )
        site_lat = np.degrees(np.arctan2(p_xyz[2], np.hypot(p_xyz[0], p_xyz[1])))
        site_lon = np.degrees(np.arctan2(p_xyz[1], p_xyz[0]))
        site_elevation = float(np.linalg.norm(p_xyz) - geometry.spheroid_radius(site_lat))
        local = None
        if triangle_frame == "local" or rupture_rx == "gc2":
            local = _local_geometry(path)
        if need_rjb or need_rrup:
            tri_segment_id = models.load(path, "tri_segment_id.npy")
        if need_rjb:
            if triangle_frame == "local":
                rjb_all = geometry.local_triangle_distance(
                    local["tri_rjb"], site_lat, site_lon, 0.0, tri_segment_id
                )
            else:
                # (copies: point_triangle_distance changes degenerate triangles in place)
                tri_rjb = models.load(path, "tri_rjb.npy").copy()
                rjb_all = geometry.point_triangle_distance(tri_rjb, p_xyz, tri_segment_id)
        if need_rrup:
            rect_segment_id = models.load(path, "rect_segment_id.npy")
            rect = models.load(path, "rect_rjb.npy")
            rx_all, rx1_all, ry0_all = geometry.get_Rx_Rx1_Ry0(
                rect, p_xyz, rect_segment_id
            )
            if triangle_frame == "local":
                rrup_all = geometry.local_triangle_distance(
                    local["tri_rrup"], site_lat, site_lon, site_elevation, tri_segment_id
                )
            else:
                tri_rrup = models.load(path, "tri_rrup.npy").copy()
                rrup_all = geometry.point_triangle_distance(
                    tri_rrup, p_xyz, tri_segment_id
                )
        boundaries = models.cached(
            ("boundaries", str(path)),
            lambda: np.r_[0, np.where(np.diff(ruptures_index) != 0)[0] + 1],
        )
        if need_rrup:
            rrup_seg = rrup_all[segment_index]
            rrup = np.minimum.reduceat(rrup_seg, boundaries)
            if rupture_rx in ("closest_section", "gc2"):
                # gc2: rx and ry0 are replaced below; rx - rx1 (the horizontal width of the
                # closest segment) is kept for rx1
                closest = _closest_segment(rrup_seg, rrup, segment_index, boundaries)
                rx = rx_all[closest]
                rx1 = rx1_all[closest]
            else:
                rx = np.minimum.reduceat(rx_all[segment_index], boundaries)
                rx1 = np.minimum.reduceat(rx1_all[segment_index], boundaries)
            ry0 = np.minimum.reduceat(ry0_all[segment_index], boundaries)
        else:
            rrup = empty_array
            rx = empty_array
            rx1 = empty_array
            ry0 = empty_array

        if need_rjb:
            rjb = np.minimum.reduceat(rjb_all[segment_index], boundaries)
        else:
            rjb = empty_array

        filter = _filter_mask(len(m), rjb if need_rjb else rrup, m, dist_cutoff, m_min)
        if need_rrup and rupture_rx == "gc2":
            # GC2 rx and ry0 of the ruptures that are kept
            keep = np.flatnonzero(filter)
            width = rx - rx1
            rx, rx1, ry0 = (np.array(a, dtype=np.float64) for a in (rx, rx1, ry0))
            if len(keep):
                rx_gc2, ry0_gc2 = _gc2_rx_ry0(
                    local, site_lat, site_lon, keep, boundaries, segment_index
                )
                rx[keep] = rx_gc2
                rx1[keep] = rx_gc2 - width[keep]
                ry0[keep] = ry0_gc2
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

    elif source_type == "point_source_models":
        ruptures = models.load(path, "ruptures.npz")
        rate = ruptures["rate"]
        m = ruptures["m"]
        fault_type = ruptures["style"]
        node_index = ruptures["node_index"]
        nodes_index = models.load(path, "node_index.npy")
        points = models.load(path, "points.npy")

        repi_all = np.sqrt(
            (points[:, 0] - p_xyz[0]) ** 2
            + (points[:, 1] - p_xyz[1]) ** 2
            + (points[:, 2] - p_xyz[2]) ** 2
        )
        node_of_rupture = models.cached(
            ("node_of_rupture", str(path)),
            lambda: np.searchsorted(
                nodes_index, node_index, sorter=np.argsort(nodes_index)
            ),
        )
        repi = repi_all[node_of_rupture]
        rjb = repi / (1.0 + np.exp(-1.05 * (np.log(repi) - 1.037 * m + 4.2776)))

        filter = _filter_mask(len(rjb), rjb, m, dist_cutoff, m_min)
        dip, w, zhyp, ztor, zbor = models.cached(
            ("crustal_point_source", str(path)),
            lambda: _crustal_point_source_geometry(m, fault_type),
        )
        # Expression for Rrup is fit to default parameters (first row in table 1) in Thompson and Worden (2017)
        d = np.radians(80)
        rrup1 = np.sqrt((repi - 0.6 * w * np.cos(d)) ** 2 + ztor**2)
        filt = repi < 0.6 * w * np.cos(d) - ztor * np.tan(d)
        rrup1[filt] = (zhyp[filt] - repi[filt] * np.tan(d)) / (
            np.cos(d) + np.sin(d) * np.tan(d)
        )
        rrup2 = zhyp.copy()
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


def _crustal_point_source_geometry(m, fault_type):
    """Dip, width, hypocentre depth, ztor, and zbor of the ucla_plha crustal point source
    ruptures (site-independent; see get_source_data)."""
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
    zbor = ztor + w * np.sin(dip * np.pi / 180.0)
    return dip, w, zhyp, ztor, zbor


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
        gmm (string): Ground motion model. One of "ask14", "bssa14", "cb14", "cy14"
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
    else:
        raise ValueError(
            f'incorrect ground motion model "{gmm}", '
            'expected one of "ask14", "bssa14", "cb14", "cy14"'
        )
    return [mu_ln_pga, sigma_ln_pga]


def get_liquefaction_cdfs(
    m, mu_ln_pga, sigma_ln_pga, fsl, liquefaction_model, config, need_eps=True,
    truncation_level=None,
):
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
        need_eps (bool): if False, the epsilons are not needed (no disaggregation), and models
            that compute them only for the disaggregation (ngl_smt_2024) return None for eps
        truncation_level (float): upper truncation of the ground motion distribution in standard
            deviations (config constraints.truncation_level), used by ngl_smt_2024 with
            "pga_uncertainty" "shared" (the other models are untruncated). None: untruncated.

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
        # Optional input parameter uncertainty (Cetin et al. 2018b); off by default
        parameter_uncertainty = c.get("parameter_uncertainty", False)
        uncertainties = {
            key: c.get(key, 0.0)
            for key in ["sigma_n160", "sigma_fc", "cov_sigmavp", "cov_csr"]
        }
        if parameter_uncertainty and not any(uncertainties.values()):
            warnings.warn(
                'cetin_et_al_2018: "parameter_uncertainty" is true, but "sigma_n160", "sigma_fc", '
                '"cov_sigmavp", and "cov_csr" are all zero, so only the model uncertainty is included',
                UserWarning,
                stacklevel=2,
            )
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
            parameter_uncertainty=parameter_uncertainty,
            **uncertainties,
        )
    elif liquefaction_model == "moss_et_al_2006":
        c = config["liquefaction_models"]["moss_et_al_2006"]
        _reject_removed_options(config["liquefaction_models"])
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
            pga_uncertainty=_pga_uncertainty(liquefaction_model, c),
            truncation_level=truncation_level,
            rd_amax_cap=c.get("rd_amax_cap", moss_et_al_2006.RD_AMAX_CAP),
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
            pga_uncertainty=_pga_uncertainty(liquefaction_model, c),
            truncation_level=truncation_level,
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
            return_eps=need_eps,
            pga_uncertainty=_pga_uncertainty(liquefaction_model, c),
            truncation_level=truncation_level,
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

        Bins include their lower edge and exclude their upper edge. Events outside the bins
        (below the first or at or above the last edge of any of m, r, or eps) are not binned,
        so the bins then sum to less than the hazard (the "residual" of the nshmp-lib
        disaggregation). Use wide outer edges (e.g. -inf and inf for epsilon) to bin all events.
    """
    # use Numpy digitize function to assign bin numbers
    m_hazard = np.digitize(m, m_bin_edges)
    r_hazard = np.digitize(r, r_bin_edges)
    eps_hazard = np.digitize(eps, eps_bin_edges)
    in_bins = disagg_in_bins(m, r, eps, m_bin_edges, r_bin_edges, eps_bin_edges)

    # compute number of bins per intensity measure value
    Nbins = (len(m_bin_edges) - 1) * (len(r_bin_edges) - 1) * (len(eps_bin_edges) - 1)

    # use tensor rank reduction to create L x Nbins length array of indices
    bin_indices = (
        eps_hazard
        - 1
        + (r_hazard - 1) * (len(eps_bin_edges) - 1)
        + (m_hazard - 1) * (len(eps_bin_edges) - 1) * (len(r_bin_edges) - 1)
    )
    # events outside the bins go to bin 0 with zero weight
    bin_indices = np.where(in_bins, bin_indices, 0)
    hazards = np.where(in_bins, hazards, 0.0)

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


def disagg_in_bins(m, r, eps, m_bin_edges, r_bin_edges, eps_bin_edges):
    """Whether each event falls inside the disaggregation bins (see get_disagg), shape = L x N."""

    def inside(x, edges):
        return (x >= edges[0]) & (x < edges[-1])

    return (
        inside(np.asarray(eps), eps_bin_edges)
        & inside(np.asarray(m), m_bin_edges)
        & inside(np.asarray(r), r_bin_edges)
    )


def disagg_distance(distance_metric, gmms, rjb, rrup):
    """Distance used to bin a source model's ruptures in the disaggregation.

    distance_metric: "default" (or None): rjb if any of the ground motion models gmms uses rjb
    (or there are no models), otherwise rrup, as in ucla_plha 1.x; "rrup": rrup, as the USGS
    NSHM disaggregation (nshmp-lib Disaggregator bins and averages HazardInput.rRup for every
    source type, and its web service labels the distance "Closest Distance, rRup"); "rjb": rjb.
    """
    if distance_metric in (None, "default"):
        return rjb if "rjb" in _distance_needs(gmms) or not len(gmms) else rrup
    if distance_metric == "rrup":
        return rrup
    if distance_metric == "rjb":
        if len(rjb) != len(rrup):
            raise ValueError(
                'disaggregation distance_metric "rjb" needs a ground motion model that uses rjb'
            )
        return rjb
    raise ValueError(f'unknown disaggregation distance_metric "{distance_metric}"')


def _disagg_sums(hazards, m, r, eps, in_bins):
    """Binned hazard and its m, r, and eps moments, shape = 4 x L (see get_hazard "means")."""
    h = np.where(in_bins, hazards, 0.0)
    return np.stack([h.sum(axis=1), h @ m, h @ r, np.sum(h * eps, axis=1)])


def _disagg_means(sums, total):
    """Contribution (% of total) and means of binned hazard sums (_disagg_sums)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        binned = sums[0]
        out = {
            "contribution": (100.0 * binned / total).tolist(),
            "mean_magnitude": (sums[1] / binned).tolist(),
            "mean_distance": (sums[2] / binned).tolist(),
            "mean_epsilon": (sums[3] / binned).tolist(),
        }
    return {k: [None if not np.isfinite(x) else x for x in v] for k, v in out.items()}


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


def _normalize_weights(models):
    """Normalizes the "weight" of the entries of a {name: {"weight": w}} dict in place so that
    they sum to one (if their sum is positive)."""
    total = sum(entry.get("weight", 0.0) for entry in models.values())
    if total > 0:
        for entry in models.values():
            entry["weight"] /= total


# Liquefaction model options that were removed, with the reason; a config that gives one is
# rejected so that it does not silently change meaning
REMOVED_LIQUEFACTION_OPTIONS = {
    ("moss_et_al_2006", "duration_weighting"): (
        "the duration weighting factor (CSR / DWF) of Moss et al. (2006) was removed: their "
        "Eq. 20 and 21 have their own magnitude term and are compared with the CSR of the event "
        "(Eq. 21 note and Table 1), so the DWF double counted the magnitude effect. The 2.x Moss "
        "et al. (2006) results, which divided CSR by the DWF, can't be reproduced with this "
        "version; remove \"duration_weighting\" from the config"
    ),
}


def _reject_removed_options(liquefaction_models):
    """Raises ValueError if a liquefaction model entry gives a removed option."""
    for (model, key), reason in REMOVED_LIQUEFACTION_OPTIONS.items():
        if key in (liquefaction_models.get(model) or {}):
            raise ValueError(f'{model}: "{key}" is not an option: {reason}')


LIQUEFACTION_MODELS = [
    "boulanger_idriss_2012",
    "boulanger_idriss_2016",
    "cetin_et_al_2018",
    "moss_et_al_2006",
    "ngl_smt_2024",
]


def _normalize_liquefaction_weights(liquefaction_models):
    """Divides the weights of the liquefaction models (config "liquefaction_models") by their sum."""
    liquefaction_model_weight_sum = 0.0
    for liquefaction_model in LIQUEFACTION_MODELS:
        liquefaction_model_weight_sum += liquefaction_models.get(liquefaction_model, {}).get(
            "weight", 0.0
        )
    if liquefaction_model_weight_sum > 0:
        for liquefaction_model in LIQUEFACTION_MODELS:
            if liquefaction_models.get(liquefaction_model, {}):
                liquefaction_models[liquefaction_model]["weight"] /= liquefaction_model_weight_sum


# Liquefaction models with a treatment of the peak acceleration conditioned on PGA: the value of
# their "pga_uncertainty" for that treatment, and the default value (change a default here, e.g.
# to "closed_form", to make the earlier treatment the default of a model)
CONDITIONAL_PGA_UNCERTAINTY = {
    "boulanger_idriss_2016": "conditional",
    "moss_et_al_2006": "conditional",
    "ngl_smt_2024": "shared",
}
DEFAULT_PGA_UNCERTAINTY = {
    "boulanger_idriss_2016": "conditional",
    "moss_et_al_2006": "conditional",
    "ngl_smt_2024": "shared",
}
_NGL_LAYER_KEYS = ("ztop", "zbot", "qc1ncs", "ic", "sigmav", "sigmavp", "ksat")


def _pga_uncertainty(liquefaction_model, c):
    """The "pga_uncertainty" of a liquefaction model entry c (DEFAULT_PGA_UNCERTAINTY if not given)."""
    return c.get("pga_uncertainty", DEFAULT_PGA_UNCERTAINTY[liquefaction_model])


def _conditional_cdfs(liquefaction_model, c, fsl):
    """conditional_cdfs(m, ln_pga) of a liquefaction model entry c: P(FS < fsl | m, PGA), K x L."""
    pa = float(c.get("pa", 101.325))
    if liquefaction_model == "ngl_smt_2024":
        layers = [np.asarray(c[k], dtype=float) for k in _NGL_LAYER_KEYS]

        def conditional_cdfs(m, ln_pga):
            return ngl_smt_2024.get_conditional_fsl_cdfs(m, ln_pga, *layers, fsl, pa)

    elif liquefaction_model == "boulanger_idriss_2016":

        def conditional_cdfs(m, ln_pga):
            return boulanger_idriss_2016.get_conditional_fsl_cdfs(
                m, ln_pga, c["sigmav"], c["sigmavp"], c["depth"], c["qc1ncs"], fsl, pa
            )

    else:

        def conditional_cdfs(m, ln_pga):
            return moss_et_al_2006.get_conditional_fsl_cdfs(
                m, ln_pga, c["sigmav"], c["sigmavp"], c["depth"], c["qc"], c["fs"], fsl, pa,
                c.get("rd_amax_cap", moss_et_al_2006.RD_AMAX_CAP),
            )

    return conditional_cdfs


def _liquefaction_tables(profiles, fsl, truncation_level):
    """Rate grids and conditional probability tables of the liquefaction models conditioned on
    PGA (CONDITIONAL_PGA_UNCERTAINTY) and integrated on the magnitude x PGA grid
    ("pga_integration" "grid", the default).

    Sets profile["tables"] ({liquefaction model: pga_grid.ConditionalTable} of the models of
    the profile that use the grid) and returns the rate grids (pga_grid.RateGrid), one for each
    set of grid parameters ("pga_grid"), so that all profiles and liquefaction models with the
    same parameters share one grid.
    """
    grids = {}
    for profile in profiles:
        profile["tables"] = {}
        for liquefaction_model, conditional in CONDITIONAL_PGA_UNCERTAINTY.items():
            c = profile["models"].get(liquefaction_model)
            if (
                not c
                or c.get("weight", 0.0) <= 0
                or _pga_uncertainty(liquefaction_model, c) != conditional
                or c.get("pga_integration", "grid") != "grid"
            ):
                continue
            params = pga_grid.grid_parameters(c.get("pga_grid"))
            key = tuple(sorted(params.items()))
            if key not in grids:
                grids[key] = pga_grid.RateGrid(truncation_level, **params)
            profile["tables"][liquefaction_model] = pga_grid.ConditionalTable(
                grids[key], _conditional_cdfs(liquefaction_model, c, fsl)
            )
    return grids


def get_hazard(config_file, models=None):
    """Reads config file and runs PSHA and PLHA

    Inputs:
        config_file (string or dict): Filename, including path, of config file, or a config
            dict (as read from a config file; it is not changed). Must follow the schema
            defined by ucla_plha_schema.json. See documentation for more thorough documentation
            of the config file.
        models (SourceModels): source models read and prepared once (load_source_models), to
            run many sites without reading the source model files again. None (default): the
            source models are read for this call. The results are the same.

    Returns:
        output (dict): Python dictionary containing output of analysis. The output contains all
            of the inputs for preservation, along with the hazard curve(s) and any requested
            disaggregation data. See documentation for more thorough description of output.

    Notes:
        The weights of the fault source models, of the point source models, of the ground
        motion models, and of the liquefaction models (of each soil profile) are each
        normalized to sum to one. The hazard is the sum of the weighted fault source hazard
        and the weighted point source hazard. With "profiles", the PLHA of each soil profile
        is in output "profiles" (a list of {"name", "plha"}) instead of "plha".
    """
    # Validate config_file against schema. If ngl_smt_2024 liquefaction model is used, the cpt_data file is
    # validated in the get_liquefaction_hazards function.

    validator = _schema_validator()
    config = _read_config(config_file)

    # validate config file and return messageg if errors are encountered
    try:
        error = jsonschema.exceptions.best_match(validator.iter_errors(config))
        if error is not None:
            raise error
    except jsonschema.ValidationError as e:
        print("Config File Error:", e.message)
        return
    if "profiles" in config and "liquefaction_models" in config:
        print(
            'Config File Error: a config with "profiles" gives the liquefaction models of each '
            'profile, not "liquefaction_models"'
        )
        return
    if isinstance(config_file, dict) and config["output"].get("outputfile") == "default":
        raise ValueError(
            'outputfile "default" needs a config file name; give the output file name for '
            "a config dict"
        )
    if models is None:
        models = SourceModels()

    # normalize weights in config file
    for type_models in config["source_models"].values():
        _normalize_weights(type_models)
    _normalize_weights(config["ground_motion_models"])

    # Soil profiles: the liquefaction models of the config, or of each of its "profiles" (the
    # source data, ground motions, and PSHA are computed once for all profiles)
    if "profiles" in config:
        profiles = [
            {
                "name": profile.get("name", str(i)),
                "models": profile["liquefaction_models"],
                "config": profile,
            }
            for i, profile in enumerate(config["profiles"])
        ]
    else:
        profiles = [
            {"name": None, "models": config.get("liquefaction_models", {}), "config": config}
        ]
    for profile in profiles:
        _reject_removed_options(profile["models"])
        _normalize_liquefaction_weights(profile["models"])

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
    dist_cutoff = config.get("constraints", {}).get("dist_cutoff", None)
    m_min = config.get("constraints", {}).get("m_min", None)
    truncation_level = config.get("constraints", {}).get("truncation_level", None)
    rupture_rx = config.get("constraints", {}).get("rupture_rx", "gc2")
    triangle_frame = config.get("constraints", {}).get("triangle_frame", "local")

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
            psha_distance_metric = config["output"]["psha"]["disaggregation"].get(
                "distance_metric", "default"
            )
            psha_disagg_means = config["output"]["psha"]["disaggregation"].get(
                "means", False
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
            psha_distance_metric = "default"
            psha_disagg_means = False
        seismic_hazard = np.zeros(len(pga))
        # binned hazard and its m, r, eps moments (total and per source model)
        disagg_sums = np.zeros((4, len(pga)))
        source_disagg_sums = {}
    else:
        output_psha = False
        output_psha_disaggregation = False
        output_source_hazard = False
        psha_distance_metric = "default"
        psha_disagg_means = False

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
            plha_distance_metric = config["output"]["plha"]["disaggregation"].get(
                "distance_metric", "default"
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
            for profile in profiles:
                profile["disagg"] = np.zeros(
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
            plha_distance_metric = "default"
        for profile in profiles:
            profile["hazard"] = np.zeros(len(fsl))
    else:
        output_plha = False
        output_plha_disaggregation = False
        plha_distance_metric = "default"

    # Liquefaction models conditioned on PGA (NGL SMT 2024 with the PGA shared by the layers of the
    # profile, Boulanger and Idriss 2016, Moss et al. 2006; the defaults): the rates of the ground
    # motions of the events are accumulated on a magnitude x ln(PGA) grid (pga_grid.RateGrid),
    # built once for all profiles and liquefaction models, and the hazard of each profile and
    # model is the sum over the grid of the rates times its conditional probabilities (added
    # after the loops). The disaggregation uses per-event probabilities from the grid nodes.
    liquefaction_grids = (
        _liquefaction_tables(profiles, fsl, truncation_level) if output_plha else {}
    )
    on_grid = bool(liquefaction_grids) and not output_plha_disaggregation

    # all ground motion models determine the distance types, and rrup is also computed if the
    # disaggregation is binned by rrup
    gmms = list(config["ground_motion_models"].keys())
    extra_distances = (
        {"rrup"} if "rrup" in (psha_distance_metric, plha_distance_metric) else set()
    )

    # Loop over source models. We have fault_source_models and point_source_models, so there are two loops
    for source_model in config["source_models"].keys():
        for fault_source_model in config["source_models"][source_model].keys():
            source_model_weight = config["source_models"][source_model][
                fault_source_model
            ]["weight"]
            # move on to next source model if weight is less than or equal to zero
            if source_model_weight <= 0:
                continue
            m, fault_type, rate, rjb, rrup, rx, rx1, ry0, dip, ztor, zbor = get_source_data(
                source_model,
                fault_source_model,
                p_xyz,
                dist_cutoff,
                m_min,
                gmms,
                rupture_rx=rupture_rx,
                triangle_frame=triangle_frame,
                extra_distances=extra_distances,
                models=models,
            )
            # distances used for disaggregation (see disagg_distance)
            if output_psha_disaggregation:
                r_disagg = disagg_distance(psha_distance_metric, gmms, rjb, rrup)
            if output_plha_disaggregation:
                r_disagg_plha = disagg_distance(plha_distance_metric, gmms, rjb, rrup)

            # Loop over ground motion models.
            for gmm in gmms:
                ground_motion_model_weight = config["ground_motion_models"][gmm]["weight"]
                # move on to next ground motion model if weight is less than or equal to zero
                if ground_motion_model_weight <= 0:
                    continue
                mu_ln_pga, sigma_ln_pga = get_ground_motion_data(
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
                )
                # Compute seismic hazard if requested in config file
                if output_psha:
                    eps, p = get_exceedance(pga, mu_ln_pga, sigma_ln_pga, truncation_level)
                    seismic_hazards = p * rate
                    curve = np.sum(seismic_hazards, axis=1)
                    weight = source_model_weight * ground_motion_model_weight
                    seismic_hazard += weight * curve
                    if output_source_hazard:
                        source_hazard[fault_source_model] = (
                            source_hazard.get(fault_source_model, 0.0) + weight * curve
                        )
                    # Compute seismic hazard disaggregation if requested in config file
                    if output_psha_disaggregation:
                        psha_disagg += weight * get_disagg(
                            seismic_hazards,
                            m,
                            r_disagg,
                            eps,
                            psha_magnitude_bin_edges,
                            psha_distance_bin_edges,
                            psha_epsilon_bin_edges,
                        )
                        if psha_disagg_means:
                            in_bins = disagg_in_bins(
                                m,
                                r_disagg,
                                eps,
                                psha_magnitude_bin_edges,
                                psha_distance_bin_edges,
                                psha_epsilon_bin_edges,
                            )
                            sums = _disagg_sums(
                                weight * seismic_hazards, m, r_disagg, eps, in_bins
                            )
                            disagg_sums += sums
                            source_disagg_sums[fault_source_model] = (
                                source_disagg_sums.get(fault_source_model, 0.0) + sums
                            )
                if liquefaction_grids:
                    sigma_events = np.broadcast_to(sigma_ln_pga, np.shape(mu_ln_pga))
                if on_grid:
                    for grid in liquefaction_grids.values():
                        grid.add(
                            m,
                            mu_ln_pga,
                            sigma_events,
                            source_model_weight * ground_motion_model_weight * rate,
                        )
                # per-event probabilities of the models on the grid (for the disaggregation), with
                # the quadrature weights of the events shared by the profiles and models of a grid
                grid_event_cdfs = {}
                if output_plha and not on_grid:
                    grid_tables = [
                        (i, model, table)
                        for i, profile in enumerate(profiles)
                        for model, table in profile["tables"].items()
                    ]
                    if grid_tables:
                        cdfs = pga_grid.shared_event_cdfs(
                            [table for _, _, table in grid_tables], m, mu_ln_pga, sigma_events
                        )
                        grid_event_cdfs = {
                            (i, model): p for (i, model, _), p in zip(grid_tables, cdfs)
                        }
                # Compute liquefaction hazard if requested in config file, for each soil profile
                for i, profile in enumerate(profiles):
                    for liquefaction_model in profile["models"].keys():
                        liquefaction_model_weight = profile["models"][liquefaction_model][
                            "weight"
                        ]
                        if liquefaction_model_weight <= 0:
                            continue
                        if output_plha:
                            table = profile["tables"].get(liquefaction_model)
                            if table is not None and on_grid:
                                # the events are on the rate grid
                                continue
                            if table is not None:
                                liquefaction_hazards = grid_event_cdfs[(i, liquefaction_model)]
                                eps = (
                                    ngl_smt_2024.get_eps(liquefaction_hazards)
                                    if liquefaction_model == "ngl_smt_2024"
                                    else ndtri(liquefaction_hazards)
                                )
                            else:
                                liquefaction_hazards, eps = get_liquefaction_cdfs(
                                    m,
                                    mu_ln_pga,
                                    sigma_ln_pga,
                                    fsl,
                                    liquefaction_model,
                                    profile["config"],
                                    need_eps=output_plha_disaggregation,
                                    truncation_level=truncation_level,
                                )
                            liquefaction_hazards *= rate[:, np.newaxis]
                            profile["hazard"] += (
                                source_model_weight
                                * ground_motion_model_weight
                                * liquefaction_model_weight
                                * np.sum(liquefaction_hazards, axis=0)
                            )
                            # Compute liquefaction hazard disaggregation if requested in config file
                            if output_plha_disaggregation:
                                eps = eps.T
                                liquefaction_hazards = liquefaction_hazards.T
                                profile["disagg"] += (
                                    source_model_weight
                                    * ground_motion_model_weight
                                    * liquefaction_model_weight
                                    * get_disagg(
                                        liquefaction_hazards,
                                        m,
                                        r_disagg_plha,
                                        eps,
                                        plha_magnitude_bin_edges,
                                        plha_distance_bin_edges,
                                        plha_epsilon_bin_edges,
                                    )
                                )
    # Liquefaction hazard from the rate grid
    if on_grid:
        for profile in profiles:
            for liquefaction_model, table in profile["tables"].items():
                profile["hazard"] += (
                    profile["models"][liquefaction_model]["weight"] * table.curve()
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
            if psha_disagg_means:
                output["output"]["psha"]["disaggregation_means"] = _disagg_means(
                    disagg_sums, seismic_hazard
                )
        else:
            output["output"]["psha"] = {
                "PGA": pga.tolist(),
                "annual_rate_of_exceedance": seismic_hazard.tolist(),
            }
        if output_source_hazard:
            output["output"]["psha"]["source_model_hazard"] = {}
            for key, curve in source_hazard.items():
                output["output"]["psha"]["source_model_hazard"][key] = {
                    "annual_rate_of_exceedance": np.asarray(curve).tolist(),
                }
                if psha_disagg_means:
                    output["output"]["psha"]["source_model_hazard"][key][
                        "disaggregation_means"
                    ] = _disagg_means(
                        source_disagg_sums.get(key, np.zeros((4, len(pga)))),
                        seismic_hazard,
                    )
    if output_plha:
        for profile in profiles:
            liquefaction_hazard = profile["hazard"]
            if output_plha_disaggregation:
                plha_disagg = profile["disagg"]
                for i in range(len(fsl)):
                    plha_disagg[i] = plha_disagg[i] / liquefaction_hazard[i] * 100.0
                plha_output = {
                    "FSL": fsl.tolist(),
                    "annual_rate_of_nonexceedance": liquefaction_hazard.tolist(),
                    "disaggregation": plha_disagg.tolist(),
                }
            else:
                plha_output = {
                    "FSL": fsl.tolist(),
                    "annual_rate_of_nonexceedance": liquefaction_hazard.tolist(),
                }
            if profile["name"] is None:
                output["output"]["plha"] = plha_output
            else:
                # one entry per soil profile, in the order of the config "profiles"
                output["output"].setdefault("profiles", []).append(
                    {"name": profile["name"], "plha": plha_output}
                )

    if "outputfile" in config["output"].keys():
        if config["output"]["outputfile"] == "default":
            outputfilename = config_file.split(".json")[0] + "_output.json"
        else:
            outputfilename = config["output"]["outputfile"]
        with open(outputfilename, "w") as outputfile:
            json.dump(output, outputfile, indent=4)

    return output

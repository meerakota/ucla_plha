"""Tectonic regions, source model information, and ground motion model logic trees.

Every source model directory has a ``source_info.json`` file that gives its tectonic region,
the USGS NSHM hazard component it belongs to, and the settings needed to compute distances
(point sources) and cluster hazard (cluster models). Ground motion models are specified per
tectonic region in the config file, and each source model uses only the ground motion models
of its region.
"""

import json
import os
import re
import warnings

from ucla_plha import pygmm_gmms

TECTONIC_REGIONS = [
    "active_crust",
    "stable_crust",
    "subduction_interface",
    "subduction_slab",
]

NSHM_COMPONENTS = [
    "Fault",
    "FaultCluster",
    "FaultSystem",
    "Grid",
    "Interface",
    "Slab",
    "Zone",
]

#: Ground motion model logic trees of the current USGS NSHM (nshm-conus 6.2.0 gmm-tree.json
#: files; the stable crust weights 0.3333 and 0.1667 are 1/3 and 1/6). Models that are not in
#: pygmm (AM_09, ZHAO_06) are removed when the tree is used, and the remaining weights are
#: renormalized (with a warning). NGA-West2 ids that pygmm does not provide yet are replaced
#: by the ucla_plha models (with a warning).
DEFAULT_GMM_TREES = {
    "active_crust": {
        "ASK_14_BASIN": 0.25,
        "BSSA_14_BASIN": 0.25,
        "CB_14_BASIN": 0.25,
        "CY_14_BASIN": 0.25,
    },
    "stable_crust": {
        "NGA_EAST_2026": 1.0 / 3.0,
        "NGA_EAST_2026_ADJUSTED": 1.0 / 3.0,
        "NGA_EAST_SEEDS_2026": 1.0 / 6.0,
        "NGA_EAST_SEEDS_2026_ADJUSTED": 1.0 / 6.0,
    },
    "subduction_interface": {
        "AM_09_INTERFACE_BASIN": 0.125,
        "ZHAO_06_INTERFACE_BASIN": 0.125,
        "AG_20_CASCADIA_INTERFACE_ADJUSTED_BASIN": 0.25,
        "KBCG_20_CASCADIA_INTERFACE_BASIN": 0.25,
        "PSBAH_20_CASCADIA_INTERFACE_BASIN": 0.25,
    },
    "subduction_slab": {
        "ZHAO_06_SLAB_BASIN": 0.25,
        "AG_20_CASCADIA_SLAB_BASIN": 0.0825,
        "AG_20_CASCADIA_SLAB_ADJUSTED_BASIN": 0.1675,
        "KBCG_20_CASCADIA_SLAB_BASIN": 0.25,
        "PSBAH_20_CASCADIA_SLAB_BASIN": 0.25,
    },
}

#: USGS NSHM maximum distances (nshm-conus 6.2.0 gmm-config.json max-distance, km)
NSHM_MAX_DISTANCE = {
    "active_crust": 300.0,
    "stable_crust": 1000.0,
    "subduction_interface": 1000.0,
    "subduction_slab": 300.0,
}


#: nshmp-lib grid optimization settings (grid-config.json opt-distance-bin, smoothing-density,
#: smoothing-limit, grid-spacing) of NSHM grid models whose source_info.json describes them
#: only in text. nshmp-lib applies them by default (optimizeGrids, smoothGrids).
NSHMP_GRID_OPTIMIZATION = {
    "nshm23_ceus_grid": {
        "distance_bin": 5.0,
        "smoothing": {"density": 4, "limit": 40.0, "grid_spacing": 0.1},
    },
    "nshm23_ceus_grid_system": {
        "distance_bin": 5.0,
        "smoothing": {"density": 4, "limit": 40.0, "grid_spacing": 0.1},
    },
}

# ---------------------------------------------------------------------------------------
# Source model information
# ---------------------------------------------------------------------------------------


def default_source_info(source_type, name):
    """source_info of a source model directory without a source_info.json file.

    Source models without the file are treated as active crust models with the ucla_plha
    crustal point-source treatment (the behavior before source_info.json existed).
    """
    return {
        "name": name,
        "tectonic_region": "active_crust",
        "nshm_component": "FaultSystem" if source_type == "fault_source_models" else "Grid",
        "cluster": False,
    }


def _first(d, *keys, default=None):
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return default


def normalize_point_source(ps):
    """Normalize the point_source settings of a source_info.json.

    Accepts the ucla_plha names (point_source_type, rupture_scaling, max_depth, max_width or width_km,
    grid_depth_map, smoothing, distance_bin, depth) and the nshmp-lib config names
    (point-source-type, rupture-scaling, max-depth, max-width, grid-depth-map,
    smoothing-density, smoothing-limit, grid-spacing, opt-distance-bin).

    Returns a dict with keys:
        distance_method: "ucla_plha_crustal" (Kaklamanos et al. 2011 and Thompson and
            Worden 2017 approximations used for the UCERF3 and NSHM23 WUS grids) or "nshmp"
        type: "point" or "finite"
        rupture_scaling: lower-case nshmp-lib RuptureScaling name
        max_depth, max_width: finite rupture limits (km), one of them is None
        depth: where the depth to the top of rupture comes from: "rupture" (ruptures.npz
            "depth"), "node" (node depth, depth.npy in the model directory), "depth_map", or
            "auto" (the first of these that is available)
        grid_depth_map: list of {"m_min", "m_max", "depths", "weights"} or None
        smoothing: {"density", "limit", "grid_spacing"} or None
        distance_bin: distance bin width (km) of the nshmp-lib grid optimization, or None
    """
    ps = dict(ps or {})
    source_type = _first(ps, "point_source_type", "type", "point-source-type")
    method = _first(ps, "distance_method")
    if method is None:
        method = "nshmp" if source_type is not None else "ucla_plha_crustal"
    method = method.lower()
    out = {"distance_method": method}
    if method == "ucla_plha_crustal":
        return out
    if method != "nshmp":
        raise ValueError(f'unknown point source distance_method "{method}"')
    # The first word of descriptive values is used, e.g. "FINITE (GridSourceFinite)"
    out["type"] = re.split(r"[\s(]", str(source_type or "finite").strip())[0].lower()
    if out["type"] not in ("point", "finite", "fixed_strike"):
        raise ValueError(f'unsupported point source type "{out["type"]}"')
    out["rupture_scaling"] = _first(ps, "rupture_scaling", "rupture-scaling", default="none").lower()
    depth_info = ps.get("depth") if isinstance(ps.get("depth"), dict) else {}
    out["max_depth"] = _first(ps, "max_depth", "max-depth", "max_depth_km") or _first(
        depth_info, "max_depth_km", "max_depth", "max-depth"
    )
    out["max_width"] = _first(ps, "max_width", "max-width", "width_km")

    depth_map = _first(ps, "grid_depth_map", "grid-depth-map") or _first(
        depth_info, "depth_map", "grid_depth_map", "grid-depth-map"
    )
    if isinstance(depth_map, dict):
        # nshmp-lib format: {"label": {"mMin": .., "mMax": .., "depth-tree": [{"weight", "value"}]}}
        entries = []
        for entry in depth_map.values():
            tree = entry.get("depth-tree", entry.get("depth_tree"))
            if tree is None and "ztor_km_and_weight" in entry:
                tree = [{"value": d, "weight": w} for d, w in entry["ztor_km_and_weight"]]
            entries.append(
                {
                    "m_min": float(_first(entry, "mMin", "mmin", "m_min")),
                    "m_max": float(_first(entry, "mMax", "mmax", "m_max")),
                    "depths": [float(b["value"]) for b in tree],
                    "weights": [float(b["weight"]) for b in tree],
                }
            )
        depth_map = entries
    out["grid_depth_map"] = depth_map

    # "rupture", "node", "depth_map", or "auto" (ruptures.npz "depth" if present, else
    # depth.npy, else the depth map); descriptive text is treated as "auto"
    depth = str(_first(ps, "depth", default="auto")).lower()
    if depth not in ("rupture", "node", "depth_map"):
        depth = "auto"
    out["depth"] = depth

    smoothing = _first(ps, "smoothing")
    density = _first(ps, "smoothing_density", "smoothing-density")
    if smoothing is None and density is not None:
        smoothing = {
            "density": density,
            "limit": _first(ps, "smoothing_limit", "smoothing-limit", default=40.0),
            "grid_spacing": _first(ps, "grid_spacing", "grid-spacing", default=0.1),
        }
    if smoothing is not None and smoothing is not False:
        smoothing = {
            "density": int(_first(smoothing, "density", "smoothing-density")),
            "limit": float(_first(smoothing, "limit", "smoothing-limit", default=40.0)),
            "grid_spacing": float(_first(smoothing, "grid_spacing", "grid-spacing", default=0.1)),
        }
    else:
        smoothing = None
    distance_bin = _first(ps, "distance_bin", "opt_distance_bin", "opt-distance-bin")
    text = ps.get("optimization")
    if distance_bin is None and isinstance(text, str):
        # descriptive nshmp-lib grid optimization settings, e.g. "opt-distance-bin 5.0 km",
        # "smoothing-limit 40.0 km", "smoothing-density 4"
        found = re.search(r"opt-distance-bin\s+([0-9.]+)", text)
        if found:
            distance_bin = float(found.group(1))
        dens = re.search(r"smoothing-density\s+([0-9]+)", text)
        if smoothing is None and dens:
            limit = re.search(r"smoothing-limit\s+([0-9.]+)", text)
            spacing = re.search(r"grid spacing\s+([0-9.]+)", text)
            smoothing = {
                "density": int(dens.group(1)),
                "limit": float(limit.group(1)) if limit else 40.0,
                "grid_spacing": float(spacing.group(1)) if spacing else 0.1,
            }
    if out["type"] == "fixed_strike":
        # nshmp-lib only optimizes GRID rupture sets (not zones)
        distance_bin, smoothing = None, None
    out["smoothing"] = smoothing
    out["distance_bin"] = distance_bin
    return out


def read_source_info(path, source_type, name):
    """Read source_info.json in a source model directory (path), with defaults.

    The returned dict always has "name", "tectonic_region" (lower case), "nshm_component",
    "cluster" (bool), "cluster_nshm_component", "cluster_rupture_rate", and, for point
    source models, "point_source" (normalized by :func:`normalize_point_source`).
    """
    filename = os.path.join(str(path), "source_info.json")
    info = default_source_info(source_type, name)
    if os.path.exists(filename):
        with open(filename) as f:
            info.update(json.load(f))
    info["tectonic_region"] = info["tectonic_region"].lower()
    if info["tectonic_region"] not in TECTONIC_REGIONS:
        raise ValueError(
            f'source model "{name}" has unknown tectonic_region "{info["tectonic_region"]}"'
        )
    info["cluster"] = bool(info.get("cluster", False))
    info.setdefault("cluster_nshm_component", "FaultCluster")
    info["cluster_rupture_rate"] = str(info.get("cluster_rupture_rate", "conditional")).lower()
    if source_type == "point_source_models":
        info["point_source"] = normalize_point_source(info.get("point_source"))
        ps = info["point_source"]
        if (
            ps["distance_method"] == "nshmp"
            and ps.get("distance_bin") is None
            and name in NSHMP_GRID_OPTIMIZATION
        ):
            ps.update(NSHMP_GRID_OPTIMIZATION[name])
    return info


def weight_group(info, source_type):
    """Key of the group of alternative source models whose weights are normalized together.

    Source models in the same group are alternative logic tree branches (their weights are
    normalized to sum to one, and their weighted hazards are summed); hazards of different
    groups are summed. The default group is the source type, tectonic region, and NSHM
    component, which reproduces the earlier normalization of all fault source models (all
    active crust FaultSystem models) and all point source models (all active crust Grid
    models). source_info.json can set "logic_tree_group" explicitly.
    """
    group = info.get("logic_tree_group")
    if group:
        return (source_type, str(group).lower())
    return (source_type, info["tectonic_region"], str(info["nshm_component"]).lower())


# ---------------------------------------------------------------------------------------
# Ground motion model logic trees
# ---------------------------------------------------------------------------------------


class GmmBranch:
    """A ground motion model branch of a tectonic region."""

    def __init__(self, key, spec, weight, entry):
        self.key = key  # name in the config (lower case)
        self.spec = spec  # pygmm_gmms.GmmSpec
        self.weight = weight  # normalized weight
        self.entry = entry  # config dict of the model (weight normalized in place)

    def __repr__(self):
        return f"GmmBranch({self.key!r}, {self.spec.kind}:{self.spec.name}, {self.weight:.4g})"


def is_region_config(gmm_config):
    """True for the per-region format, False for the earlier flat format.

    Raises ValueError if the keys mix tectonic regions and ground motion models.
    """
    keys = list(gmm_config.keys())
    regions = [k for k in keys if k in TECTONIC_REGIONS]
    if regions and len(regions) != len(keys):
        raise ValueError(
            "ground_motion_models must either contain only tectonic regions "
            f"({', '.join(TECTONIC_REGIONS)}) or only ground motion models"
        )
    return bool(regions)


def _region_entries(gmm_config, region, notes):
    """Config entries {name: {"weight": ..}} of a region, or the default tree."""
    value = gmm_config.get(region)
    if isinstance(value, str):
        if value.lower() not in ("default", "nshm", "nshm23", "usgs"):
            raise ValueError(f'ground_motion_models["{region}"] must be an object or "default"')
        value = None
    if value is None:
        notes.append(f"{region}: using the default (USGS NSHM 2023) ground motion models")
        value = {k.lower(): {"weight": w} for k, w in DEFAULT_GMM_TREES[region].items()}
        gmm_config[region] = value
    return value


def build_region_tree(region, entries, notes):
    """Resolve and normalize the ground motion models of a region.

    Unavailable models (not in pygmm) are removed and the remaining weights renormalized,
    with a warning; their weights are set to 0 in the config entries. Weights are normalized
    in place in the config entries (as for the earlier flat format).
    """
    branches = []
    dropped = []
    for key, entry in entries.items():
        try:
            spec = pygmm_gmms.resolve(key, entry.get("options"), entry.get("scenario"))
        except pygmm_gmms.UnavailableGmmError as e:
            dropped.append((key, entry, e.args[0]))
            continue
        if spec.note:
            notes.append(f"{region}: {spec.note}")
        branches.append(GmmBranch(key, spec, entry.get("weight", 0.0), entry))
    for key, entry, message in dropped:
        if entry.get("weight", 0.0) > 0:
            notes.append(
                f"{region}: {message}; its weight ({entry['weight']:.4g}) is removed and the "
                "remaining weights are renormalized"
            )
        entry["weight"] = 0.0
    total = sum(b.weight for b in branches)
    if dropped and total <= 0 and any(e.get("weight", 0.0) for _, e, _ in dropped):
        raise ValueError(f"{region}: none of the ground motion models are available")
    if total > 0:
        for b in branches:
            b.weight = b.weight / total
            b.entry["weight"] = b.weight
    return branches


def parse_ground_motion_models(config, regions_used):
    """Build the ground motion model logic tree of every tectonic region that is used.

    Args:
        config (dict): config (lower case); config["ground_motion_models"] is updated in place
            with normalized weights (and with the default trees of regions that were not
            specified, in the per-region format)
        regions_used (iterable): tectonic regions of the source models with non-zero weight

    Returns:
        (trees, notes): trees is {region: [GmmBranch, ...]}; notes are warning messages
    """
    gmm_config = config["ground_motion_models"]
    notes = []
    trees = {}
    if is_region_config(gmm_config):
        for region in TECTONIC_REGIONS:
            if region in gmm_config or region in regions_used:
                entries = _region_entries(gmm_config, region, notes)
                trees[region] = build_region_tree(region, entries, notes)
    else:
        # Earlier flat format: the models are the active crust models
        trees["active_crust"] = build_region_tree("active_crust", gmm_config, notes)
        defaults = {}
        for region in TECTONIC_REGIONS:
            if region != "active_crust" and region in regions_used:
                entries = _region_entries(defaults, region, notes)
                trees[region] = build_region_tree(region, entries, notes)
    for note in notes:
        warnings.warn(note, UserWarning, stacklevel=3)
    return trees, notes


def source_model_gmm_entries(entry, info):
    """Ground motion model entries {name: {"weight": w}} specific to a source model, or None.

    The config entry of the source model can have "ground_motion_models" (same format as a
    tectonic region); otherwise a source_info.json "gmm_tree" list of {"id", "weight"} is
    used. Descriptive (string) gmm_tree values are ignored.
    """
    if isinstance(entry.get("ground_motion_models"), dict):
        return entry["ground_motion_models"]
    tree = info.get("gmm_tree")
    if isinstance(tree, list):
        entries = {}
        for branch in tree:
            key = str(branch["id"]).lower()
            entries[key] = {"weight": entries.get(key, {}).get("weight", 0.0) + float(branch["weight"])}
        return entries
    return None


def region_value(value, region):
    """A constraint given as a number or as {region: number}."""
    if isinstance(value, dict):
        return value.get(region)
    return value

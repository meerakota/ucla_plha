"""Ground motion models from pygmm, evaluated vectorized over the ruptures of a source model.

A ground motion model in a ucla_plha config is one of

* a ucla_plha model ("ask14", "bssa14", "cb14", "cy14", "idriss14"), evaluated by
  :func:`ucla_plha.plha.get_ground_motion_data` exactly as before;
* an nshmp-lib ``Gmm`` id (case insensitive), e.g. "NGA_EAST_2026",
  "NGA_EAST_SEEDS_2026_ADJUSTED", "KBCG_20_CASCADIA_INTERFACE_BASIN",
  "AG_20_CASCADIA_SLAB_ADJUSTED_BASIN", or "PSBAH_20_CASCADIA_SLAB_BASIN". The ids are read
  from the ``GMM_IDS`` attribute of the pygmm classes (e.g. ``NgaEast.GMM_IDS``); ids of
  pygmm classes that do not have a ``GMM_IDS`` attribute yet (AbrahamsonGulerce2020,
  KuehnEtAl2020, ParkerEtAl2020 at pygmm c69c417) come from :data:`FALLBACK_GMM_IDS`, which
  follows the tables in their docstrings. NGA-West2 ids that pygmm does not provide yet
  (e.g. "ASK_14_BASIN") are replaced by the ucla_plha model (:data:`SUBSTITUTES`) with a
  warning, and are used from pygmm automatically once pygmm provides them in ``GMM_IDS``.
  Ids of models that the installed pygmm does not provide (AM_09 and ZHAO_06 before pygmm
  1b30f5b; :data:`UNAVAILABLE`) are removed from the logic tree with a warning, and the
  remaining weights of the tectonic region are renormalized;
* a pygmm class name (case insensitive), e.g. "ParkerEtAl2020" or
  "AbrahamsonSilvaKamai2014" (the published pygmm models), with optional constructor
  ``options`` and ``scenario`` values in the config.

pygmm is imported only when a pygmm model is used, so ucla_plha does not require it
otherwise. Only PGA is computed (``ims=["pga"]``); the natural log of the median and the
total standard deviation are returned.
"""

import inspect
import logging
import warnings
from dataclasses import dataclass, field

import numpy as np

NATIVE_GMMS = ["ask14", "bssa14", "cb14", "cy14", "idriss14"]

# Distances used by the ucla_plha models (as in plha.get_source_data)
NATIVE_DISTANCES = {
    "bssa14": {"rjb"},
    "ask14": {"rrup", "rx"},
    "cb14": {"rjb", "rrup", "rx"},
    "cy14": {"rjb", "rrup", "rx"},
    "idriss14": {"rrup"},
}

# pygmm scenario keys of the distances computed by ucla_plha
PYGMM_DISTANCES = {"dist_jb": "rjb", "dist_rup": "rrup", "dist_x": "rx", "dist_y0": "ry0"}

# nshmp-lib Gmm ids of the NSHM trees that older versions of pygmm do not provide (pygmm added
# them in 1b30f5b). If the installed pygmm does not provide one of them, it is removed from
# the logic tree (with a warning) and the remaining weights of the tectonic region are
# renormalized.
UNAVAILABLE = {
    gmm_id: f"{name} is not provided by the installed pygmm"
    for name, ids in [
        (
            "Atkinson and Macias (2009)",
            [
                "AM_09_INTERFACE",
                "AM_09_INTERFACE_BASIN",
                "AM_09_INTERFACE_BASIN_M9",
                "AM_09_INTERFACE_BASIN_SITE_FIX",
                "AM_09_INTERFACE_BASIN_M9_SITE_FIX",
            ],
        ),
        (
            "Zhao et al. (2006)",
            [
                "ZHAO_06_INTERFACE",
                "ZHAO_06_INTERFACE_BASIN",
                "ZHAO_06_INTERFACE_BASIN_M9",
                "ZHAO_06_SLAB",
                "ZHAO_06_SLAB_BASIN",
            ],
        ),
    ]
    for gmm_id in ids
}

# nshmp-lib NGA-West2 ids that are replaced by the ucla_plha models when pygmm does not
# provide them (the nshmp-lib *_BASIN variants are being added to pygmm separately).
SUBSTITUTES = {
    f"{prefix}_14{suffix}": native
    for prefix, native in [("ASK", "ask14"), ("BSSA", "bssa14"), ("CB", "cb14"), ("CY", "cy14")]
    for suffix in ["", "_BASIN"]
}


def _subduction_ids(prefix, class_name, variants):
    """Build nshmp-lib ids for an NGA-Subduction model.

    variants: list of (id suffix after the event type, region, options); the id is
    f"{prefix}_{REGION}_{TYPE}{suffix}".
    """
    ids = {}
    for region, options, suffix, types in variants:
        for type_name in types:
            event_type = "interface" if type_name == "INTERFACE" else "intraslab"
            gmm_id = f"{prefix}_{region.upper()}_{type_name}{suffix}"
            ids[gmm_id] = (class_name, dict(options), dict(event_type=event_type, region=region))
    return ids


BOTH = ["INTERFACE", "SLAB"]
#: nshmp-lib ids of pygmm classes without a GMM_IDS attribute: id -> (class name, constructor
#: options, scenario values). From the docstrings of the pygmm classes (pygmm c69c417).
FALLBACK_GMM_IDS = {
    **_subduction_ids(
        "AG_20",
        "AbrahamsonGulerce2020",
        [
            ("global", {}, "", BOTH),
            ("global", {"epistemic": False}, "_NO_EPI", BOTH),
            ("global", {"ak_adjusted": True}, "_AK_ADJUSTED", ["INTERFACE"]),
            ("cascadia", {}, "", BOTH),
            ("cascadia", {"basin": True}, "_BASIN", BOTH),
            ("cascadia", {"adjusted": True}, "_ADJUSTED", BOTH),
            ("cascadia", {"adjusted": True, "basin": True}, "_ADJUSTED_BASIN", BOTH),
            ("alaska", {}, "", BOTH),
            ("alaska", {"adjusted": True}, "_ADJUSTED", BOTH),
            ("prvi", {}, "", BOTH),
        ],
    ),
    **_subduction_ids(
        "KBCG_20",
        "KuehnEtAl2020",
        [
            ("global", {}, "", BOTH),
            ("cascadia", {}, "", BOTH),
            ("alaska", {}, "", BOTH),
            ("prvi", {}, "", BOTH),
            ("global", {"epistemic": False}, "_NO_EPI", BOTH),
            ("global", {"ak_adjusted": True}, "_AK_ADJUSTED", ["INTERFACE"]),
            ("cascadia", {"basin": True}, "_BASIN", BOTH),
            ("cascadia", {"basin": True, "seattle_basin": True}, "_SEATTLE_BASIN", ["INTERFACE"]),
            (
                "cascadia",
                {"basin": True, "seattle_basin": True, "m9": True},
                "_SEATTLE_BASIN_M9",
                ["INTERFACE"],
            ),
            ("cascadia", {"seattle_basin": True}, "_SEATTLE_BASIN", ["SLAB"]),
        ],
    ),
    **_subduction_ids(
        "PSBAH_20",
        "ParkerEtAl2020",
        [
            ("global", {}, "", BOTH),
            ("global", {"epistemic": False}, "_NO_EPI", BOTH),
            ("global", {"ak_adjusted": True}, "_AK_ADJUSTED", ["INTERFACE"]),
            ("alaska", {}, "", BOTH),
            ("cascadia", {}, "", BOTH),
            ("cascadia", {"basin": True}, "_BASIN", BOTH),
            ("cascadia", {"basin": True, "m9": True}, "_BASIN_M9", ["INTERFACE"]),
            ("prvi", {}, "", BOTH),
        ],
    ),
}

# Scenario values implied by the tectonic region of the source model
REGION_EVENT_TYPE = {"subduction_interface": "interface", "subduction_slab": "intraslab"}

# ucla_plha fault type (1 reverse, 2 normal, 3 strike slip) -> pygmm mechanism
MECHANISM = np.array(["U", "RS", "NS", "SS"])


@dataclass
class GmmSpec:
    """A resolved ground motion model."""

    key: str  # name as given in the config (lower case)
    kind: str  # "native" or "pygmm"
    name: str  # native model name, nshmp-lib id, or pygmm class name
    cls: object = None  # pygmm class
    options: dict = field(default_factory=dict)  # constructor keyword arguments
    scenario: dict = field(default_factory=dict)  # fixed scenario values
    note: str = ""  # substitution note

    @property
    def distances(self):
        """Distances ("rjb", "rrup", "rx", "ry0") used by the model."""
        if self.kind == "native":
            return NATIVE_DISTANCES[self.name]
        names = {p.name for p in self.cls.PARAMS}
        return {PYGMM_DISTANCES[k] for k in names if k in PYGMM_DISTANCES}


class UnavailableGmmError(KeyError):
    """An nshmp-lib Gmm id of a model that is not implemented in pygmm."""


_REGISTRY = None


def _pygmm():
    try:
        import pygmm
    except ImportError as e:  # pragma: no cover - depends on the environment
        raise ImportError(
            "pygmm is required for the nshmp-lib and pygmm ground motion models "
            "(pip install pygmm, or add the pygmm source directory to PYTHONPATH)"
        ) from e
    return pygmm


def _split_options(cls, values):
    """Split a GMM_IDS entry into constructor options and scenario values."""
    params = set(inspect.signature(cls.__init__).parameters)
    options = {k: v for k, v in values.items() if k in params}
    scenario = {k: v for k, v in values.items() if k not in params}
    return options, scenario


def registry():
    """Return {"ids": {ID: (cls, options, scenario)}, "classes": {lower name: cls}}."""
    global _REGISTRY
    if _REGISTRY is not None:
        return _REGISTRY
    pygmm = _pygmm()
    classes = {}
    ids = {}
    for name in getattr(pygmm, "__all__", dir(pygmm)):
        cls = getattr(pygmm, name, None)
        if not inspect.isclass(cls) or not hasattr(cls, "PARAMS") or name == "Scenario":
            continue
        classes[name.lower()] = cls
        for gmm_id, values in (getattr(cls, "GMM_IDS", None) or {}).items():
            options, scenario = _split_options(cls, dict(values))
            ids[gmm_id.upper()] = (cls, options, scenario)
    for gmm_id, (class_name, options, scenario) in FALLBACK_GMM_IDS.items():
        cls = getattr(pygmm, class_name, None)
        if cls is not None and gmm_id not in ids:
            ids[gmm_id] = (cls, options, scenario)
    _REGISTRY = {"ids": ids, "classes": classes}
    return _REGISTRY


def _scenario_from_id(gmm_id, cls):
    """Event type and region implied by an nshmp-lib id (for GMM_IDS without them)."""
    names = {p.name for p in cls.PARAMS}
    scenario = {}
    if "event_type" in names:
        if "_INTERFACE" in gmm_id:
            scenario["event_type"] = "interface"
        elif "_SLAB" in gmm_id:
            scenario["event_type"] = "intraslab"
    if "region" in names:
        for region in ["cascadia", "alaska", "prvi", "global"]:
            if f"_{region.upper()}_" in gmm_id + "_":
                scenario["region"] = region
                break
    return scenario


def resolve(key, options=None, scenario=None):
    """Resolve a ground motion model name from a config into a :class:`GmmSpec`.

    Raises:
        UnavailableGmmError: for nshmp-lib ids of NSHM models that the installed pygmm does
            not provide
        ValueError: for unknown names
    """
    key = key.lower()
    if key in NATIVE_GMMS:
        return GmmSpec(key=key, kind="native", name=key)
    gmm_id = key.upper()
    try:
        reg = registry()
    except ImportError:
        if gmm_id in SUBSTITUTES or gmm_id in UNAVAILABLE:
            reg = {"ids": {}, "classes": {}}
        else:
            raise
    if gmm_id in reg["ids"]:
        cls, opts, scen = reg["ids"][gmm_id]
        scen = {**_scenario_from_id(gmm_id, cls), **scen, **(scenario or {})}
        return GmmSpec(
            key=key, kind="pygmm", name=gmm_id, cls=cls, options={**opts, **(options or {})},
            scenario=scen,
        )
    if gmm_id in UNAVAILABLE:
        raise UnavailableGmmError(f"{gmm_id}: {UNAVAILABLE[gmm_id]}")
    if gmm_id in SUBSTITUTES:
        native = SUBSTITUTES[gmm_id]
        return GmmSpec(
            key=key, kind="native", name=native,
            note=f"{gmm_id} is not available in pygmm; using the ucla_plha model {native}",
        )
    if key in reg["classes"]:
        return GmmSpec(
            key=key, kind="pygmm", name=reg["classes"][key].__name__, cls=reg["classes"][key],
            options=dict(options or {}), scenario=dict(scenario or {}),
        )
    raise ValueError(
        f'unknown ground motion model "{key}": expected a ucla_plha model ({", ".join(NATIVE_GMMS)}), '
        "an nshmp-lib Gmm id provided by pygmm, or a pygmm class name"
    )


def build_scenario(spec, rupture, site, tectonic_region=None):
    """Map ucla_plha rupture and site arrays to the pygmm scenario values used by a model.

    Args:
        spec (GmmSpec): resolved pygmm model
        rupture (dict): arrays "m", "fault_type", "rjb", "rrup", "rx", "ry0", "dip", "ztor",
            "zbor" (length N)
        site (dict): "vs30", "measured_vs30", "z1p0", "z2p5", "zsed" (km; None if unknown)
        tectonic_region (str): tectonic region of the source model, used for the event type
            of the NGA-Subduction models when the model id does not imply it

    Returns:
        dict of scenario values
    """
    names = {p.name for p in spec.cls.PARAMS}
    m = np.asarray(rupture["m"], dtype=float)
    values = {"mag": m, "v_s30": float(site["vs30"])}
    for key, name in PYGMM_DISTANCES.items():
        values[key] = rupture.get(name)
    dip = rupture.get("dip")
    ztor = rupture.get("ztor")
    zbor = rupture.get("zbor")
    values["dip"] = dip
    values["depth_tor"] = ztor
    values["depth_bor"] = zbor
    if dip is not None and ztor is not None and zbor is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            values["width"] = (np.asarray(zbor) - np.asarray(ztor)) / np.sin(np.radians(dip))
    if rupture.get("fault_type") is not None:
        values["mechanism"] = MECHANISM[np.asarray(rupture["fault_type"], dtype=int)]
    if rupture.get("rx") is not None:
        values["on_hanging_wall"] = np.asarray(rupture["rx"]) >= 0.0
    values["vs_source"] = "measured" if site.get("measured_vs30") else "inferred"
    # Site terms: None means the model default (reference basin depth, no coastal plain)
    values["depth_1_0"] = site.get("z1p0")
    values["depth_2_5"] = site.get("z2p5")
    values["depth_sed"] = site.get("zsed")
    if tectonic_region in REGION_EVENT_TYPE:
        values["event_type"] = REGION_EVENT_TYPE[tectonic_region]
    values.update(spec.scenario)
    return {k: v for k, v in values.items() if k in names and v is not None}


def get_ground_motion(spec, rupture, site, tectonic_region=None):
    """Mean and standard deviation of ln(PGA) of a pygmm model for every rupture.

    Returns:
        (mu_ln_pga, sigma_ln_pga): arrays of length N
    """
    pygmm = _pygmm()
    n = len(rupture["m"])
    if n == 0:
        return np.empty(0), np.empty(0)
    values = build_scenario(spec, rupture, site, tectonic_region)
    # pygmm warns (with warnings and logging) about values outside of the recommended limits
    # of the models, which is expected for the full set of ruptures of a source model
    previous = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            model = spec.cls(pygmm.Scenario(**values), ims=["pga"], **spec.options)
            mu = np.asarray(model.ln_pga, dtype=float)
            sigma = np.asarray(model.ln_std_pga, dtype=float)
    finally:
        logging.disable(previous)
    mu = np.broadcast_to(mu, (n,)).astype(float)
    sigma = np.broadcast_to(sigma, (n,)).astype(float)
    return mu, sigma


def list_gmm_ids():
    """All nshmp-lib ids available from pygmm, sorted."""
    return sorted(registry()["ids"])


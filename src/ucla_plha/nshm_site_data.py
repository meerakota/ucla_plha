"""Location-dependent site data of the USGS NSHM (conterminous U.S. 2023) stable crust models.

nshmp-lib (``model.SiteData`` and ``calc.Sites``) gives every site of a hazard calculation two
location-dependent values from the nshm-conus ``site-data/`` directory, including the sites of
the USGS hazard web service, which take only a longitude, latitude, and Vs30:

- the Gulf and Atlantic coastal plain sediment thickness, zSed (km; Boyd, 2023), which the
  NGA-East models use for the coastal plain taper of the 2023 adjustment and for the Chapman and
  Guo (2021) coastal plain amplification, and
- the "Coastal Plain CPA region" ground motion model logic tree, which replaces the stable crust
  logic trees of all stable crust sources at sites inside the region and adds the coastal plain
  amplification models (NGA_EAST_2026_CPA and NGA_EAST_SEEDS_2026_CPA).

The data are in ``site_data/coastal_plain.npz`` (made by
``utilities/convert_nshm23_site_data.py`` from nshm-conus 6.2.0).
"""

import functools
import json
from importlib.resources import files

import numpy as np


@functools.lru_cache(maxsize=1)
def _data():
    with np.load(str(files("ucla_plha").joinpath("site_data", "coastal_plain.npz"))) as d:
        return {
            "zsed_m": d["zsed_m"],
            "lon0": float(d["lon0"]),
            "lat0": float(d["lat0"]),
            "spacing": float(d["spacing"]),
            "margin_polygon": d["margin_polygon"],
            "region_polygon": d["region_polygon"],
            "region_tree": json.loads(str(d["region_tree"])),
            "info": json.loads(str(d["info"])),
        }


def _contains(polygon, lon, lat):
    """Even-odd (ray casting) test of whether a polygon (n, 2) contains a point."""
    x, y = polygon[:, 0], polygon[:, 1]
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    crosses = (y > lat) != (y2 > lat)
    with np.errstate(divide="ignore", invalid="ignore"):
        x_cross = x + (lat - y) * (x2 - x) / (y2 - y)
    return bool(np.count_nonzero(crosses & (lon < x_cross)) % 2)


def _snap(value, spacing):
    # nshmp-lib SiteData.snapToGrid: Math.round(value / spacing) * spacing, rounded HALF_UP
    # to the decimals of the spacing
    return round(float(np.floor(value / spacing + 0.5)) * spacing, 2)


def coastal_plain_zsed(longitude, latitude):
    """Coastal plain sediment thickness (km) of the NSHM at a site, or None off the plain.

    As nshmp-lib: the site is snapped to the 0.05 degree grid, and the value is used if the
    snapped location is inside the margin polygon and has data.
    """
    d = _data()
    s = d["spacing"]
    lon, lat = _snap(longitude, s), _snap(latitude, s)
    if not _contains(d["margin_polygon"], lon, lat):
        return None
    i = int(round((lat - d["lat0"]) / s))
    j = int(round((lon - d["lon0"]) / s))
    grid = d["zsed_m"]
    if not (0 <= i < grid.shape[0] and 0 <= j < grid.shape[1]) or grid[i, j] < 0:
        return None
    return float(grid[i, j]) / 1000.0


def in_coastal_plain_region(longitude, latitude):
    """True if the site is inside the NSHM "Coastal Plain CPA region" (gmm-region polygon)."""
    return _contains(_data()["region_polygon"], float(longitude), float(latitude))


def coastal_plain_stable_crust_tree():
    """{nshmp-lib Gmm id: weight} of the stable crust tree of the Coastal Plain CPA region."""
    tree = _data()["region_tree"]["stable-crust"]["all"]
    return {branch["id"]: float(branch["weight"]) for branch in tree}


def region_name():
    return _data()["info"]["region"]["name"]

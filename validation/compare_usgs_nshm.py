"""
Compare ucla_plha PGA hazard curves with the USGS NSHM (conterminous U.S. 2023) hazard web
service, per NSHM source component (Fault, FaultCluster, FaultSystem, Grid, Interface, Slab,
Zone).

The USGS service (https://earthquake.usgs.gov/ws/nshmp/conus-2023/dynamic/hazard/<lon>/<lat>/<vs30>)
returns the hazard curves of each source component computed by nshmp-lib. ucla_plha computes
the hazard of every source model whose source_info.json gives its nshm_component, with the
nshmp-lib settings that ucla_plha supports: upper truncation of the ground motion
distribution at 3 sigma, the NSHM maximum distance of each tectonic region, and the NSHM
ground motion model logic trees (or the alternatives below for the active crust).

Usage (from the repository root, with pygmm on PYTHONPATH):

    python validation/compare_usgs_nshm.py --site seattle --active-gmms nshmp
    python validation/compare_usgs_nshm.py --lon -118.25 --lat 34.05 --active-gmms ucla_plha
    python validation/compare_usgs_nshm.py --site all --active-gmms all --out results.json

Options:
    --active-gmms: active crust ground motion models: "ucla_plha" (ucla_plha ask14, bssa14,
        cb14, cy14, 0.25 each), "pygmm" (the published pygmm NGA-West2 models,
        AbrahamsonSilvaKamai2014 etc.), "nshmp" (the nshmp-lib ASK/BSSA/CB/CY_14_BASIN models
        from pygmm, the NSHM tree), or "all" (all three)
    --grid-distances: "source_info" (default; the treatment in each model's source_info.json,
        i.e. the nshmp-lib FINITE point sources for nshm23_wus_grid) or "ucla_plha" (use the
        ucla_plha crustal approximations recorded in ucla_plha_point_source of
        nshm23_wus_grid, the treatment before the nshmp-lib settings became the default)
    --models: source models to use (default: every source model with a name starting with
        nshm23_)
    --zsed: coastal plain sediment thickness (km); by default ucla_plha uses the NSHM value at
        the site, as the USGS service does
    --no-nshm-site-data: do not use the NSHM site data (zSed and the Coastal Plain CPA region
        ground motion models of the stable crust)
    --cache: directory for the service responses (default validation/usgs_nshm_cache)

The service responses are cached, so the comparison can be repeated offline.
"""

import argparse
import contextlib
import json
import os
import subprocess
import sys
import tempfile
import time
import warnings
from importlib.resources import files

import numpy as np

from ucla_plha import plha, tectonic_regions

SERVICE = "https://earthquake.usgs.gov/ws/nshmp/conus-2023/dynamic/hazard/{lon}/{lat}/{vs30}"

SITES = {
    "los_angeles": (-118.25, 34.05),
    "bay_area": (-122.25, 37.80),
    "seattle": (-122.30, 47.60),
    "portland": (-122.68, 45.52),
    "coastal_oregon": (-124.05, 44.63),
    "northern_california_coast": (-124.16, 40.80),
    "memphis": (-90.05, 35.15),
    "new_madrid": (-89.60, 36.60),
    "charleston": (-79.94, 32.78),
    "central_virginia": (-77.90, 37.80),
    "quiet_ceus": (-99.00, 41.00),
}

ACTIVE_GMMS = {
    "ucla_plha": {"ask14": 0.25, "bssa14": 0.25, "cb14": 0.25, "cy14": 0.25},
    "pygmm": {
        "abrahamsonsilvakamai2014": 0.25,
        "boorestewartseyhanatkinson2014": 0.25,
        "campbellbozorgnia2014": 0.25,
        "chiouyoungs2014": 0.25,
    },
    "nshmp": {k.lower(): w for k, w in tectonic_regions.DEFAULT_GMM_TREES["active_crust"].items()},
}

RETURN_PERIODS = [475.0, 2475.0]


def fetch_service(lon, lat, vs30, cache):
    """USGS hazard service response for a site (cached in cache/)."""
    os.makedirs(cache, exist_ok=True)
    filename = os.path.join(cache, f"conus2023_{lon:.3f}_{lat:.3f}_{vs30:.0f}.json")
    if not os.path.exists(filename):
        url = SERVICE.format(lon=f"{lon:g}", lat=f"{lat:g}", vs30=f"{vs30:g}")
        text = subprocess.run(
            ["curl", "-s", "-L", url], capture_output=True, text=True, check=True
        ).stdout
        data = json.loads(text)
        if data.get("status") != "success":
            raise RuntimeError(f"USGS service error for {url}: {text[:500]}")
        # keep only the PGA curves (the comparison uses PGA)
        data["response"]["hazardCurves"] = [
            c for c in data["response"]["hazardCurves"] if c["imt"]["value"] == "PGA"
        ]
        with open(filename, "w") as f:
            json.dump(data, f)
    with open(filename) as f:
        data = json.load(f)
    for curves in data["response"]["hazardCurves"]:
        if curves["imt"]["value"] == "PGA":
            return {
                c["component"]: (
                    np.array(c["values"]["xs"], dtype=float),
                    np.array(c["values"]["ys"], dtype=float),
                )
                for c in curves["data"]
            }
    raise RuntimeError("no PGA curves in the service response")


def available_models(prefix="nshm23_"):
    """{name: (source_type, source_info)} of the source models in the package."""
    root = str(files("ucla_plha").joinpath("source_models"))
    models = {}
    for source_type in ["fault_source_models", "point_source_models"]:
        for name in sorted(os.listdir(os.path.join(root, source_type))):
            if name.startswith(prefix) and os.path.isdir(os.path.join(root, source_type, name)):
                models[name] = (source_type, plha.get_source_info(source_type, name))
    return models


@contextlib.contextmanager
def ucla_plha_grid_distances(enabled):
    """Temporarily use the ucla_plha crustal point source distances of models that record them."""
    if not enabled:
        yield
        return
    original = plha.get_source_info

    def patched(source_type, source_model):
        info = original(source_type, source_model)
        if "ucla_plha_point_source" in info:
            info = dict(info)
            info["point_source"] = tectonic_regions.normalize_point_source(
                info["ucla_plha_point_source"]
            )
        return info

    plha.get_source_info = patched
    try:
        yield
    finally:
        plha.get_source_info = original


def nshm_config(lon, lat, vs30, pga, models, active_gmms, z1p0=None, z2p5=None, zsed=None,
                nshm_site_data=True):
    """ucla_plha config with the NSHM settings (see the module docstring)."""
    site = {"latitude": lat, "longitude": lon, "elevation": 0.0, "vs30": vs30}
    for key, value in (("z1p0", z1p0), ("z2p5", z2p5), ("zsed", zsed)):
        if value is not None:
            site[key] = value
    if not nshm_site_data:
        site["nshm_site_data"] = False
    source_models = {"fault_source_models": {}, "point_source_models": {}}
    for name, (source_type, info) in models.items():
        source_models[source_type][name] = {"weight": 1.0}
        # model-specific NSHM maximum distance (e.g. 300 km for the stable crust system grid)
        cutoff = info.get("gmm_max_distance_km")
        if cutoff is not None:
            source_models[source_type][name]["dist_cutoff"] = float(cutoff)
    source_models = {k: v for k, v in source_models.items() if v}
    gmms = {"active_crust": {k: {"weight": w} for k, w in ACTIVE_GMMS[active_gmms].items()}}
    for region in tectonic_regions.TECTONIC_REGIONS[1:]:
        gmms[region] = "default"
    config = {
        "site": site,
        "constraints": {
            "dist_cutoff": dict(tectonic_regions.NSHM_MAX_DISTANCE),
            "truncation_level": 3.0,
        },
        "source_models": source_models,
        "ground_motion_models": gmms,
        "output": {"psha": {"pga": list(map(float, pga)), "source_model_hazard": True}},
    }
    return config


def run_ucla_plha(lon, lat, vs30, pga, models, active_gmms, grid_distances="source_info",
                  z1p0=None, z2p5=None, zsed=None, nshm_site_data=True):
    """ucla_plha hazard per NSHM component.

    Returns ({component: curve}, {source model key: curve}, notes, seconds)
    """
    config = nshm_config(lon, lat, vs30, pga, models, active_gmms, z1p0, z2p5, zsed,
                         nshm_site_data)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "config.json")
        with open(path, "w") as f:
            json.dump(config, f)
        start = time.time()
        with warnings.catch_warnings(), ucla_plha_grid_distances(grid_distances == "ucla_plha"):
            warnings.simplefilter("ignore")
            out = plha.get_hazard(path)
        seconds = time.time() - start
    by_model = out["output"]["psha"]["source_model_hazard"]
    components = {}
    for key, value in by_model.items():
        curve = np.array(value["annual_rate_of_exceedance"])
        comp = value["nshm_component"]
        components[comp] = components.get(comp, 0.0) + curve
    return components, by_model, out.get("notes", []), seconds


def ground_motion_at(xs, ys, return_period):
    """PGA (g) at an annual rate of 1/return_period, log-log interpolation; NaN if outside."""
    rate = 1.0 / return_period
    ok = ys > 0
    xs, ys = xs[ok], ys[ok]
    if len(ys) < 2 or rate > ys[0] or rate < ys[-1]:
        return float("nan")
    return float(np.exp(np.interp(np.log(rate), np.log(ys[::-1]), np.log(xs[::-1]))))


def compare_site(name, lon, lat, vs30, active_gmms, models, cache, grid_distances,
                 z1p0=None, z2p5=None, zsed=None, nshm_site_data=True):
    usgs = fetch_service(lon, lat, vs30, cache)
    xs = usgs["Total"][0]
    components, by_model, notes, seconds = run_ucla_plha(
        lon, lat, vs30, xs, models, active_gmms, grid_distances, z1p0, z2p5, zsed,
        nshm_site_data,
    )
    present = {info["nshm_component"] for _, info in models.values()} | {
        info["cluster_nshm_component"] for _, info in models.values() if info["cluster"]
    }
    rows = []
    for comp in sorted(set(usgs) - {"Total"}):
        if comp not in present:
            continue
        ys_usgs = usgs[comp][1]
        ys_ucla = np.asarray(components.get(comp, np.zeros(len(xs))))
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = ys_ucla / ys_usgs
        row = {
            "component": comp,
            "xs": xs.tolist(),
            "usgs": ys_usgs.tolist(),
            "ucla_plha": ys_ucla.tolist(),
            "ratio": [float(r) if np.isfinite(r) else None for r in ratio],
        }
        for rp in RETURN_PERIODS:
            g_usgs = ground_motion_at(xs, ys_usgs, rp)
            g_ucla = ground_motion_at(xs, ys_ucla, rp)
            row[f"pga_{int(rp)}yr_usgs"] = g_usgs
            row[f"pga_{int(rp)}yr_ucla_plha"] = g_ucla
        rows.append(row)
    return {
        "site": name,
        "longitude": lon,
        "latitude": lat,
        "vs30": vs30,
        "active_gmms": active_gmms,
        "grid_distances": grid_distances,
        "models": sorted(models),
        "seconds": seconds,
        "notes": notes,
        "components": rows,
        "source_models": by_model,
    }


def print_result(result):
    print(
        f"\n{result['site']} ({result['longitude']}, {result['latitude']}), Vs30 = {result['vs30']}, "
        f"active crust GMMs: {result['active_gmms']}, grid distances: {result['grid_distances']} "
        f"({result['seconds']:.1f} s)"
    )
    header = f"  {'component':<13}{'PGA 475 yr (USGS / ucla)':>28}{'PGA 2475 yr (USGS / ucla)':>30}   ratio of rates at PGA = 0.01, 0.1, 0.4, 1.0 g"
    print(header)
    for row in result["components"]:
        xs = np.array(row["xs"])
        ratios = []
        for g in [0.01, 0.1, 0.4, 1.0]:
            i = int(np.argmin(np.abs(np.log(xs / g))))
            r = row["ratio"][i]
            ratios.append("   -" if r is None else f"{r:5.2f}")
        a = f"{row['pga_475yr_usgs']:.3f} / {row['pga_475yr_ucla_plha']:.3f} ({row['pga_475yr_ucla_plha'] / row['pga_475yr_usgs'] - 1:+.0%})"
        b = f"{row['pga_2475yr_usgs']:.3f} / {row['pga_2475yr_ucla_plha']:.3f} ({row['pga_2475yr_ucla_plha'] / row['pga_2475yr_usgs'] - 1:+.0%})"
        print(f"  {row['component']:<13}{a:>28}{b:>30}   {'  '.join(ratios)}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--site", default="seattle", help="site name (" + ", ".join(SITES) + ") or all")
    parser.add_argument("--lon", type=float)
    parser.add_argument("--lat", type=float)
    parser.add_argument("--vs30", type=float, default=760.0)
    parser.add_argument("--z1p0", type=float)
    parser.add_argument("--z2p5", type=float)
    parser.add_argument("--zsed", type=float)
    parser.add_argument("--no-nshm-site-data", action="store_true")
    parser.add_argument("--active-gmms", default="nshmp", choices=list(ACTIVE_GMMS) + ["all"])
    parser.add_argument("--grid-distances", default="source_info", choices=["source_info", "ucla_plha"])
    parser.add_argument("--models", nargs="*")
    parser.add_argument("--cache", default=os.path.join(os.path.dirname(__file__), "usgs_nshm_cache"))
    parser.add_argument("--out", help="write the results to this JSON file")
    args = parser.parse_args(argv)

    if args.lon is not None and args.lat is not None:
        sites = {"custom": (args.lon, args.lat)}
    elif args.site == "all":
        sites = SITES
    else:
        sites = {args.site: SITES[args.site]}
    models = available_models()
    if args.models:
        models = {k: v for k, v in models.items() if k in args.models}
    gmm_sets = list(ACTIVE_GMMS) if args.active_gmms == "all" else [args.active_gmms]

    results = []
    for name, (lon, lat) in sites.items():
        for gmm_set in gmm_sets:
            result = compare_site(
                name, lon, lat, args.vs30, gmm_set, models, args.cache, args.grid_distances,
                args.z1p0, args.z2p5, args.zsed, not args.no_nshm_site_data,
            )
            print_result(result)
            results.append(result)
    if args.out:
        with open(args.out, "w") as f:
            json.dump(results, f, indent=1)
    return results


if __name__ == "__main__":
    sys.exit(0 if main() else 1)

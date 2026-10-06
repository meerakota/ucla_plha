"""
Compare the ucla_plha PGA disaggregation with the USGS NSHM (conterminous U.S. 2023)
disaggregation web service, per NSHM source component and per source (contributor).

The USGS service
(https://earthquake.usgs.gov/ws/nshmp/conus-2023/dynamic/disagg/<lon>/<lat>/<vs30>/<returnPeriod>?imt=PGA,
or .../<lon>/<lat>/<vs30>?PGA=<iml> for a ground motion) returns the nshmp-lib
disaggregation (calc/Disaggregator.java): the contribution and the mean magnitude, distance,
and epsilon of every source type component and of every source contributing 1% or more, and
the magnitude-distance-epsilon distribution (m 4.4-9.4 by 0.2, rRup 0-1000 km by 20 km, and
epsilon bins of 0.5 from -2.5 to 2.5 with open outer bins).

nshmp-lib conventions (Disaggregator.processSource): every rupture of every source type is
binned and averaged by its rRup (HazardInput.rRup), epsilon is (ln IML - mu) / sigma of each
ground motion model branch (untruncated, while the rate uses the truncated exceedance
probability), and the contributions of the ruptures of a cluster are scaled to the cluster
rate at the IML (the cluster curve of each ground motion model, so the hazard of a cluster is
shared over all branches of the model). ucla_plha uses the same IML (the USGS target of the
return period), the NSHM settings of compare_usgs_nshm.py, the NSHM bins, and the disaggregation
options "distance_metric": "rrup" (--distance-metric) and "cluster_attribution": "gmm"
(--cluster-attribution). The ucla_plha distribution is binned here as nshmp-lib bins (see
nshmp_matrix); the get_hazard bins (lower edges included) are compared for the total too.

Contributors: the service lists the sources contributing 1% or more. ucla_plha groups
nshm23_ceus_fault by fault (NSHM source tree, all leaves together), the clusters of
nshm23_ceus_fault_cluster by cluster (the recurrence-rate branches together), and
nshm23_ceus_zone by zone (zones with merged nodes together); the USGS contributors are merged
the same way, so a ucla_plha group can be larger than the USGS one when some of its leaves or
branches are below 1%.

Usage (from the repository root, with pygmm on PYTHONPATH):

    python validation/compare_usgs_disagg.py --site new_madrid
    python validation/compare_usgs_disagg.py --site new_madrid --iml 0.8
    python validation/compare_usgs_disagg.py --site ceus --out results.json
    python validation/compare_usgs_disagg.py --site bay_area --grid-distances nshmp

--site: a site of compare_usgs_nshm.SITES, "ceus" (the CEUS sites), or "all"; --return-periods
(default 475 2475) or --iml; --models (default: the nshm23_ models). The service responses are
reduced and cached in validation/usgs_nshm_cache/ (disagg_*.json), so the comparison can be
repeated offline.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_usgs_nshm as cun  # noqa: E402

from ucla_plha import plha  # noqa: E402

SERVICE = "https://earthquake.usgs.gov/ws/nshmp/conus-2023/dynamic/disagg/{lon}/{lat}/{vs30}"
CEUS_SITES = ["new_madrid", "memphis", "charleston", "central_virginia", "quiet_ceus"]

# nshmp-lib default disaggregation bins (DisaggDataset / service "Discretization")
M_EDGES = np.round(np.arange(4.4, 9.4 + 1e-9, 0.2), 2)
R_EDGES = np.arange(0.0, 1000.0 + 1e-9, 20.0)
EPS_EDGES = np.r_[-1.0e9, np.arange(-2.5, 2.5 + 1e-9, 0.5), 1.0e9]

COMPONENTS = ["Fault", "FaultCluster", "FaultSystem", "Grid", "Zone", "Interface", "Slab"]


# ---------------------------------------------------------------------------------------
# USGS service
# ---------------------------------------------------------------------------------------


def _summary(component, name):
    for block in component["summary"]:
        if block["name"].startswith(name):
            return {item["name"]: item["value"] for item in block["data"]}
    return {}


def _reduce(data):
    """Keep the parts of a service response used here (rounded), to keep the cache small."""
    disagg = next(d for d in data["response"]["disaggs"] if d["imt"]["value"] == "PGA")
    out = {"url": data["url"], "request": data["request"], "components": {}}
    for comp in disagg["data"]:
        name = comp["component"].replace("Source Type: ", "")
        target = _summary(comp, "Disaggregation targets")
        totals = _summary(comp, "Totals")
        mean = _summary(comp, "Mean")
        mode = _summary(comp, "Mode (largest m-r bin)")
        out["iml"] = target.get("PGA ground motion")
        out["rate"] = target.get("Exceedance rate")
        out["recovered_rate"] = _summary(comp, "Recovered targets").get("Exceedance rate")
        bins = []
        # the Total bins are the sums of the bins of the source type components
        for rm in (comp["data"] or []) if name != "Total" else []:
            for e in rm["εdata"]:
                if e["value"] >= 1e-4:
                    bins.append([rm["r"], rm["m"], e["εbin"], float(f"{e['value']:.4g}")])
        entry = {
            "binned": totals.get("Binned"),
            "residual": totals.get("Residual"),
            "mean": [mean.get("m"), mean.get("r"), mean.get("ε₀")],
            "mode_mr": [mode.get("m"), mode.get("r"), mode.get("ε₀"), mode.get("Contribution")],
            "bins": bins,
        }
        if name == "Total":
            entry["sources"] = [
                [s["name"], s["source"], s["type"], s["contribution"],
                 float(f"{s['m']:.4g}"), float(f"{s['r']:.4g}"), float(f"{s['ε']:.3g}"),
                 s["latitude"], s["longitude"]]
                for s in comp["sources"]
                if s["type"] == "SET"
            ]
        out["components"][name] = entry
    return out


def fetch_disagg(lon, lat, vs30, cache, return_period=None, iml=None):
    """Reduced USGS disaggregation response (cached in cache/)."""
    os.makedirs(cache, exist_ok=True)
    target = f"rp{return_period:g}" if return_period else f"pga{iml:g}"
    filename = os.path.join(cache, f"disagg_{lon:.3f}_{lat:.3f}_{vs30:.0f}_{target}.json")
    if not os.path.exists(filename):
        url = SERVICE.format(lon=f"{lon:g}", lat=f"{lat:g}", vs30=f"{vs30:g}")
        if return_period:
            url += f"/{return_period:g}?imt=PGA&out=SOURCE,DISAGG_DATA"
        else:
            url += f"?PGA={iml:g}&out=SOURCE,DISAGG_DATA"
        text = subprocess.run(
            ["curl", "-s", "-L", url], capture_output=True, check=True
        ).stdout.decode("utf-8")
        data = json.loads(text)
        if data.get("status") != "success":
            raise RuntimeError(f"USGS service error for {url}: {text[:500]}")
        with open(filename, "w", encoding="utf-8") as f:
            json.dump(_reduce(data), f, ensure_ascii=False, separators=(",", ":"))
    with open(filename, encoding="utf-8") as f:
        return json.load(f)


def usgs_matrix(bins):
    """M x R x E array (%) of the reduced service bins."""
    out = np.zeros((len(M_EDGES) - 1, len(R_EDGES) - 1, len(EPS_EDGES) - 1))
    for r, m, e, value in bins:
        i = int(np.searchsorted(M_EDGES, m) - 1)
        j = int(np.searchsorted(R_EDGES, r) - 1)
        out[i, j, e] += value
    return out


def usgs_total_bins(usgs):
    return [b for c, v in usgs["components"].items() if c != "Total" for b in v["bins"]]


def usgs_groups(usgs, models):
    """{(component, group): [contribution %, m, r, eps]} of the USGS contributors (sources
    contributing 1% or more), merged as ucla_plha merges them (see _Recorder.groups):
    fault source tree leaves into their fault ("New Madrid - USGS (center)" -> "New Madrid"),
    cluster recurrence-rate branches (": R1-500-yr") into their cluster, and zones with
    shared nodes ("Charleston (regional)" and "(local)"). Grid sources are not grouped."""
    faults, zones = [], {}
    for name, (_, info) in models.items():
        if info.get("nshm_component") == "Fault" and "sources" in info:
            faults += [v["path"].split("/")[-1] for v in info["sources"].values()]
        if "zones" in info and isinstance(info["zones"], list):
            zones.update(_zone_labels(info["zones"]))
    acc = {}
    for name, source, typ, contribution, m, r, eps, _, _ in usgs["components"]["Total"]["sources"]:
        if source == "Grid":
            continue
        group = name
        if source == "FaultCluster":  # strip the recurrence branch, e.g. " : R1-500-yr"
            group = re.sub(r" : R\d+-[\d.]+-yr$", "", name)
        elif source == "Fault":
            for f in sorted(faults, key=len, reverse=True):
                if name.startswith(f) or name.startswith(f.rstrip(")") + ","):
                    group = f
                    break
        elif source == "Zone":
            group = zones.get(name, name)
        a = acc.setdefault((source, group), np.zeros(4))
        a += contribution * np.array([1.0, m, r, eps])
    return {k: [a[0], a[1] / a[0], a[2] / a[0], a[3] / a[0]] for k, a in acc.items()}


# ---------------------------------------------------------------------------------------
# ucla_plha
# ---------------------------------------------------------------------------------------


def nshmp_matrix(h, m, r, eps):
    """M x R x E sums of h binned as nshmp-lib bins (DisaggDataset): the bin index is
    (int) ((x - min) / delta), so a value on a bin edge, such as the magnitudes 6.6, 6.8, ...
    of the 0.2 magnitude bins starting at 4.4, often goes to the lower bin ((6.6 - 4.4) / 0.2 =
    10.999999999999996), and the outer epsilon bins are open."""
    n_m, n_r, n_e = len(M_EDGES) - 1, len(R_EDGES) - 1, len(EPS_EDGES) - 1
    mi = np.trunc((m - M_EDGES[0]) / 0.2).astype(int)
    ri = np.trunc((r - R_EDGES[0]) / 20.0).astype(int)
    ei = np.where(eps < -3.0, 0, np.where(eps >= 3.0, n_e - 1, np.trunc((eps + 3.0) / 0.5)))
    ei = ei.astype(int)
    ok = (mi >= 0) & (mi < n_m) & (ri >= 0) & (ri < n_r) & (h != 0)
    flat = (mi[ok] * n_r + ri[ok]) * n_e + ei[ok]
    return np.bincount(flat, weights=h[ok], minlength=n_m * n_r * n_e).reshape(n_m, n_r, n_e)


def _zone_labels(zones):
    """{zone name: group label}: nshm23_ceus_zone merges the nodes of zones with the same
    strike (e.g. Charleston regional and local), so their ruptures are grouped together."""
    by_strike = {}
    for z in zones:
        by_strike.setdefault(z.get("strike"), []).append(z["name"])
    return {n: " / ".join(v) for v in by_strike.values() for n in v}


class _Recorder:
    """Records the weighted per-rupture contributions that get_hazard averages (plha._disagg_sums)
    by contributor group and component, and with rjb and rrup.

    plha.get_source_data is wrapped to keep the rupture index of each source model, and
    plha._disagg_sums to accumulate the contributions (it is called for every source model and
    ground motion branch, for cluster models first with the independent ruptures and then with
    the cluster ruptures).
    """

    def __init__(self, models):
        self.models = models
        self.current = None
        self.sums = {}  # (component, group) -> [binned rate, m, rrup, rjb, eps] sums
        self.matrices = {}  # component -> M x R x E rates
        self.total = 0.0

    def groups(self, name, info, index):
        """Group names of the ruptures of a source model, and component of each part."""
        root = os.path.join(str(cun.files("ucla_plha").joinpath("source_models")),
                            self.models[name][0], name)
        rup = np.load(os.path.join(root, "ruptures.npz"))
        comp = info["nshm_component"]
        if info.get("cluster"):
            cid = np.asarray(rup["cluster_id"])[index]
            names = {c["cluster_id"]: c["name"] for c in json.load(open(os.path.join(root, "clusters.json")))} \
                if os.path.exists(os.path.join(root, "clusters.json")) else {
                    c["cluster_id"]: c.get("name", str(c["cluster_id"])) for c in info["clusters"]}
            indep = cid < 0
            sid = np.asarray(rup["source_id"])[index] if "source_id" in rup.files else np.zeros(len(index), int)
            g_indep = np.array([self._source_name(info, s) for s in sid[indep]], dtype=object)
            g_clus = np.array([names[c] for c in cid[~indep]], dtype=object)
            return [(comp, np.flatnonzero(indep), g_indep),
                    (info["cluster_nshm_component"], np.flatnonzero(~indep), g_clus)]
        if "source_id" in rup.files and "sources" in info:
            g = np.array([self._source_name(info, s) for s in np.asarray(rup["source_id"])[index]], dtype=object)
        elif os.path.exists(os.path.join(root, "strike.npy")) and "zones" in info:
            strike = np.load(os.path.join(root, "strike.npy"))[np.asarray(rup["node_index"])[index]]
            zones = json.load(open(os.path.join(root, "zones.json")))
            labels = _zone_labels(zones)
            label = {(-1.0 if z["strike"] is None else float(z["strike"])): labels[z["name"]] for z in zones}
            g = np.array([label[-1.0 if not np.isfinite(s) else float(s)] for s in strike], dtype=object)
        else:
            g = np.full(len(index), name, dtype=object)
        return [(comp, np.arange(len(index)), g)]

    @staticmethod
    def _source_name(info, source_id):
        s = info.get("sources", {}).get(str(int(source_id)))
        return s["path"].split("/")[-1] if s else str(int(source_id))

    def get_source_data(self, original):
        def wrapped(source_type, source_model, *args, extras=False, **kwargs):
            arrays, ex = original(source_type, source_model, *args, extras=True, **kwargs)
            info = kwargs.get("source_info") or plha.get_source_info(source_type, source_model)
            self.current = {
                "parts": self.groups(source_model, info, ex["index"]),
                "rjb": arrays[3],
                "rrup": arrays[4],
                "calls": 0,
            }
            return (arrays, ex) if extras else arrays

        return wrapped

    def disagg_sums(self, original):
        def wrapped(hazards, m, r, eps, in_bins):
            cur = self.current
            comp, sel, groups = cur["parts"][cur["calls"] % len(cur["parts"])]
            cur["calls"] += 1
            h = np.where(in_bins, hazards, 0.0)[0]  # one IML
            rrup = cur["rrup"][sel] if len(cur["rrup"]) else np.full(len(sel), np.nan)
            rjb = cur["rjb"][sel] if len(cur["rjb"]) else np.full(len(sel), np.nan)
            names, inv = np.unique(groups, return_inverse=True)
            for k, values in enumerate([h, h * m, h * rrup, h * rjb, h * eps[0]]):
                s = np.bincount(inv, weights=values, minlength=len(names))
                for g, x in zip(names, s):
                    self.sums.setdefault((comp, g), np.zeros(5))[k] += x
            d = nshmp_matrix(h, m, r, eps[0])
            self.matrices[comp] = self.matrices.get(comp, 0.0) + d
            return original(hazards, m, r, eps, in_bins)

        return wrapped


def run_ucla_plha(lon, lat, vs30, iml, models, distance_metric="rrup",
                  cluster_attribution="gmm", grid_distances="source_info", nshm_site_data=True):
    config = cun.nshm_config(lon, lat, vs30, [iml], models, "nshmp",
                             nshm_site_data=nshm_site_data)
    config["output"]["psha"]["disaggregation"] = {
        "magnitude_bin_edges": M_EDGES.tolist(),
        "distance_bin_edges": R_EDGES.tolist(),
        "epsilon_bin_edges": EPS_EDGES.tolist(),
        "distance_metric": distance_metric,
        "cluster_attribution": cluster_attribution,
        "means": True,
    }
    rec = _Recorder(models)
    originals = plha.get_source_data, plha._disagg_sums
    plha.get_source_data = rec.get_source_data(originals[0])
    plha._disagg_sums = rec.disagg_sums(originals[1])
    try:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "config.json")
            with open(path, "w") as f:
                json.dump(config, f)
            start = time.time()
            with warnings.catch_warnings(), cun.nshmp_grid_distances(grid_distances == "nshmp"):
                warnings.simplefilter("ignore")
                out = plha.get_hazard(path)
            seconds = time.time() - start
    finally:
        plha.get_source_data, plha._disagg_sums = originals
    psha = out["output"]["psha"]
    rate = psha["annual_rate_of_exceedance"][0]
    # the recorded contributions are those get_hazard bins and averages
    recorded = sum(s[0] for s in rec.sums.values())
    binned = np.sum(psha["disaggregation"]) / 100.0 * rate
    assert abs(recorded - binned) <= 1e-6 * rate, (recorded, binned)
    return {
        "rate": rate,
        "matrix": np.array(psha["disaggregation"])[0],
        "nshmp_matrix": 100.0 * sum(rec.matrices.values()) / rate,
        "means": psha["disaggregation_means"],
        "groups": rec.sums,
        "component_matrices": {k: 100.0 * v / rate for k, v in rec.matrices.items()},
        "seconds": seconds,
    }


# ---------------------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------------------


def _mean(s, rate):
    """[contribution %, m, rrup, rjb, eps] of recorded sums."""
    return [100.0 * s[0] / rate, s[1] / s[0], s[2] / s[0], s[3] / s[0], s[4] / s[0]]


def _distribution_stats(a, b):
    """Half L1 distances (% points) of the m, r, eps marginals and of the m-r-eps bins."""
    out = {}
    for name, axes in (("m", (1, 2)), ("r", (0, 2)), ("eps", (0, 1)), ("mr", (2,))):
        out[name] = 0.5 * float(np.abs(a.sum(axis=axes) - b.sum(axis=axes)).sum())
    out["mre"] = 0.5 * float(np.abs(a - b).sum())
    return out


def _mode(a):
    mr = a.sum(axis=2)
    i, j = np.unravel_index(np.argmax(mr), mr.shape)
    return float(0.5 * (M_EDGES[i] + M_EDGES[i + 1])), float(0.5 * (R_EDGES[j] + R_EDGES[j + 1])), float(mr[i, j])


def compare(name, lon, lat, vs30, models, cache, return_period=None, iml=None,
            distance_metric="rrup", cluster_attribution="gmm", grid_distances="source_info"):
    usgs = fetch_disagg(lon, lat, vs30, cache, return_period, iml)
    ucla = run_ucla_plha(lon, lat, vs30, usgs["iml"], models, distance_metric,
                         cluster_attribution, grid_distances)
    rate = ucla["rate"]
    comps = []
    for comp in ["Total"] + COMPONENTS:
        u = usgs["components"].get(comp)
        if comp == "Total":
            s = sum(ucla["groups"].values())
            mat = ucla["nshmp_matrix"]
        else:
            g = [v for (c, _), v in ucla["groups"].items() if c == comp]
            s = sum(g) if g else None
            mat = ucla["component_matrices"].get(comp)
        if u is None and (s is None or s[0] <= 0):
            continue
        row = {"component": comp}
        if u is not None:
            row["usgs"] = [u["binned"]] + u["mean"]
            row["usgs_mode"] = u["mode_mr"]
        if s is not None and s[0] > 0:
            row["ucla"] = _mean(s, rate)
            row["ucla_mode"] = _mode(mat)
        bins = usgs_total_bins(usgs) if comp == "Total" else (u or {}).get("bins")
        if u is not None and mat is not None and bins:
            row["distribution"] = _distribution_stats(usgs_matrix(bins), mat)
            if comp == "Total":  # ucla_plha bins (bin edges included in the upper bin)
                row["distribution_ucla_bins"] = _distribution_stats(
                    usgs_matrix(bins), ucla["matrix"])
        comps.append(row)
    ug = usgs_groups(usgs, models)
    groups = []
    keys = set(ug) | {k for k, v in ucla["groups"].items() if 100 * v[0] / rate >= 1.0}
    for key in keys:
        if key[0] == "Grid":
            continue
        row = {"component": key[0], "group": key[1]}
        if key in ug:
            row["usgs"] = ug[key]
        if key in ucla["groups"] and ucla["groups"][key][0] > 0:
            row["ucla"] = _mean(ucla["groups"][key], rate)
        groups.append(row)
    groups.sort(key=lambda r: -max(r.get("usgs", [0])[0], r.get("ucla", [0])[0]))
    return {
        "site": name, "longitude": lon, "latitude": lat, "vs30": vs30,
        "return_period": return_period, "iml": usgs["iml"],
        "usgs_rate": usgs["rate"], "usgs_recovered_rate": usgs["recovered_rate"],
        "ucla_rate": rate, "distance_metric": distance_metric,
        "cluster_attribution": cluster_attribution, "grid_distances": grid_distances,
        "seconds": ucla["seconds"], "components": comps, "groups": groups,
    }


def _fmt(values, n):
    if values is None:
        return "-".rjust(n)
    return values


def print_result(res):
    target = f"{res['return_period']:g} yr" if res["return_period"] else "IML"
    print(
        f"\n{res['site']} ({res['longitude']}, {res['latitude']}), Vs30 {res['vs30']:g}, {target}: "
        f"PGA {res['iml']:.4g} g, rate USGS {res['usgs_recovered_rate']:.4g} / ucla_plha {res['ucla_rate']:.4g} "
        f"({res['ucla_rate'] / res['usgs_recovered_rate'] - 1:+.1%}) [{res['seconds']:.0f} s; "
        f"{res['distance_metric']}, cluster attribution {res['cluster_attribution']}, "
        f"grid distances {res['grid_distances']}]"
    )
    print(f"  {'component':<14}{'% USGS/ucla':>16}{'mean M':>14}{'mean rRup (ucla rJB)':>26}{'mean eps':>16}"
          f"{'mode M, R (USGS bin mean/ucla bin)':>28}   dist. diff % (M, R, eps, MR, MRe)")
    for row in res["components"]:
        u = row.get("usgs") or [np.nan] * 4
        c = row.get("ucla") or [np.nan] * 5
        um = row.get("usgs_mode") or [np.nan] * 4
        cm = row.get("ucla_mode") or (np.nan, np.nan, np.nan)
        dist = row.get("distribution")
        d = "" if dist is None else ", ".join(f"{dist[k]:.1f}" for k in ("m", "r", "eps", "mr", "mre"))
        if "distribution_ucla_bins" in row:
            dist = row["distribution_ucla_bins"]
            d += " (ucla_plha bins: " + ", ".join(f"{dist[k]:.1f}" for k in ("m", "r", "eps", "mr", "mre")) + ")"
        print(f"  {row['component']:<14}{u[0]:7.2f} /{c[0]:6.2f}{u[1]:7.2f} /{c[1]:5.2f}"
              f"{u[2]:8.1f} /{c[2]:6.1f} ({c[3]:5.1f}){u[3]:8.2f} /{c[4]:5.2f}"
              f"{um[0]:9.2f},{um[1]:5.0f} /{cm[0]:5.2f},{cm[1]:4.0f}   {d}")
    print(f"  {'source (contributions >= 1%)':<56}{'% USGS/ucla':>14}{'mean M':>14}{'mean rRup':>16}{'mean eps':>15}")
    for row in res["groups"]:
        u = row.get("usgs") or [np.nan] * 4
        c = row.get("ucla") or [np.nan] * 5
        label = f"{row['component']}: {row['group']}"[:55]
        print(f"  {label:<56}{u[0]:7.2f} /{c[0]:6.2f}{u[1]:7.2f} /{c[1]:5.2f}"
              f"{u[2]:8.1f} /{c[2]:6.1f}{u[3]:8.2f} /{c[4]:5.2f}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--site", default="new_madrid")
    parser.add_argument("--lon", type=float)
    parser.add_argument("--lat", type=float)
    parser.add_argument("--vs30", type=float, default=760.0)
    parser.add_argument("--return-periods", type=float, nargs="*", default=[475.0, 2475.0])
    parser.add_argument("--iml", type=float, nargs="*", help="PGA (g) instead of return periods")
    parser.add_argument("--distance-metric", default="rrup", choices=["rrup", "rjb", "default"])
    parser.add_argument("--cluster-attribution", default="gmm", choices=["gmm", "branch"])
    parser.add_argument("--grid-distances", default="source_info", choices=["source_info", "nshmp"],
                        help="see compare_usgs_nshm.py (affects nshm23_wus_grid only)")
    parser.add_argument("--models", nargs="*")
    parser.add_argument("--cache", default=os.path.join(os.path.dirname(__file__), "usgs_nshm_cache"))
    parser.add_argument("--out", help="write the results to this JSON file")
    args = parser.parse_args(argv)

    if args.lon is not None and args.lat is not None:
        sites = {"custom": (args.lon, args.lat)}
    elif args.site == "all":
        sites = cun.SITES
    elif args.site == "ceus":
        sites = {k: cun.SITES[k] for k in CEUS_SITES}
    else:
        sites = {args.site: cun.SITES[args.site]}
    models = cun.available_models()
    if args.models:
        models = {k: v for k, v in models.items() if k in args.models}
    targets = [(None, x) for x in args.iml] if args.iml else [(rp, None) for rp in args.return_periods]
    results = []
    for name, (lon, lat) in sites.items():
        for rp, iml in targets:
            res = compare(name, lon, lat, args.vs30, models, args.cache, rp, iml,
                          args.distance_metric, args.cluster_attribution, args.grid_distances)
            print_result(res)
            results.append(res)
    if args.out:
        with open(args.out, "w") as f:
            json.dump(results, f, indent=1, default=float)
    return results


if __name__ == "__main__":
    sys.exit(0 if main() else 1)

"""Hazard curves at many sites (e.g. on a cluster node), with the source models loaded once.

The sites are given in a table (CSV, or Parquet if pandas is installed) with the columns

* ``site_id``: name of the site (unique)
* ``latitude``, ``longitude``, ``vs30``: as in the config "site"
* optional: ``elevation`` (default 0), ``measured_vs30``, ``z1p0``, ``z2p5``
* optional: ``profiles_file``: JSON file (path relative to the table) with the soil profiles of
  the site, a list of {"name": ..., "liquefaction_models": {...}} as the config "profiles"

and a base config file gives everything else (source models, ground motion models,
constraints, output "psha" and/or "plha" values, and "profiles" or "liquefaction_models" for
the sites without a profiles_file; "liquefaction_models" is one profile named "0").
Disaggregation and output files of the base config are not used. The curves are exactly those
of plha.get_hazard for each site and profile.

The source models are read and prepared once (plha.load_source_models) and the sites are run
by a pool of worker processes (fork on Linux, so the workers share the source model arrays;
spawn elsewhere, where each worker loads the source models). The sites are split into chunks
of --chunk-size sites; each chunk is written by its worker to one file in --out-dir,
chunk_<first site index>.npz (or .h5 with --format hdf5 and h5py), with

* ``site_id``: site ids of the chunk (strings)
* ``pga``, ``fsl``: the ground motion and factor of safety values
* ``psha``: annual rate of exceedance, shape = sites x pga (if the config has output psha)
* ``plha``: annual rate of non-exceedance, shape = sites x profiles x fsl (if output plha)
* ``profile``: names of the profiles (sites x profiles; "" for missing profiles of sites
  with fewer profiles, whose plha values are NaN)

Chunks whose file exists are skipped, so an interrupted run can be restarted.

Usage::

    python -m ucla_plha.batch --sites sites.csv --config base.json --out-dir results \\
        --processes 128 [--start 0 --max-sites 100000] [--chunk-size 100]

On a cluster node (e.g. TACC Stampede3 or Vista), run one process per core with
--processes equal to the number of cores, and split a large table over the nodes of a SLURM
job array with --start and --max-sites (e.g. --start $((SLURM_ARRAY_TASK_ID * 100000))
--max-sites 100000), so that each array task writes its own chunk files. Set
OMP_NUM_THREADS=1 (numpy then uses one thread per process). For example, a SLURM job array
script (one node per array task)::

    #SBATCH -N 1 -n 1 --array=0-9
    export OMP_NUM_THREADS=1
    python -m ucla_plha.batch --sites sites.csv --config base.json --out-dir results \\
        --processes $(nproc) --start $((SLURM_ARRAY_TASK_ID * 100000)) --max-sites 100000

mpi4py could distribute the chunks over several nodes in one job instead, but is not needed.
read_results(out_dir) reads the chunk files.
"""

import argparse
import csv
import json
import multiprocessing
import os
import sys
import time
import warnings

import numpy as np

from ucla_plha import plha

SITE_COLUMNS = {
    "latitude": float,
    "longitude": float,
    "vs30": float,
    "elevation": float,
    "z1p0": float,
    "z2p5": float,
}


def read_sites(path):
    """Rows of the sites table (list of dicts with the values as strings or numbers)."""
    if path.lower().endswith((".parquet", ".pq")):
        import pandas

        return pandas.read_parquet(path).to_dict("records")
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _is_missing(value):
    return value is None or (isinstance(value, str) and value.strip() == "") or (
        isinstance(value, float) and np.isnan(value)
    )


def site_config(base, row, table_dir, profile_cache=None):
    """Config of the site of a row of the sites table (see the module docstring)."""
    config = {k: v for k, v in base.items() if k not in ("site", "profiles", "liquefaction_models")}
    site = {k: v for k, v in base.get("site", {}).items() if k not in ("latitude", "longitude")}
    site.setdefault("elevation", 0.0)
    for key, kind in SITE_COLUMNS.items():
        if key in row and not _is_missing(row[key]):
            site[key] = kind(row[key])
    if "measured_vs30" in row and not _is_missing(row["measured_vs30"]):
        value = row["measured_vs30"]
        site["measured_vs30"] = value if isinstance(value, bool) else str(value).lower() in ("1", "true")
    config["site"] = site
    output = {}
    for kind, values in (("psha", "pga"), ("plha", "fsl")):
        if kind in base["output"]:
            output[kind] = {values: base["output"][kind][values]}
    config["output"] = output
    profiles_file = row.get("profiles_file")
    if not _is_missing(profiles_file):
        path = os.path.join(table_dir, str(profiles_file))
        if profile_cache is None or path not in profile_cache:
            with open(path) as f:
                profiles = json.load(f)
            if profile_cache is not None:
                profile_cache[path] = profiles
        else:
            profiles = profile_cache[path]
        config["profiles"] = profiles
    elif "profiles" in base:
        config["profiles"] = base["profiles"]
    elif "liquefaction_models" in base:
        config["profiles"] = [{"name": "0", "liquefaction_models": base["liquefaction_models"]}]
    return config


def run_sites(base, rows, table_dir, models):
    """Hazard curves of the sites (rows of the table). Returns a dict of arrays (see module)."""
    profile_cache = {}
    psha, plha_curves, names = [], [], []
    pga = np.asarray(base["output"].get("psha", {}).get("pga", []), dtype=float)
    fsl = np.asarray(base["output"].get("plha", {}).get("fsl", []), dtype=float)
    for row in rows:
        config = site_config(base, row, table_dir, profile_cache)
        out = plha.get_hazard(config, models=models)
        if out is None:
            raise ValueError(f"site {row['site_id']}: config error (see the message above)")
        if "psha" in out["output"]:
            psha.append(out["output"]["psha"]["annual_rate_of_exceedance"])
        profiles = out["output"].get("profiles", [])
        plha_curves.append([p["plha"]["annual_rate_of_nonexceedance"] for p in profiles])
        names.append([p["name"] for p in profiles])
    result = {"site_id": np.array([str(r["site_id"]) for r in rows]), "pga": pga, "fsl": fsl}
    if "psha" in base["output"]:
        result["psha"] = np.array(psha, dtype=float).reshape(len(rows), len(pga))
    if "plha" in base["output"]:
        n_profiles = max([len(n) for n in names] + [0])
        curves = np.full((len(rows), n_profiles, len(fsl)), np.nan)
        profile = np.full((len(rows), n_profiles), "", dtype=object)
        for i, (c, n) in enumerate(zip(plha_curves, names)):
            if c:
                curves[i, : len(c)] = c
                profile[i, : len(n)] = n
        result["plha"] = curves
        result["profile"] = profile.astype(str)
    return result


def write_result(filename, result):
    """Writes the arrays of a chunk to an npz or (.h5) HDF5 file, through a temporary file."""
    tmp = filename + ".tmp"
    if filename.endswith(".h5"):
        import h5py

        with h5py.File(tmp, "w") as f:
            for key, value in result.items():
                if value.dtype.kind == "U":
                    value = value.astype(h5py.string_dtype())
                f.create_dataset(key, data=value)
    else:
        with open(tmp, "wb") as f:
            np.savez(f, **result)
    os.replace(tmp, filename)


def read_results(out_dir):
    """{site_id: {"psha": curve, "plha": {profile: curve}}} of all chunk files in out_dir."""
    results = {}
    for name in sorted(os.listdir(out_dir)):
        path = os.path.join(out_dir, name)
        if name.startswith("chunk_") and name.endswith(".npz"):
            data = dict(np.load(path))
        elif name.startswith("chunk_") and name.endswith(".h5"):
            import h5py

            with h5py.File(path, "r") as f:
                data = {k: f[k][()] for k in f}
            for k in ("site_id", "profile"):
                if k in data:
                    data[k] = data[k].astype(str)
        else:
            continue
        for i, site_id in enumerate(data["site_id"]):
            entry = {}
            if "psha" in data:
                entry["psha"] = data["psha"][i]
            if "plha" in data:
                entry["plha"] = {
                    p: data["plha"][i, j] for j, p in enumerate(data["profile"][i]) if p != ""
                }
            results[str(site_id)] = entry
    return results


_WORKER = {}


def _init_worker(base, table_dir, models, warmup=None):
    """Worker set-up: the source models of the parent (fork), or loaded here for the warm-up
    config (spawn)."""
    warnings.simplefilter("ignore")
    if models is None:
        models = plha.load_source_models(warmup)
    _WORKER.update(base=base, table_dir=table_dir, models=models)


def _run_chunk(task):
    filename, rows = task
    start = time.perf_counter()
    result = run_sites(_WORKER["base"], rows, _WORKER["table_dir"], _WORKER["models"])
    write_result(filename, result)
    return filename, len(rows), time.perf_counter() - start


def _warmup_config(base, row, table_dir):
    config = site_config(base, row, table_dir)
    config.pop("profiles", None)
    return config


def run(sites, config, out_dir, processes=1, start=0, max_sites=None, chunk_size=100,
        file_format="npz", start_method=None, verbose=True):
    """Runs the sites of a table (see the module docstring). Returns the chunk files written."""
    rows = read_sites(sites)
    rows = rows[start:] if max_sites is None else rows[start : start + max_sites]
    table_dir = os.path.dirname(os.path.abspath(sites))
    with open(config) as f:
        base = json.loads(f.read().lower())
    for key in ("outputfile",):
        base.get("output", {}).pop(key, None)
    os.makedirs(out_dir, exist_ok=True)
    extension = ".h5" if file_format == "hdf5" else ".npz"
    tasks = []
    for i in range(0, len(rows), chunk_size):
        filename = os.path.join(out_dir, f"chunk_{start + i:09d}{extension}")
        if not os.path.exists(filename):
            tasks.append((filename, rows[i : i + chunk_size]))
    if not tasks:
        return []
    if start_method is None:
        # fork shares the source models of this process with the workers (Linux); spawn
        # elsewhere (fork is not available on Windows and not safe on macOS)
        start_method = "fork" if sys.platform.startswith("linux") else "spawn"
    processes = min(processes, len(tasks))
    t0 = time.perf_counter()
    models = None
    if processes <= 1 or start_method == "fork":
        # loaded once in this process (and shared by the forked workers)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            models = plha.load_source_models(_warmup_config(base, tasks[0][1][0], table_dir))
    if verbose:
        print(f"{len(rows)} sites in {len(tasks)} chunks; source models loaded in "
              f"{time.perf_counter() - t0:.1f} s", flush=True)
    written = []
    if processes <= 1:
        _init_worker(base, table_dir, models)
        results = map(_run_chunk, tasks)
        for filename, n, seconds in results:
            written.append(filename)
            if verbose:
                print(f"{filename}: {n} sites, {seconds / n:.2f} s per site", flush=True)
    else:
        warmup = None if models is not None else _warmup_config(base, tasks[0][1][0], table_dir)
        ctx = multiprocessing.get_context(start_method)
        initargs = (base, table_dir, models, warmup)
        with ctx.Pool(processes, initializer=_init_worker, initargs=initargs) as pool:
            for filename, n, seconds in pool.imap_unordered(_run_chunk, tasks):
                written.append(filename)
                if verbose:
                    print(f"{filename}: {n} sites, {seconds / n:.2f} s per site", flush=True)
    if verbose:
        print(f"done in {time.perf_counter() - t0:.1f} s", flush=True)
    return sorted(written)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sites", required=True, help="sites table (CSV or Parquet)")
    parser.add_argument("--config", required=True, help="base config file")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--processes", type=int, default=os.cpu_count())
    parser.add_argument("--start", type=int, default=0, help="first site (row) of the table")
    parser.add_argument("--max-sites", type=int, default=None, help="number of sites to run")
    parser.add_argument("--chunk-size", type=int, default=100, help="sites per output file")
    parser.add_argument("--format", default="npz", choices=["npz", "hdf5"])
    parser.add_argument("--start-method", default=None, choices=["fork", "spawn", "forkserver"])
    args = parser.parse_args(argv)
    run(args.sites, args.config, args.out_dir, args.processes, args.start, args.max_sites,
        args.chunk_size, args.format, args.start_method)


if __name__ == "__main__":
    main()

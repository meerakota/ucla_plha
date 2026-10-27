"""Builds the on-disk cache of the site-independent GC2 set-up of the fault source models
(ucla_plha.gc2_cache), so that later processes read it instead of computing it.

Run once after installing ucla_plha (e.g. on a cluster, before a job array):

    python -m ucla_plha.precompute [--cache-dir DIR] [--force] [model ...]

The models are the fault source models of the package (default: all of them). --cache-dir sets
the cache directory (default: UCLA_PLHA_CACHE_DIR or the user cache directory, see
ucla_plha.gc2_cache); the jobs must then use the same UCLA_PLHA_CACHE_DIR. --force recomputes
files that are already in the cache.
"""

import argparse
import os
import sys
import time

import numpy as np

from ucla_plha import gc2_cache, plha


def fault_source_models():
    """Names of the fault source models of the package."""
    root = plha._source_model_path("fault_source_models", "x").parent
    return sorted(
        name
        for name in os.listdir(str(root))
        if os.path.exists(str(root.joinpath(name).joinpath("rect_rjb.npy")))
    )


def precompute(model, directory=None, force=False):
    """Computes the GC2 set-up of a fault source model and writes it to the cache. Returns
    (file name or None, seconds, whether it was computed)."""
    path = plha._source_model_path("fault_source_models", model)
    name = "ruptures_segments.npy.npz"
    if not os.path.exists(str(path.joinpath(name))):
        name = "ruptures_segments.npz"
    with np.load(str(path.joinpath(name))) as rs:
        segment_index = rs["segment_index"]
        rupture_index = rs["rupture_index"]
    boundaries = np.r_[0, np.where(np.diff(rupture_index) != 0)[0] + 1]
    local = plha._LocalGeometry(path)
    model_key = gc2_cache.key(path)
    expected = {
        "seg_length": len(local["piece_first"]),
        "piece_before": len(local["top_lat"]),
        "flip": len(segment_index),
        "offset": len(segment_index),
        "u_max": len(boundaries),
    }
    start = time.perf_counter()
    if not force and gc2_cache.load(model, model_key, expected, directory) is not None:
        return gc2_cache.file_name(model, model_key, directory or gc2_cache.cache_dir()), 0.0, False
    setup = plha._compute_gc2_setup(local, segment_index, boundaries)
    out = gc2_cache.save(model, model_key, setup, directory)
    return out, time.perf_counter() - start, True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("models", nargs="*", help="fault source models (default: all)")
    parser.add_argument("--cache-dir", default=None, help="cache directory")
    parser.add_argument("--force", action="store_true", help="recompute cached files")
    args = parser.parse_args(argv)
    directory = args.cache_dir or gc2_cache.cache_dir()
    if directory is None:
        print("the cache is disabled (UCLA_PLHA_CACHE_DIR is empty or none)", file=sys.stderr)
        return 1
    models = args.models or fault_source_models()
    for model in models:
        out, seconds, computed = precompute(model, directory, args.force)
        if out is None:
            print(f"{model}: could not write to {directory}", file=sys.stderr)
            return 1
        size = os.path.getsize(out) / 1e6
        state = f"computed in {seconds:.1f} s" if computed else "already cached"
        print(f"{model}: {out} ({size:.1f} MB), {state}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

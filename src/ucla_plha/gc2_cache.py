"""On-disk cache of the site-independent GC2 set-up of the fault source models.

The GC2 set-up of a fault source model (plha._gc2_setup: geometry.gc2_rupture_setup of all its
ruptures, and the segment and piece lengths) does not depend on the site, but it takes several
seconds per model (UCERF3 FM3.1 about 6 s, FM3.2 7 s, NSHM23 WUS 9 s) and is too large to ship
with the package (75 to 140 MB compressed per model, against 9 to 19 MB of source model data).
It is therefore computed once and written to a cache directory, from which later processes read
it (about 0.3 s per model):

    UCLA_PLHA_CACHE_DIR (environment variable): the cache directory; "" or "none" disables the
        cache. Default: $XDG_CACHE_HOME/ucla_plha (~/.cache/ucla_plha) on Linux and macOS,
        %LOCALAPPDATA%/ucla_plha/Cache on Windows.

python -m ucla_plha.precompute builds the cache ahead of a cluster run. Files are named after a
key, the SHA-256 of the format version, the numpy version, GC2_RUN_GAP, and the contents of the
input geometry files (rect_rjb.npy, rect_segment_id.npy, ruptures_segments.npz), so a changed
model or version uses a new file. Each file is an uncompressed npz (the CRC-32 of every member is
checked when it is read) that also holds the key; it is written to a temporary file in the cache
directory and renamed, so concurrent processes never read a partial file. A missing, stale, or
corrupt file, or an unusable cache directory, falls back to computing the set-up; the results are
bit-identical with and without the cache.
"""

import hashlib
import os
import tempfile
import time
import zipfile

import numpy as np

from ucla_plha.geometry import geometry

FORMAT_VERSION = 1
INPUT_FILES = (
    "rect_rjb.npy", "rect_segment_id.npy", "ruptures_segments.npz", "ruptures_segments.npy.npz",
)
ARRAYS = ("seg_length", "piece_before", "flip", "offset", "u_max")
_DTYPES = {"seg_length": "f8", "piece_before": "f8", "flip": "?", "offset": "f8", "u_max": "f8"}
# attempts to read a file that another process is replacing (PermissionError on Windows)
_READ_ATTEMPTS = 10


def cache_dir():
    """The cache directory (UCLA_PLHA_CACHE_DIR, or the default user cache directory), or None if
    the cache is disabled."""
    value = os.environ.get("UCLA_PLHA_CACHE_DIR")
    if value is not None:
        if value.strip().lower() in ("", "none"):
            return None
        return value
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), "AppData", "Local")
        return os.path.join(base, "ucla_plha", "Cache")
    base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    return os.path.join(base, "ucla_plha")


def key(path):
    """Cache key of the GC2 set-up of the fault source model in directory path (hex SHA-256)."""
    h = hashlib.sha256()
    h.update(f"gc2 format {FORMAT_VERSION}; numpy {np.__version__}; run_gap {geometry.GC2_RUN_GAP!r}".encode())
    for name in INPUT_FILES:
        source = str(path.joinpath(name)) if hasattr(path, "joinpath") else os.path.join(path, name)
        h.update(b"\0" + name.encode() + b"\0")
        if os.path.exists(source):
            with open(source, "rb") as f:
                for block in iter(lambda: f.read(1 << 20), b""):
                    h.update(block)
        else:
            h.update(b"missing")
    return h.hexdigest()


def file_name(model, model_key, directory):
    """Path of the cache file of a model with the given key."""
    return os.path.join(directory, f"gc2_{model}_{model_key[:24]}.npz")


def load(model, model_key, expected, directory=None):
    """The cached GC2 set-up of a model, or None if it is not in the cache or not valid.

    expected: {array name: length} of the set-up (number of segments, pieces, rupture entries,
    and ruptures); a file whose arrays do not have these lengths and dtypes, or whose key differs,
    is not used.
    """
    directory = cache_dir() if directory is None else directory
    if directory is None:
        return None
    name = file_name(model, model_key, directory)
    if not os.path.exists(name):
        return None
    for attempt in range(_READ_ATTEMPTS):
        try:
            with np.load(name, allow_pickle=False) as data:
                if str(data["key"][()]) != model_key:
                    return None
                setup = {k: np.array(data[k]) for k in ARRAYS}
            break
        except PermissionError:
            # Windows: the file is being replaced by another process (os.replace); try again
            if attempt == _READ_ATTEMPTS - 1:
                return None
            time.sleep(0.02 * (attempt + 1))
        except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile):
            return None
    for k, v in setup.items():
        if v.ndim != 1 or len(v) != expected[k] or v.dtype != np.dtype(_DTYPES[k]):
            return None
    return setup


def save(model, model_key, setup, directory=None):
    """Writes the GC2 set-up of a model to the cache (atomically: a temporary file in the cache
    directory, then renamed). Returns the file name, or None if the cache is disabled or the
    directory cannot be written (the set-up is then computed again by the next process)."""
    directory = cache_dir() if directory is None else directory
    if directory is None:
        return None
    name = file_name(model, model_key, directory)
    tmp = None
    try:
        os.makedirs(directory, exist_ok=True)
        fd, tmp = tempfile.mkstemp(prefix=os.path.basename(name) + ".", suffix=".tmp", dir=directory)
        with os.fdopen(fd, "wb") as f:
            np.savez(f, key=np.array(model_key), **{k: setup[k] for k in ARRAYS})
        os.replace(tmp, name)
        tmp = None
        return name
    except OSError:
        return None
    finally:
        if tmp is not None:
            try:
                os.remove(tmp)
            except OSError:
                pass

"""On-disk cache of the site-independent GC2 set-up (ucla_plha.gc2_cache, ucla_plha.precompute)."""

import copy
import os
import threading

import numpy as np
import pytest

from ucla_plha import gc2_cache, plha, precompute
from pga_grid_helpers import CONFIG

# A small fault source model (GC2 set-up in milliseconds): the ruptures of UCERF3 FM3.1 on its
# first SMALL_SEGMENTS segments, written to a temporary directory (see small_model)
SMALL = "small_fault"
SMALL_SEGMENTS = 60


@pytest.fixture(scope="session")
def small_model_dir(tmp_path_factory):
    source = plha._source_model_path("fault_source_models", "ucerf3_fm31")
    name = "ruptures_segments.npy.npz"
    if not os.path.exists(str(source.joinpath(name))):
        name = "ruptures_segments.npz"
    with np.load(str(source.joinpath(name))) as rs:
        segment_index = rs["segment_index"]
        rupture_index = rs["rupture_index"]
    rect = np.load(str(source.joinpath("rect_rjb.npy")))
    rect_segment_id = np.load(str(source.joinpath("rect_segment_id.npy")))
    # ruptures whose segments are all among the first SMALL_SEGMENTS
    outside = np.zeros(int(rupture_index.max()) + 1, dtype=bool)
    outside[rupture_index[segment_index >= SMALL_SEGMENTS]] = True
    keep = ~outside[rupture_index]
    _, new_rupture = np.unique(rupture_index[keep], return_inverse=True)
    directory = tmp_path_factory.mktemp("models") / SMALL
    directory.mkdir()
    np.savez(
        directory / "ruptures_segments.npz",
        segment_index=segment_index[keep],
        rupture_index=new_rupture.astype(rupture_index.dtype),
    )
    in_small = rect_segment_id < SMALL_SEGMENTS
    np.save(directory / "rect_rjb.npy", rect[in_small])
    np.save(directory / "rect_segment_id.npy", rect_segment_id[in_small])
    assert len(np.unique(new_rupture)) > 10
    return directory


@pytest.fixture(autouse=True)
def small_model(small_model_dir, monkeypatch):
    """plha._source_model_path gives the temporary directory for the small model."""
    original = plha._source_model_path

    def source_model_path(source_type, source_model):
        if source_type == "fault_source_models" and source_model == SMALL:
            return small_model_dir
        return original(source_type, source_model)

    monkeypatch.setattr(plha, "_source_model_path", source_model_path)


def _inputs(model):
    path = plha._source_model_path("fault_source_models", model)
    with np.load(str(path.joinpath("ruptures_segments.npz"))) as rs:
        segment_index = rs["segment_index"]
        rupture_index = rs["rupture_index"]
    boundaries = np.r_[0, np.where(np.diff(rupture_index) != 0)[0] + 1]
    return path, segment_index, boundaries


def _setup(model, monkeypatch):
    """GC2 set-up through plha._gc2_setup in a fresh session (no in-process copy)."""
    monkeypatch.setattr(plha, "_LOCAL_GEOMETRY", {})
    path, segment_index, boundaries = _inputs(model)
    return plha._gc2_setup(plha._local_geometry(path), segment_index, boundaries)


def _computed(model):
    path, segment_index, boundaries = _inputs(model)
    return plha._compute_gc2_setup(plha._LocalGeometry(path), segment_index, boundaries)


def _equal(a, b):
    assert set(a) == set(b) == set(gc2_cache.ARRAYS)
    for k in a:
        assert a[k].dtype == b[k].dtype and np.array_equal(a[k], b[k]), k


@pytest.fixture
def cache(tmp_path, monkeypatch):
    directory = tmp_path / "cache"
    monkeypatch.setenv("UCLA_PLHA_CACHE_DIR", str(directory))
    return directory


def _file(model, directory):
    path = plha._source_model_path("fault_source_models", model)
    return gc2_cache.file_name(model, gc2_cache.key(path), str(directory))


def test_miss_then_hit(cache, monkeypatch):
    expected = _computed(SMALL)
    first = _setup(SMALL, monkeypatch)  # miss: computed and written
    assert os.path.exists(_file(SMALL, cache))
    _equal(first, expected)
    calls = []
    monkeypatch.setattr(plha, "_compute_gc2_setup", lambda *a: calls.append(1))
    second = _setup(SMALL, monkeypatch)  # hit: read, not computed
    assert calls == []
    _equal(second, expected)
    assert [p.name for p in cache.iterdir()] == [os.path.basename(_file(SMALL, cache))]


def test_disabled_and_unwritable(tmp_path, monkeypatch):
    expected = _computed(SMALL)
    for value in ("", "none", "None"):
        monkeypatch.setenv("UCLA_PLHA_CACHE_DIR", value)
        assert gc2_cache.cache_dir() is None
        _equal(_setup(SMALL, monkeypatch), expected)
    blocker = tmp_path / "file"
    blocker.write_text("not a directory")
    monkeypatch.setenv("UCLA_PLHA_CACHE_DIR", str(blocker / "cache"))
    _equal(_setup(SMALL, monkeypatch), expected)
    assert gc2_cache.save(SMALL, "0" * 64, expected) is None


def test_default_directory(monkeypatch):
    monkeypatch.delenv("UCLA_PLHA_CACHE_DIR", raising=False)
    assert "ucla_plha" in gc2_cache.cache_dir()


def test_stale_key(cache, monkeypatch):
    expected = _computed(SMALL)
    path = plha._source_model_path("fault_source_models", SMALL)
    model_key = gc2_cache.key(path)
    # a file under this name with another key inside (e.g. a hash prefix collision) is not used
    os.makedirs(cache)
    wrong = {k: v.copy() for k, v in expected.items()}
    wrong["offset"] = wrong["offset"] + 1.0
    gc2_cache.save(SMALL, model_key[:24] + "f" * 40, wrong)
    os.replace(gc2_cache.file_name(SMALL, model_key[:24] + "f" * 40, str(cache)), _file(SMALL, cache))
    sizes = {k: len(v) for k, v in expected.items()}
    assert gc2_cache.load(SMALL, model_key, sizes) is None
    _equal(_setup(SMALL, monkeypatch), expected)  # recomputed and rewritten
    _equal(gc2_cache.load(SMALL, model_key, sizes), expected)
    # a new format version or changed inputs give a new key (and file)
    monkeypatch.setattr(gc2_cache, "FORMAT_VERSION", gc2_cache.FORMAT_VERSION + 1)
    assert gc2_cache.key(path) != model_key
    # arrays of other lengths are not used
    assert gc2_cache.load(SMALL, model_key, dict(sizes, u_max=sizes["u_max"] + 1)) is None


def test_key_depends_on_inputs(tmp_path):
    path = plha._source_model_path("fault_source_models", SMALL)
    copy_dir = tmp_path / SMALL
    copy_dir.mkdir()
    for name in gc2_cache.INPUT_FILES:
        source = path.joinpath(name)
        if os.path.exists(str(source)):
            (copy_dir / name).write_bytes(source.read_bytes())
    assert gc2_cache.key(str(copy_dir)) == gc2_cache.key(path)
    data = bytearray((copy_dir / "rect_rjb.npy").read_bytes())
    data[-1] ^= 1
    (copy_dir / "rect_rjb.npy").write_bytes(bytes(data))
    assert gc2_cache.key(str(copy_dir)) != gc2_cache.key(path)


@pytest.mark.parametrize("damage", ["truncate", "flip_byte", "garbage", "empty"])
def test_corrupt_file(cache, monkeypatch, damage):
    expected = _computed(SMALL)
    _setup(SMALL, monkeypatch)
    name = _file(SMALL, cache)
    data = bytearray(open(name, "rb").read())
    if damage == "truncate":
        data = data[: len(data) // 2]
    elif damage == "flip_byte":
        data[len(data) // 2] ^= 0xFF  # inside an array: CRC-32 of the npz member
    elif damage == "garbage":
        data = bytearray(b"not an npz file" * 10)
    else:
        data = bytearray()
    with open(name, "wb") as f:
        f.write(bytes(data))
    path = plha._source_model_path("fault_source_models", SMALL)
    sizes = {k: len(v) for k, v in expected.items()}
    assert gc2_cache.load(SMALL, gc2_cache.key(path), sizes) is None
    _equal(_setup(SMALL, monkeypatch), expected)  # computed again (and rewritten)
    _equal(gc2_cache.load(SMALL, gc2_cache.key(path), sizes), expected)


def test_atomic_write(cache, monkeypatch):
    expected = _computed(SMALL)
    path = plha._source_model_path("fault_source_models", SMALL)
    model_key = gc2_cache.key(path)

    # a failure while writing leaves no file (and no temporary file) behind
    def failing_savez(f, **arrays):
        f.write(b"partial")
        raise OSError("disk full")

    savez = gc2_cache.np.savez
    monkeypatch.setattr(gc2_cache.np, "savez", failing_savez)
    assert gc2_cache.save(SMALL, model_key, expected) is None
    assert list(cache.iterdir()) == []
    monkeypatch.setattr(gc2_cache.np, "savez", savez)

    # the file appears only by the rename of a complete temporary file
    seen = []
    replace = os.replace

    def checking_replace(src, dst):
        seen.append((os.path.dirname(src), os.path.exists(dst)))
        with np.load(src) as data:
            assert str(data["key"][()]) == model_key
        replace(src, dst)

    monkeypatch.setattr(gc2_cache.os, "replace", checking_replace)
    assert gc2_cache.save(SMALL, model_key, expected) == _file(SMALL, cache)
    assert seen == [(str(cache), False)]
    assert [p.name for p in cache.iterdir()] == [os.path.basename(_file(SMALL, cache))]


def test_concurrent_writers(cache):
    expected = _computed(SMALL)
    path = plha._source_model_path("fault_source_models", SMALL)
    model_key = gc2_cache.key(path)
    sizes = {k: len(v) for k, v in expected.items()}
    results = []

    def write_and_read():
        out = gc2_cache.save(SMALL, model_key, expected)
        got = gc2_cache.load(SMALL, model_key, sizes)
        results.append((out, got))

    threads = [threading.Thread(target=write_and_read) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    # every reader got a complete file (a writer may fail to replace a file that is open on
    # Windows; it then leaves the existing file)
    assert any(out is not None for out, _ in results)
    for _, got in results:
        assert got is not None
        _equal(got, expected)
    assert [p.name for p in cache.iterdir()] == [os.path.basename(_file(SMALL, cache))]


def test_precompute_command(cache, capsys):
    assert precompute.main([SMALL]) == 0
    assert "computed" in capsys.readouterr().out
    assert precompute.main([SMALL]) == 0
    assert "already cached" in capsys.readouterr().out
    assert precompute.main([SMALL, "--force"]) == 0
    assert "computed" in capsys.readouterr().out
    assert precompute.fault_source_models() == ["ucerf3_fm31", "ucerf3_fm32"]


def test_hazard_identical_with_and_without_cache(tmp_path, monkeypatch):
    # UCERF3 FM3.1 with the default GC2 Rx and Ry0: cache disabled, cache miss, and cache hit
    config = copy.deepcopy(CONFIG)
    config["output"] = copy.deepcopy(CONFIG["output"])
    config["output"]["psha"] = {"pga": [0.05, 0.1, 0.3, 0.6]}
    config["ground_motion_models"] = {"ask14": {"weight": 1.0}}  # uses Rx and Ry0
    curves = []
    for value in ("none", str(tmp_path / "cache"), str(tmp_path / "cache")):
        monkeypatch.setenv("UCLA_PLHA_CACHE_DIR", value)
        monkeypatch.setattr(plha, "_LOCAL_GEOMETRY", {})
        out = plha.get_hazard(copy.deepcopy(config))
        curves.append(
            (out["output"]["psha"]["annual_rate_of_exceedance"], out["output"]["plha"]["annual_rate_of_nonexceedance"])
        )
    assert len(list((tmp_path / "cache").iterdir())) == 1
    assert curves[0] == curves[1] == curves[2]

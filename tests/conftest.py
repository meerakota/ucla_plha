import os

import pytest


@pytest.fixture(autouse=True, scope="session")
def _gc2_cache_dir(tmp_path_factory):
    """The tests write the GC2 cache (ucla_plha.gc2_cache) to a temporary directory, not to the
    user cache directory (tests that need another directory set UCLA_PLHA_CACHE_DIR themselves)."""
    previous = os.environ.get("UCLA_PLHA_CACHE_DIR")
    os.environ["UCLA_PLHA_CACHE_DIR"] = str(tmp_path_factory.mktemp("gc2_cache"))
    yield
    if previous is None:
        os.environ.pop("UCLA_PLHA_CACHE_DIR", None)
    else:
        os.environ["UCLA_PLHA_CACHE_DIR"] = previous

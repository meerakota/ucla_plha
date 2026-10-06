"""
Tests of the location-dependent USGS NSHM site data (coastal plain sediment thickness zSed and
the Coastal Plain CPA region ground motion model tree, nshm-conus 6.2.0 site-data/).
"""

import pytest

from ucla_plha import nshm_site_data, plha

# (lon, lat): zSed (km) of the nshm-conus coastal-plain.csv at the snapped 0.05 degree grid
# point (nshmp-lib SiteData), or None off the coastal plain
SITES = {
    "new_madrid": ((-89.6, 36.6), 0.903, True),  # 902.5 m, rounded HALF_UP to 0.903 km
    "memphis": ((-90.05, 35.15), 1.222, True),
    "charleston": ((-79.94, 32.78), 0.854, True),
    # inside the CPA region polygon, but outside the margin data (zSed NaN in nshmp-lib)
    "central_virginia": ((-77.9, 37.8), None, True),
    "quiet_ceus": ((-99.0, 41.0), None, False),
    "los_angeles": ((-118.25, 34.05), None, False),
}


@pytest.mark.parametrize("name", SITES)
def test_coastal_plain_zsed_and_region(name):
    (lon, lat), zsed, in_region = SITES[name]
    assert nshm_site_data.coastal_plain_zsed(lon, lat) == zsed
    assert nshm_site_data.in_coastal_plain_region(lon, lat) is in_region


def test_snapping_to_the_grid():
    # -89.6 and 36.6 are grid points; a site within half a grid spacing snaps to them
    assert nshm_site_data.coastal_plain_zsed(-89.6 + 0.024, 36.6 - 0.024) == 0.903


def test_coastal_plain_tree():
    tree = nshm_site_data.coastal_plain_stable_crust_tree()
    assert tree == {
        "NGA_EAST_2026": 0.1667,
        "NGA_EAST_2026_CPA": 0.1667,
        "NGA_EAST_2026_ADJUSTED": 0.3333,
        "NGA_EAST_SEEDS_2026": 0.0833,
        "NGA_EAST_SEEDS_2026_CPA": 0.0833,
        "NGA_EAST_SEEDS_2026_ADJUSTED": 0.1667,
    }


def test_nshm_site_parameters():
    site = {"latitude": 36.6, "longitude": -89.6}
    zsed, trees, notes = plha.nshm_site_parameters(site, {"stable_crust"})
    assert zsed == 0.903 and set(trees) == {"stable_crust"} and notes
    # explicit values override the NSHM value
    assert plha.nshm_site_parameters({**site, "zsed": 2.0}, {"stable_crust"})[0] == 2.0
    assert plha.nshm_site_parameters({**site, "zsed": None}, {"stable_crust"})[0] is None
    # no site data
    assert plha.nshm_site_parameters({**site, "nshm_site_data": False}, {"stable_crust"}) == (
        None, {}, [],
    )
    # the region tree is only used for the stable crust; no notes without stable crust models
    zsed, trees, notes = plha.nshm_site_parameters(site, {"active_crust"})
    assert zsed == 0.903 and trees == {} and notes == []

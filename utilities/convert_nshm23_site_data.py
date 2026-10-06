'''
Convert the site data of the USGS NSHM conterminous U.S. model (nshm-conus release 6.2.0) that
the stable crust ground motion models use into ucla_plha/site_data/coastal_plain.npz.

nshmp-lib (commit 44728a7d, gov.usgs.earthquake.nshmp.model.SiteData and calc.Sites) gives every
site of a hazard calculation, including the sites of the USGS hazard web service, which take only
a longitude, latitude, and Vs30, two location-dependent values from the model's site-data/
directory:

    site-data/margin/coastal-plain.{geojson,csv}   Gulf and Atlantic coastal plain sediment
        thickness zSed (Boyd, 2023, CPABasementDepth_v220517, https://doi.org/10.5066/P9EBOWU8),
        in m on a 0.05 degree grid inside a polygon. nshmp-lib snaps the site to the grid
        (Math.round(x / 0.05) * 0.05, rounded HALF_UP to 2 decimals), uses the value if the
        snapped location is inside the polygon and has data, and converts it to km rounded to
        3 decimals (HALF_UP). Sites without a value have zSed = NaN (off the coastal plain).
    site-data/gmm-region/coastal-plain.geojson   the "Coastal Plain CPA region": a polygon and
        the stable crust ground motion model logic tree used for all stable crust sources (fault,
        grid, and zone, including the system grid with its own tree) at sites inside it. The
        tree adds the Chapman and Guo (2021) coastal plain amplification models NGA_EAST_2026_CPA
        and NGA_EAST_SEEDS_2026_CPA.

The output file has:
    zsed_m: int16 array (n_lat, n_lon) of zSed in m rounded as nshmp-lib rounds the km value
        (-1 where the csv has no value)
    lon0, lat0, spacing: longitude and latitude of zsed_m[0, 0] and the grid spacing (degrees)
    margin_polygon, region_polygon: (n, 2) lon, lat vertices
    region_tree: JSON text of the gmm-region trees ({"stable-crust": {"all": [{"id", "weight"}]}})
    info: JSON text with the names, references, and sources

Usage (from the repository root): python utilities/convert_nshm23_site_data.py
'''
import json
import os
from decimal import ROUND_HALF_UP, Decimal

import numpy as np

from nshm23_ceus_common import NSHM_REF, NSHMP_LIB_COMMIT, fetch

OUTPUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'ucla_plha', 'site_data',
    'coastal_plain.npz')


def main():
    margin = json.load(open(fetch('site-data/margin/coastal-plain.geojson')))
    region = json.load(open(fetch('site-data/gmm-region/coastal-plain.geojson')))
    csv = fetch('site-data/margin/coastal-plain.csv')
    spacing = float(margin['properties']['spacing'])
    with open(csv) as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    assert lines[0].split(',') == ['lon', 'lat', 'zsed']
    data = np.array([line.split(',') for line in lines[1:]])
    lon = np.array([float(v) for v in data[:, 0]])
    lat = np.array([float(v) for v in data[:, 1]])
    # nshmp-lib: Maths.round(zsed / 1000.0, 3) (HALF_UP) km, stored here in m
    zsed_m = np.array([
        int(Decimal(repr(float(v) / 1000.0)).quantize(Decimal('0.001'), ROUND_HALF_UP) * 1000)
        for v in data[:, 2]
    ])
    lon0, lat0 = lon.min(), lat.min()
    i = np.rint((lat - lat0) / spacing).astype(int)
    j = np.rint((lon - lon0) / spacing).astype(int)
    assert np.allclose(lat0 + i * spacing, lat) and np.allclose(lon0 + j * spacing, lon)
    grid = np.full((i.max() + 1, j.max() + 1), -1, dtype=np.int16)
    assert zsed_m.max() < np.iinfo(np.int16).max
    grid[i, j] = zsed_m
    info = {
        'source': 'nshm-conus %s site-data/margin/coastal-plain.{geojson,csv} and '
                  'site-data/gmm-region/coastal-plain.geojson, as read by nshmp-lib %s '
                  '(model.SiteData, calc.Sites)' % (NSHM_REF, NSHMP_LIB_COMMIT),
        'margin': margin['properties'],
        'region': {'name': region['properties']['name']},
    }
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    np.savez_compressed(
        OUTPUT,
        zsed_m=grid,
        lon0=lon0,
        lat0=lat0,
        spacing=spacing,
        margin_polygon=np.array(margin['geometry']['coordinates'][0], dtype=float),
        region_polygon=np.array(region['geometry']['coordinates'][0], dtype=float),
        region_tree=json.dumps(region['properties']['gmm-trees']),
        info=json.dumps(info),
    )
    print('%s: %d zSed values, grid %s' % (OUTPUT, len(zsed_m), grid.shape))


if __name__ == '__main__':
    main()

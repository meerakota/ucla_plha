"""
Convert the USGS nshmp-lib point source rJB correction tables into the Numpy file used by
ucla_plha (src/ucla_plha/geometry/nshmp_rjb_tables.npz).

nshmp-lib (https://code.usgs.gov/ghsc/nshmp/nshmp-lib, commit 44728a7d) represents a point
source of unknown strike by the average Joyner-Boore distance of finite ruptures centered on
the point (RuptureScaling.pointSourceDistance). The averages are tabulated for magnitudes
6.05 to 8.55 (step 0.1) and epicentral distances 0 to 1000 km (step 1 km) in
src/main/resources/fault/surface/rjb_{wc94length,geomatrix,somerville}.dat, for the
NSHM_POINT_WC94_LENGTH, NSHM_SUB_GEOMAT_LENGTH, and NSHM_SOMERVILLE rupture scaling
relations. Each table is saved as a 26 x 1001 array (rows: magnitude, columns: distance).

Usage: python convert_nshmp_rjb_tables.py  (downloads the tables with curl)
"""

import os
import subprocess

import numpy as np

COMMIT = "44728a7d"
URL = (
    "https://code.usgs.gov/api/v4/projects/1356/repository/files/"
    "src%2Fmain%2Fresources%2Ffault%2Fsurface%2F{name}.dat/raw?ref=" + COMMIT
)
TABLES = {
    "nshm_point_wc94_length": "rjb_wc94length",
    "nshm_sub_geomat_length": "rjb_geomatrix",
    "nshm_somerville": "rjb_somerville",
}
OUTPUT = os.path.join(
    os.path.dirname(__file__), "..", "src", "ucla_plha", "geometry", "nshmp_rjb_tables.npz"
)


def read_table(text):
    """Parse an nshmp-lib rjb_*.dat file (blocks of '#Mag m' followed by 'r rjb' lines)."""
    table = np.zeros((26, 1001))
    i_mag = -1
    i_r = 0
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#Mag"):
            i_mag += 1
            i_r = 0
            continue
        if line.startswith("#"):
            continue
        r, rjb = line.split()
        assert float(r) == i_r
        table[i_mag, i_r] = float(rjb)
        i_r += 1
    assert i_mag == 25 and i_r == 1001
    return table


def main():
    tables = {}
    for key, name in TABLES.items():
        text = subprocess.run(
            ["curl", "-s", URL.format(name=name)], capture_output=True, text=True, check=True
        ).stdout
        tables[key] = read_table(text)
    tables["magnitudes"] = 6.05 + 0.1 * np.arange(26)
    tables["distances"] = np.arange(1001, dtype=float)
    np.savez_compressed(OUTPUT, **tables)


if __name__ == "__main__":
    main()

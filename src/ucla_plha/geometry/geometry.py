import numpy as np
import pandas as pd

def point_to_xyz(p):
    '''
    Convert point latitude, longitude, elevation values to Cartesian coordinates.
    p = N x 3 Numpy array, where N is the number of points. 
    Each row of p has [latitude, longitude, elevation], where latitude and longitude are in degrees and elevation is in km.
    r = Earth's radius in km.
    '''
    p_xyz = np.empty(3)
    rad = np.pi/180.0
    lat = p[0]
    lon = p[1]
    elevation = p[2]
    a = 6378.1370 # Earth's equatorial radius in km
    b = 6356.7523 # Earth's polar radius in km
    r = np.sqrt(((a**2 * np.cos(lat * rad))**2 + (b**2 * np.sin(lat * rad))**2) / ((a * np.cos(lat * rad))**2 + (b * np.sin(lat * rad))**2))
    p_xyz[0] = (
        (r - elevation)
        * np.cos(lat * rad)
        * np.cos(lon * rad)
    )
    p_xyz[1] = (
        (r - elevation)
        * np.cos(lat * rad)
        * np.sin(lon * rad)
    )
    p_xyz[2] = (r - elevation) * np.sin(lat * rad)
    return np.asarray(p_xyz)

def point_triangle_distance(tri_xyz, p_xyz, fault_id):
    '''
    Compute distance between all of the triangles for a source model (tri_xyz) and the point p_xyz.
    Return a Numpy array of distances and associated fault_id's.
    tri_xyz = N x 3 x 3 Numpy array of points defining triangle in Cartesian coordinates
    p_xyz = Numpy array defining point in Cartesian coordinates
    fault_id = N x 1 Numpy array of integers defining the fault id for each triangle
    '''
    # Perform checks that all three points are different. Add 1m to elevation for points that have the same coordinates.
    for i in range(len(tri_xyz)):
        if((tri_xyz[i,0] == tri_xyz[i,1]).all() or (tri_xyz[i,0] == tri_xyz[i,2]).all()):
            tri_xyz[i,0,2] = tri_xyz[i,0,2] + 0.001
        if(tri_xyz[i,1] == tri_xyz[i,2]).all():
            tri_xyz[i,1,2] = tri_xyz[i,1,2] + 0.001

    e0 = tri_xyz[:, 1] - tri_xyz[:, 0]
    e1 = tri_xyz[:, 2] - tri_xyz[:, 0]
    a = np.sum(np.multiply(e0, e0), axis=1)
    b = np.sum(np.multiply(e0, e1), axis=1)
    c = np.sum(np.multiply(e1, e1), axis=1)

    det = a * c - b * b
    det[det==0] = 1.0e-8

    d0 = tri_xyz[:, 0] - p_xyz

    d = np.sum(np.multiply(e0, d0), axis=1)
    e = np.sum(np.multiply(e1, d0), axis=1)
    f = np.sum(np.multiply(d0, d0), axis=1)

    s = b * e - c * d
    t = b * d - a * e

    sqrdistance = np.empty(len(tri_xyz), dtype=float)

    # Region 4
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d < 0.0) & (-d >= a)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d < 0.0) & (-d < a)
    sqrdistance[cond] = (-d * d / a + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d >= 0.0) & (e >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d >= 0.0) & (e < 0.0) & (-e >= c)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t < 0.0) & (d >= 0.0) & (e < 0.0) & (-e < c)
    sqrdistance[cond] = (-e * e / c + f)[cond]

    # Region 3
    cond = (s + t <= det) & (s < 0.0) & (t >= 0.0) & (e >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t >= 0.0) & (e < 0.0) & (-e >= c)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t <= det) & (s < 0.0) & (t >= 0.0) & (e < 0.0) & (-e < c)
    sqrdistance[cond] = (-e * e / c + f)[cond]

    # Region 5
    cond = (s + t <= det) & (s >= 0.0) & (t < 0.0) & (d >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t <= det) & (s >= 0.0) & (t < 0.0) & (d < 0.0) & (-d >= a)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t <= det) & (s >= 0.0) & (t < 0.0) & (d < 0.0) & (-d < a)
    sqrdistance[cond] = (-d * d / a + f)[cond]

    # Region 0
    invDet = 1.0 / det
    stemp = s * invDet
    ttemp = t * invDet
    cond = (s + t <= det) & (s >= 0) & (t >= 0)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]

    # Region 2
    tmp0 = b + d
    tmp1 = c + e
    numer = tmp1 - tmp0
    denom = a - 2.0 * b + c
    denom[denom==0] = 1.0e-8
    stemp = numer / denom
    ttemp = 1.0 - stemp
    cond = (s + t > det) & (s < 0.0) & (tmp1 > tmp0) & (numer >= denom)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 > tmp0) & (numer < denom)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 <= tmp0) & (tmp1 <= 0.0)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 <= tmp0) & (tmp1 > 0.0) & (e >= 0.0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t > det) & (s < 0.0) & (tmp1 <= tmp0) & (tmp1 > 0.0) & (e < 0.0)
    sqrdistance[cond] = (-e * e / c + f)[cond]

    # Region 6
    tmp0 = b + e
    tmp1 = a + d
    numer = tmp1 - tmp0
    denom = a - 2.0 * b + c
    denom[denom==0] = 1.0e-8
    ttemp = numer / denom
    stemp = 1.0 - ttemp
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 > tmp0) & (numer >= denom)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 > tmp0) & (numer < denom)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 <= tmp0) & (tmp1 <= 0)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 <= tmp0) & (tmp1 > 0) & (d >= 0)
    sqrdistance[cond] = (f)[cond]
    cond = (s + t > det) & (s >= 0) & (t < 0) & (tmp1 <= tmp0) & (tmp1 > 0) & (d < 0)
    sqrdistance[cond] = (-d * d / a + f)[cond]

    # Region 1
    numer = c + e - b - d
    denom = a - 2.0 * b + c
    denom[denom==0] = 1.0e-8
    stemp = numer / denom
    ttemp = 1.0 - stemp
    cond = (s + t > det) & (s >= 0) & (t >= 0) & (numer <= 0)
    sqrdistance[cond] = (c + 2.0 * e + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t >= 0) & (numer > 0) & (numer >= denom)
    sqrdistance[cond] = (a + 2.0 * d + f)[cond]
    cond = (s + t > det) & (s >= 0) & (t >= 0) & (numer > 0) & (numer < denom)
    sqrdistance[cond] = (
        stemp * (a * stemp + b * ttemp + 2.0 * d)
        + ttemp * (b * stemp + c * ttemp + 2.0 * e)
        + f
    )[cond]

    # account for numerical round-off error
    sqrdistance[sqrdistance <= 0] = 0

    # use Pandas groupby function to find shortest distance for a particular fault rupture
    df = pd.DataFrame()
    df['fault_id'] = fault_id
    df['sqrdistance'] = sqrdistance
    sqrdistance_out = df.groupby('fault_id')['sqrdistance'].min().values

    # return distance array
    return np.sqrt(sqrdistance_out).T

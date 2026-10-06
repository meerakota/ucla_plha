"""
Convert the UCERF3 branch-averaged fault system solutions (FM3_1_branch_averaged and
FM3_2_branch_averaged, from OpenSHA) to the ucla_plha fault source model files in
src/ucla_plha/source_models/fault_source_models/ucerf3_fm31 and ucerf3_fm32.

Aseismicity: UCERF3 reduces the seismogenic area of each fault subsection by its aseismic slip
factor (Field et al. 2014, OpenSHA FaultSection.getReducedAveUpperDepth and
getReducedDownDipWidth): the upper depth moves down by
aseismic_slip_factor * (lower_depth - upper_depth), the lower depth does not change, and the
coupling coefficient reduces the slip rate. The rupture areas in ruptures/properties.csv, and
therefore the rupture magnitudes, are the sums of the reduced subsection areas, and the OpenSHA
UCERF3 ERF builds the rupture surfaces with the reduced subsections (FaultSystemRupSet
.getSurfaceForRupture with aseisReducesArea = true, the default of BaseFaultSystemSolutionERF):
StirlingGriddedSurface moves the upper edge of the subsection from its trace down-dip to the
reduced upper depth, i.e. it removes the top of the dipping fault plane, and keeps the lower
edge. The USGS NSHM (nshmp-lib, DefaultGriddedSurface, and the rupture depths and widths of the
nshm-conus UCERF3 ruptures.csv) does the same. The fault_sections.geojson files give the
original (unreduced) upper depths ("UpDepth", the depth of the trace) and the aseismic slip
factors ("AseismicSlipFactor").

With APPLY_ASEISMICITY = True (the default), the conversion applies the reduction in the same
way: the top corners of each fault segment are moved from the trace down-dip by
aseismic_slip_factor * (lower_depth - upper_depth) / tan(dip) in the dip direction, at the
reduced upper depth, the bottom corners do not move, the rupture ztor is the shallowest reduced
upper depth of the rupture's subsections, and the rupture dip is weighted by the reduced
subsection areas. With APPLY_ASEISMICITY = False (environment variable
UCERF3_APPLY_ASEISMICITY=0), the conversion uses the full (unreduced) subsection planes, as
ucla_plha versions up to 2.0.0 did.

Run from the utilities directory. Environment variables:
    UCERF3_APPLY_ASEISMICITY: 1 (default) or 0
    UCERF3_RUPTURE_DEPTHS: shallowest (default) or area_weighted (see RUPTURE_DEPTHS)
    UCERF3_OUTPUT_DIR: directory for the ucerf3_fm31 and ucerf3_fm32 directories
        (default ../src/ucla_plha/source_models/fault_source_models)
ruptures_segments.npz is made from ruptures/indices.csv (not in the repository because of its
size; it is in the OpenSHA solution zip files). If indices.csv is not present, the existing
ruptures_segments.npz in the output directory is used (it does not depend on aseismicity).
"""
import json
import os
import shutil

import numpy as np
import pandas as pd

APPLY_ASEISMICITY = os.environ.get("UCERF3_APPLY_ASEISMICITY", "1") not in ("0", "false", "False")
# Rupture depths: "shallowest" (default): ztor is the shallowest upper depth and zbor the deepest
# lower depth of the rupture's subsections. "area_weighted": the area-weighted averages, as OpenSHA
# (CompoundSurface.getAveRupTopDepth, getAveRupBottomDepth) and the nshm-conus UCERF3 ruptures.csv.
RUPTURE_DEPTHS = os.environ.get("UCERF3_RUPTURE_DEPTHS", "shallowest")
OUTPUT_DIR = os.environ.get(
    "UCERF3_OUTPUT_DIR", "../src/ucla_plha/source_models/fault_source_models"
)

rad = np.pi / 180
a_earth = 6378.1370  # Earth's equatorial radius in km
b_earth = 6356.7523  # Earth's polar radius in km


def _radius(lat):
    """Radius of the oblate spheroid (km) at latitude lat (degrees)."""
    return np.sqrt(
        ((a_earth**2 * np.cos(lat * rad)) ** 2 + (b_earth**2 * np.sin(lat * rad)) ** 2)
        / ((a_earth * np.cos(lat * rad)) ** 2 + (b_earth * np.sin(lat * rad)) ** 2)
    )


def offset_points(lat, lon, w, dip_dir):
    """
    Move points (lat, lon, degrees) a horizontal distance w (km) in the direction dip_dir
    (degrees clockwise from north). The north component is w * cos(dip_dir), and the longitude is
    then chosen so that the great circle distance from the original point is w.
    """
    wy = w * np.cos(dip_dir * rad)
    dlat = wy / _radius(lat) / rad
    lat_new = lat + dlat
    r_new = _radius(lat_new)
    invangle = 1.0 - (2.0 * np.sin(w / (2.0 * r_new)) ** 2.0 + np.cos(dlat * rad) - 1.0) / (
        np.cos(lat * rad) * np.cos(lat_new * rad)
    )
    # floating point precision may render inverse angles higher than 1.0 (or less than -1.0, which
    # doesn't happen here but is mathematically possible). So impose ranges on the inverse angles
    invangle = np.clip(invangle, -1.0, 1.0)
    # arccos is always positive, so use the sign of the east component of the dip direction to move
    # the points west for faults that dip toward the west (dip_dir between 180 and 360 degrees)
    lon_sign = np.sign(np.sin(dip_dir * rad))
    lon_new = lon + lon_sign * np.arccos(invangle) / rad
    return lat_new, lon_new


def aseismic_reduction(upper_depth, lower_depth, aseismicity, dip):
    """
    UCERF3 (OpenSHA FaultSection.getReducedAveUpperDepth) aseismic reduction of a subsection.
    Returns the reduced upper depth (km) and the horizontal distance (km) from the trace to the
    reduced upper edge in the dip direction. The distance is zero for vertical subsections.
    """
    reduced_upper = upper_depth + aseismicity * (lower_depth - upper_depth)
    shift = np.where(
        dip >= 90.0, 0.0, (reduced_upper - upper_depth) / np.tan(np.minimum(dip, 89.999) * rad)
    )
    return reduced_upper, shift


def get_lat_lon(filename, apply_aseismicity=True):
    """
    Read UCERF3 source geojson file and output a Python list containing Numpy arrays for segment_id,
    and the latitudes, longitudes, and depths of the four corners of each fault segment (the
    segments are the pieces of the subsection traces between consecutive trace points): corners 1
    and 2 are the upper edge and corners 3 and 4 the lower edge.

    Notation: UCERF3 uses the variable name "FaultID" to refer to the subsections of the fault.
    We reserve the word "fault" for a particular named fault (e.g., Airport Lake), and "section" for
    a section of the fault. We therefore have renamed "FaultID" to "segment_id" in our dataframe.

    The trace is at the (unreduced) upper depth, and the lower edge is down-dip of the trace at the
    lower depth. With apply_aseismicity, the upper edge is the trace moved down-dip to the reduced
    upper depth (see aseismic_reduction); otherwise it is the trace.
    """
    features = json.load(open(filename))["features"]

    lat1 = []
    lon1 = []
    lat2 = []
    lon2 = []
    segment_id = []
    dip = []
    dip_dir = []
    lower_depth = []
    upper_depth = []
    aseismicity = []

    for f in features:
        p = f["properties"]
        gm_array = np.asarray(f["geometry"]["coordinates"])
        for i in range(gm_array.shape[0] - 1):
            segment_id.append(p["FaultID"])
            lon1.append(gm_array[i, 0])
            lat1.append(gm_array[i, 1])
            lon2.append(gm_array[i + 1, 0])
            lat2.append(gm_array[i + 1, 1])
            dip.append(p["DipDeg"])
            dip_dir.append(p["DipDir"])
            lower_depth.append(p["LowDepth"])
            upper_depth.append(p["UpDepth"])
            aseismicity.append(p.get("AseismicSlipFactor", 0.0))

    segment_id = np.asarray(segment_id, dtype="intc")
    lat1 = np.asarray(lat1)
    lon1 = np.asarray(lon1)
    lat2 = np.asarray(lat2)
    lon2 = np.asarray(lon2)
    dip = np.asarray(dip)
    dip_dir = np.asarray(dip_dir)
    lower_depth = np.asarray(lower_depth)
    upper_depth = np.asarray(upper_depth)
    aseismicity = np.asarray(aseismicity)

    # lower edge: down-dip of the trace (which is at the original upper depth)
    w = (lower_depth - upper_depth) / np.tan(dip * rad)
    lat3, lon3 = offset_points(lat1, lon1, w, dip_dir)
    lat4, lon4 = offset_points(lat2, lon2, w, dip_dir)

    # upper edge
    if apply_aseismicity:
        top_depth, shift = aseismic_reduction(upper_depth, lower_depth, aseismicity, dip)
        moved = shift > 0.0
        lat1, lon1 = (
            np.where(moved, v, v0) for v, v0 in zip(offset_points(lat1, lon1, shift, dip_dir), (lat1, lon1))
        )
        lat2, lon2 = (
            np.where(moved, v, v0) for v, v0 in zip(offset_points(lat2, lon2, shift, dip_dir), (lat2, lon2))
        )
    else:
        top_depth = upper_depth

    return (
        segment_id,
        lat1,
        lon1,
        top_depth,
        lat2,
        lon2,
        top_depth,
        lat3,
        lon3,
        lower_depth,
        lat4,
        lon4,
        lower_depth,
        dip,
    )


def latlonel_to_xyz(geom):
    """
    Accept N x M x 3 Numpy array of lat, lon, depth data, where N is the number of geometric objects,
    M is the number of points per geometric object, and 3 are the lat, lon, depth for the point.
    Convert to Cartesian coordinates assuming the earth is an oblate spheroid.
    Return N x M x 3 Numpy array of points in Cartesian coordinates
    """
    xyz = np.empty(geom.shape)
    d = geom[:, :, 2]
    lat = geom[:, :, 0]
    lon = geom[:, :, 1]
    r = _radius(lat)
    xyz[:, :, 0] = (r - d) * np.cos(lat * rad) * np.cos(lon * rad)
    xyz[:, :, 1] = (r - d) * np.cos(lat * rad) * np.sin(lon * rad)
    xyz[:, :, 2] = (r - d) * np.sin(lat * rad)
    return xyz


def get_triangles(lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4):
    """
    Accept Numpy arrays of latitude, longitude, and depth for four points defining the segment corners
    Return a float array of triangles for computing Rrup, and a float array of triangles for computing Rjb
    """
    tri1_rrup = np.asarray([[lat1, lat2, lat4], [lon1, lon2, lon4], [d1, d2, d4]]).T
    tri2_rrup = np.asarray([[lat1, lat3, lat4], [lon1, lon3, lon4], [d1, d3, d4]]).T
    tri_rrup = np.concatenate((tri1_rrup, tri2_rrup))
    tri_rrup_xyz = latlonel_to_xyz(tri_rrup)
    zero_depth = np.zeros(len(lat1))
    tri1_rjb = np.asarray(
        [[lat1, lat2, lat4], [lon1, lon2, lon4], [zero_depth, zero_depth, zero_depth]]
    ).T
    tri2_rjb = np.asarray(
        [[lat1, lat3, lat4], [lon1, lon3, lon4], [zero_depth, zero_depth, zero_depth]]
    ).T
    tri_rjb = np.concatenate((tri1_rjb, tri2_rjb))
    tri_rjb_xyz = latlonel_to_xyz(tri_rjb)
    return (tri_rrup_xyz, tri_rjb_xyz)


def get_rectangles(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4):
    """
    Accept lat-lon array from the get_lat_lon() function, and return a float array of rectangles
    for the surface projection of the fault segments for computing Rx, Rx1, and Ry0
    """
    zero_depth = np.zeros(len(lat1))
    rect_rjb = np.asarray(
        [
            [lat1, lat2, lat3, lat4],
            [lon1, lon2, lon3, lon4],
            [zero_depth, zero_depth, zero_depth, zero_depth],
        ]
    ).T
    rect_rjb_xyz = latlonel_to_xyz(rect_rjb)
    return rect_rjb_xyz


def get_section_properties(filename, apply_aseismicity=True):
    """
    Read UCERF3 source geojson file and return Numpy arrays of upper depth, lower depth, dip, and
    down-dip area for each section, indexed by segment_id (FaultID). With apply_aseismicity, the
    upper depth and area are reduced by the aseismic slip factor (the area is then the area
    UCERF3 uses for the rupture magnitudes, "Area (m^2)" in ruptures/properties.csv).
    """
    features = json.load(open(filename))["features"]
    segment_id = np.asarray([f["properties"]["FaultID"] for f in features])
    # get_rupture_data() indexes these arrays by segment_id, so segment_id must be 0, 1, ..., N-1
    assert np.array_equal(segment_id, np.arange(len(segment_id)))
    upper_depth = np.asarray([f["properties"]["UpDepth"] for f in features])
    lower_depth = np.asarray([f["properties"]["LowDepth"] for f in features])
    dip = np.asarray([f["properties"]["DipDeg"] for f in features])
    aseismicity = np.asarray([f["properties"].get("AseismicSlipFactor", 0.0) for f in features])
    # trace length of each section from the haversine distance between consecutive points
    length = np.zeros(len(features))
    for i, f in enumerate(features):
        coords = np.asarray(f["geometry"]["coordinates"])
        lon = coords[:, 0] * rad
        lat = coords[:, 1] * rad
        r = _radius(coords[:, 1])
        hav = (
            np.sin(np.diff(lat) / 2.0) ** 2
            + np.cos(lat[:-1]) * np.cos(lat[1:]) * np.sin(np.diff(lon) / 2.0) ** 2
        )
        length[i] = np.sum(0.5 * (r[:-1] + r[1:]) * 2.0 * np.arcsin(np.sqrt(hav)))
    if apply_aseismicity:
        upper_depth, _ = aseismic_reduction(upper_depth, lower_depth, aseismicity, dip)
    area = length * (lower_depth - upper_depth) / np.sin(dip * rad)
    return (upper_depth, lower_depth, dip, area)


def get_rupture_data(
    rupture_file,
    rate_file,
    ruptures_segments_file,
    section_file,
    output_file,
    apply_aseismicity=True,
    rupture_depths="shallowest",
):
    """
    Read UCERF3 rupture data file, rate data file, and section properties. Save magnitude, rate, style of
    faulting, dip, ztor, and zbor for each rupture in compressed npz format. ztor is the shallowest upper
    depth and zbor is the deepest lower depth of the sections in the rupture, and dip is the area-weighted
    average dip of the sections in the rupture (upper depths and areas reduced by aseismicity with
    apply_aseismicity). With rupture_depths="area_weighted", ztor and zbor are the area-weighted
    averages of the section upper and lower depths instead.
    """
    rupture_df = pd.read_csv(rupture_file)
    rate_df = pd.read_csv(rate_file)
    ruptures_segments = np.load(ruptures_segments_file)
    segment_index = ruptures_segments["segment_index"]
    ruptures_index = ruptures_segments["rupture_index"]
    upper_depth, lower_depth, dip_section, area_section = get_section_properties(
        section_file, apply_aseismicity
    )
    split_indices = np.where(np.diff(ruptures_index) != 0)[0] + 1
    boundaries = np.r_[0, split_indices]
    area_all = area_section[segment_index]
    if rupture_depths == "shallowest":
        ztor = np.minimum.reduceat(upper_depth[segment_index], boundaries)
        zbor = np.maximum.reduceat(lower_depth[segment_index], boundaries)
    elif rupture_depths == "area_weighted":
        area_sum = np.add.reduceat(area_all, boundaries)
        ztor = np.add.reduceat(upper_depth[segment_index] * area_all, boundaries) / area_sum
        zbor = np.add.reduceat(lower_depth[segment_index] * area_all, boundaries) / area_sum
    else:
        raise ValueError(f'rupture_depths must be "shallowest" or "area_weighted", not "{rupture_depths}"')
    dip = np.add.reduceat(
        dip_section[segment_index] * area_all, boundaries
    ) / np.add.reduceat(area_all, boundaries)
    # averaging vertical sections can give dips slightly above 90 from floating point round-off
    dip = np.minimum(dip, 90.0)
    rate = rate_df["Annual Rate"].values
    fault_type = np.full(len(rupture_df), 1)
    rake = rupture_df["Average Rake (degrees)"].values
    m = rupture_df["Magnitude"].values
    fault_type[(rake > -150) & (rake < -30)] = 2
    fault_type[(rake >= -180) & (rake <= -150)] = 3
    fault_type[(rake >= -30) & (rake <= 30)] = 3
    fault_type[(rake >= 150) & (rake <= 180)] = 3
    np.savez_compressed(
        output_file,
        m=m,
        rate=rate,
        fault_type=fault_type,
        dip=dip,
        ztor=ztor,
        zbor=zbor,
    )

    return


def get_ruptures_segments(rupture_indices_file, output_file):
    """
    Read UCERF3 file containing the list of all of the segments associated with each rupture.
    Organize the data into a single Numpy array containing rupture_index and segment_index.
    The array is very large, so use the smallest possible integer container, and compress the file.
    """
    indices = pd.read_csv(rupture_indices_file, engine="python")
    rupture_index = indices["Rupture Index"].values
    segment_index = indices["Num Sections"].values
    indicesT = indices.drop(
        labels=["Rupture Index", "Num Sections"], axis=1
    ).transpose()
    segment_index_array = np.empty(len(rupture_index), dtype="object")
    rupture_index_array = np.empty(len(rupture_index), dtype="object")
    for i, it in indicesT.items():
        segment_index_array[i] = np.asarray(
            it.values[0 : segment_index[i]], dtype="short"
        )
        rupture_index_array[i] = np.full(segment_index[i], rupture_index[i])
    rupture_index_all = np.array(np.hstack(rupture_index_array), dtype=np.int32)
    segment_index_all = np.array(np.hstack(segment_index_array), dtype=np.int32)
    np.savez_compressed(
        output_file, rupture_index=rupture_index_all, segment_index=segment_index_all
    )


def convert(
    input_dir,
    model,
    output_dir=OUTPUT_DIR,
    apply_aseismicity=APPLY_ASEISMICITY,
    rupture_depths=RUPTURE_DEPTHS,
):
    """Convert one UCERF3 branch-averaged solution (input_dir) to the ucla_plha model files."""
    out = os.path.join(output_dir, model)
    os.makedirs(out, exist_ok=True)
    section_file = os.path.join(input_dir, "ruptures", "fault_sections.geojson")

    ### Compute triangles representing fault segments, and array of segment_id values
    segment_id, lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4, dip = (
        get_lat_lon(section_file, apply_aseismicity)
    )
    tri = get_triangles(lat1, lon1, d1, lat2, lon2, d2, lat3, lon3, d3, lat4, lon4, d4)
    np.save(os.path.join(out, "tri_segment_id.npy"), np.concatenate((segment_id, segment_id)))
    np.save(os.path.join(out, "tri_rrup.npy"), tri[0])
    np.save(os.path.join(out, "tri_rjb.npy"), tri[1])
    rect = get_rectangles(lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4)
    np.save(os.path.join(out, "rect_segment_id.npy"), segment_id)
    np.save(os.path.join(out, "rect_rjb.npy"), rect)

    ### Compute array mapping rupture and segment indices
    ruptures_segments_file = os.path.join(out, "ruptures_segments.npz")
    indices_file = os.path.join(input_dir, "ruptures", "indices.csv")
    if os.path.exists(indices_file):
        get_ruptures_segments(indices_file, ruptures_segments_file)
    else:
        default = os.path.join(
            "../src/ucla_plha/source_models/fault_source_models", model, "ruptures_segments.npz"
        )
        if not os.path.exists(ruptures_segments_file):
            shutil.copy(default, ruptures_segments_file)
        print(f"{indices_file} not found; using {ruptures_segments_file}")

    ### Compute ruptures.npz file that contains magnitude, rate, style of faulting, dip, ztor, zbor
    get_rupture_data(
        os.path.join(input_dir, "ruptures", "properties.csv"),
        os.path.join(input_dir, "solution", "rates.csv"),
        ruptures_segments_file,
        section_file,
        os.path.join(out, "ruptures.npz"),
        apply_aseismicity,
        rupture_depths,
    )


if __name__ == "__main__":
    print(
        f"apply_aseismicity = {APPLY_ASEISMICITY}, rupture_depths = {RUPTURE_DEPTHS}, "
        f"output directory {OUTPUT_DIR}"
    )
    convert("FM3_1_branch_averaged", "ucerf3_fm31")
    convert("FM3_2_branch_averaged", "ucerf3_fm32")

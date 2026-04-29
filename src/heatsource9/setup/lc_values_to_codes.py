from pathlib import Path

from heatsource9.io.input_files import read_to_list, write_input
from heatsource9.setup.constants import sheetnames
from heatsource9.setup.headers import headers_lcdata


def lc_values_to_codes(
    lcdatafile,
    output_lccodefile = None,
    output_lcdatafile = None,
    canopy_data = "CanopyCover",
    trans_count = None,
    transsample_count = None,
    heatsource8 = False,
    overwrite = False,
):
    """
    Convert an old values style land cover data file to land cover codes.

    Each land cover sample is assigned a unique code built from the NODE_ID and
    land cover header. For example, Node 200 at transect 1, sample 3 is assigned
    the code `200_T1_S3`. The output land cover codes file stores the pivoted
    vegetation attributes for those codes in the appropriate columns. The output
    land cover data file stores the codes and sample elevations.

    Land cover value columns are read positionally after the first eight land
    cover data columns, using the same transect and zone structure as the old
    land cover values headers. If canopy depth columns are missing, vegetation
    height is used as canopy depth.
    """
    if trans_count is None or transsample_count is None:
        msg = "trans_count and transsample_count are required for land cover values conversion."
        raise ValueError(msg)

    lcdata_path = Path(lcdatafile)
    output_lccode_path = (
        Path(output_lccodefile)
        if output_lccodefile is not None
        else lcdata_path.with_name(lcdata_path.stem + "_lccodes" + lcdata_path.suffix)
    )
    output_lcdata_path = (
        Path(output_lcdatafile)
        if output_lcdatafile is not None
        else lcdata_path.with_name(lcdata_path.stem + "_codes" + lcdata_path.suffix)
    )

    for path in (output_lccode_path, output_lcdata_path):
        if path.exists() and not overwrite:
            msg = "Output file already exists: {0}".format(path)
            raise FileExistsError(msg)

    rows = read_to_list(lcdata_path, skiprows=0, skipcols=0, sheetname=sheetnames["lcdatafile"])
    headers = list(rows[0])
    data_rows = rows[1:]
    base_headers = headers_lcdata({})[:8]

    trans_count = trans_count or 0
    transsample_count = transsample_count or 0
    heatsource8 = heatsource8 or False

    tran = ["NE", "E", "SE", "S", "SW", "W", "NW"] if heatsource8 else ["T{0}".format(x) for x in range(1, trans_count + 1)]
    zones = list(range(1, transsample_count + 1))

    sample_suffixes = []
    elevation_suffixes = []
    for ti, t in enumerate(tran):
        for z in zones:
            suffix = "{0}_S{1}".format(t, z)
            if ti == 0 and z == 1:
                sample_suffixes.append("T0_S0")
            sample_suffixes.append(suffix)
            elevation_suffixes.append(suffix)

    prefix = ["HT", "ELE", "LAI", "k", "OH"] if canopy_data == "LAI" else ["HT", "ELE", "CAN", "OH"]
    column_index = {}
    column = len(base_headers)
    for p in prefix:
        suffixes = elevation_suffixes if p == "ELE" else sample_suffixes
        for suffix in suffixes:
            column_index[p + "_" + suffix] = column
            column += 1

    has_canopy_depth = len(headers) >= column + len(sample_suffixes)
    if has_canopy_depth:
        for suffix in sample_suffixes:
            column_index["CD_" + suffix] = column
            column += 1

    lcdata_headers = list(base_headers)
    lcdata_headers.extend(["LC_" + suffix for suffix in sample_suffixes])
    lcdata_headers.extend(["ELE_" + suffix for suffix in elevation_suffixes])

    if canopy_data == "LAI":
        lccode_headers = ["NAME", "CODE", "HEIGHT", "LAI", "k", "OVERHANG", "CANOPY_DEPTH"]
    else:
        lccode_headers = ["NAME", "CODE", "HEIGHT", "CANOPY", "OVERHANG", "CANOPY_DEPTH"]

    lcdata_rows = []
    lccode_rows = []
    for row in data_rows:
        lcdata_row = list(row[:len(base_headers)])
        node_id = row[1]
        codes_by_suffix = {}
        for suffix in sample_suffixes:
            code = "{0}_{1}".format(node_id, suffix)
            codes_by_suffix[suffix] = code
            height = float(row[column_index["HT_" + suffix]])
            overhang = float(row[column_index["OH_" + suffix]])
            if "CD_" + suffix in column_index:
                canopy_depth = float(row[column_index["CD_" + suffix]])
            else:
                canopy_depth = height
            if canopy_data == "LAI":
                lccode_rows.append(
                    [
                        code,
                        code,
                        height,
                        float(row[column_index["LAI_" + suffix]]),
                        float(row[column_index["k_" + suffix]]),
                        overhang,
                        canopy_depth,
                    ]
                )
            else:
                lccode_rows.append(
                    [
                        code,
                        code,
                        height,
                        float(row[column_index["CAN_" + suffix]]),
                        overhang,
                        canopy_depth,
                    ]
                )
        lcdata_row.extend([codes_by_suffix[suffix] for suffix in sample_suffixes])
        lcdata_row.extend([float(row[column_index["ELE_" + suffix]]) for suffix in elevation_suffixes])
        lcdata_rows.append(lcdata_row)

    write_input(
        path=output_lcdata_path,
        headers=lcdata_headers,
        rows=lcdata_rows,
        sheetname=sheetnames["lcdatafile"],
        csv_mode=output_lcdata_path.suffix.lower() == ".csv",
    )
    write_input(
        path=output_lccode_path,
        headers=lccode_headers,
        rows=lccode_rows,
        sheetname=sheetnames["lccodefile"],
        csv_mode=output_lccode_path.suffix.lower() == ".csv",
    )

    return {
        "lcdatafile": output_lcdata_path,
        "lccodefile": output_lccode_path,
        "row_count": len(lcdata_rows),
        "code_count": len(lccode_rows),
    }

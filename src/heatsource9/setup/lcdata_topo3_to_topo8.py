from pathlib import Path

from heatsource9.io.input_files import read_to_list, write_input
from heatsource9.setup.constants import sheetnames


def lcdata_topo3_to_topo8( lcdatafile, output_lcdatafile = None, overwrite = False):
    """
    Convert a legacy lcdata file that has only three topographic shade angle directions
    columns to the current eight direction lcdata format.

    Old topo columns: TOPO_W, TOPO_S, TOPO_E
    New topo columns: TOPO_NE, TOPO_E, TOPO_SE, TOPO_S, TOPO_SW, TOPO_W, TOPO_NW, TOPO_N

    Input columns are read positionally. The first five base lcdata columns are
    preserved. Legacy TOPO_W, TOPO_S, and TOPO_E values are written to their
    matching columns in the current lcdata format, and the new topo direction
    columns are left blank. All other land cover and elevation columns after the
    legacy topo columns are preserved. If output_lcdata_path is not supplied,
    the existing lcdata file name will be used with "_topo8" added at the end.

    Note that the model will not run if the new topo direction columns are left blank.
    The model requires all eight topo columns to have data. If you don't have data, consider keeping
    the existing lcdata file with three topo columns. The model will still accept them.
    """

    lcdata_path = Path(lcdatafile)
    output_lcdata_path = (
        Path(output_lcdatafile)
        if output_lcdatafile is not None
        else lcdata_path.with_name(lcdata_path.stem + "_topo8" + lcdata_path.suffix)
    )

    if output_lcdata_path.exists() and not overwrite:
        msg = "Output file already exists: {0}".format(output_lcdata_path)
        raise FileExistsError(msg)

    rows = read_to_list(lcdata_path, skiprows=0, skipcols=0, sheetname=sheetnames["lcdatafile"])
    headers = list(rows[0])
    data_rows = rows[1:]

    output_headers = (
        list(headers[:5])
        + ["TOPO_NE", "TOPO_E", "TOPO_SE", "TOPO_S", "TOPO_SW", "TOPO_W", "TOPO_NW", "TOPO_N"]
        + list(headers[8:])
    )

    output_rows = []
    for row in data_rows:
        output_row = (
            list(row[:5])
            + [None, row[7], None, row[6], None, row[5], None, None]
            + list(row[8:])
        )
        output_rows.append(output_row)

    write_input(
        path=output_lcdata_path,
        headers=output_headers,
        rows=output_rows,
        sheetname=sheetnames["lcdatafile"],
        csv_mode=output_lcdata_path.suffix.lower() == ".csv",
    )

    return {
        "lcdatafile": output_lcdata_path,
        "row_count": len(output_rows),
    }

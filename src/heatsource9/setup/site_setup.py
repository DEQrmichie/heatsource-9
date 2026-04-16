from pathlib import Path

from heatsource9.io.input_files import read_to_dict
from heatsource9.setup.constants import sheetnames
from heatsource9.setup.headers import headers_met_sites, headers_trib_sites
from heatsource9.setup.input_setup import InputSetup


def _validate_site_file_names(site_file_name, file_names):
    """
    This checks that the file names are all unique or all the same.
    """
    unique_names = list(dict.fromkeys(file_names))
    if 1 < len(unique_names) < len(file_names):
        msg = (
            "There are some duplicate file names in {0}. Use either the same file name "
            "for all sites or a different file name for each site."
        ).format(site_file_name)
        raise ValueError(msg)


def _get_met_sites(model_path, control_params, control_path, run_type, ext):
    met_rows = []
    met_params = {}
    met_count = int(control_params.get("metsites") or 0)

    if met_count <= 0:
        result = (met_rows, met_params)
        return result

    met_path = model_path / control_params["metsitesfile"]
    if not met_path.exists():
        result = (met_rows, met_params)
        return result

    setup = InputSetup(control_path)
    headers = headers_met_sites()
    data = read_to_dict(
        path=met_path,
        colnames=headers,
        sheetname=sheetnames["metsitesfile"],
        value_check=setup._validate,
        header_check=setup.validate_headers,
    )

    col_ids = list(data.get("COLID", []))
    met_names = list(data.get("MET_NAME", []))
    file_names = list(data.get("FILE_NAME", []))
    metkm_values = list(data.get("STREAM_KM", []))
    metheight_values = list(data.get("MET_HEIGHT", []))

    if len(col_ids) != met_count:
        msg = "{0} must have exactly {1} data rows because metsites = {1} in the control file.".format(
            met_path.name, met_count
        )
        raise ValueError(msg)
    if any(col_id in (None, "") for col_id in col_ids):
        msg = "{0} is missing one or more COLID values.".format(met_path.name)
        raise ValueError(msg)
    if sorted(col_ids) != list(range(1, met_count + 1)):
        msg = "{0} must use COLID values 1 through {1}.".format(met_path.name, met_count)
        raise ValueError(msg)

    rows = list(zip(col_ids, met_names, file_names, metkm_values, metheight_values))
    rows.sort(key=lambda row: row[0])
    for col_id, met_name, file_name, stream_km, met_height in rows:
        met_rows.append(
            {
                "colid": col_id,
                "metname": met_name,
                "file_name": file_name,
                "stream_km": stream_km,
                "metheight": met_height,
            }
        )

    missing_file_names = any(file_name in (None, "") for file_name in file_names)
    missing_stream_km = any(stream_km in (None, "") for stream_km in metkm_values)
    missing_metheight = any(met_height in (None, "") for met_height in metheight_values)

    if missing_file_names:
        met_params["metfiles"] = None
    else:
        _validate_site_file_names(met_path.name, file_names)
        metfiles = []
        for row in met_rows:
            if row["file_name"] not in metfiles:
                metfiles.append(row["file_name"])
        met_params["metfiles"] = ", ".join(metfiles)

    if missing_stream_km:
        met_params["metkm"] = None
    else:
        met_params["metkm"] = ", ".join([str(row["stream_km"]) for row in met_rows])

    if missing_metheight:
        met_params["metheights"] = None
    else:
        met_params["metheights"] = ", ".join([str(row["metheight"]) for row in met_rows])

    result = (met_rows, met_params)
    return result


def _get_trib_sites(model_path, control_params, control_path, run_type, ext):
    trib_rows = []
    trib_params = {}
    trib_count = int(control_params.get("tribsites") or 0)

    if trib_count <= 0:
        result = (trib_rows, trib_params)
        return result

    trib_path = model_path / control_params["tribsitesfile"]
    if not trib_path.exists():
        result = (trib_rows, trib_params)
        return result

    setup = InputSetup(control_path)
    headers = headers_trib_sites()
    data = read_to_dict(
        path=trib_path,
        colnames=headers,
        sheetname=sheetnames["tribsitesfile"],
        value_check=setup._validate,
        header_check=setup.validate_headers,
    )

    col_ids = list(data.get("COLID", []))
    trib_names = list(data.get("TRIB_NAME", []))
    file_names = list(data.get("FILE_NAME", []))
    tribkm_values = list(data.get("STREAM_KM", []))

    if len(col_ids) != trib_count:
        msg = "{0} must have exactly {1} data rows because tribsites = {1} in the control file.".format(
            trib_path.name, trib_count
        )
        raise ValueError(msg)
    if any(col_id in (None, "") for col_id in col_ids):
        msg = "{0} is missing one or more COLID values.".format(trib_path.name)
        raise ValueError(msg)
    if sorted(col_ids) != list(range(1, trib_count + 1)):
        msg = "{0} must use COLID values 1 through {1}.".format(trib_path.name, trib_count)
        raise ValueError(msg)

    rows = list(zip(col_ids, trib_names, file_names, tribkm_values))
    rows.sort(key=lambda row: row[0])
    for col_id, trib_name, file_name, stream_km in rows:
        trib_rows.append(
            {
                "colid": col_id,
                "tribname": trib_name,
                "file_name": file_name,
                "stream_km": stream_km,
            }
        )

    missing_file_names = any(file_name in (None, "") for file_name in file_names)
    missing_stream_km = any(stream_km in (None, "") for stream_km in tribkm_values)

    if missing_file_names:
        trib_params["tribfiles"] = None
    else:
        _validate_site_file_names(trib_path.name, file_names)
        tribfiles = []
        for row in trib_rows:
            if row["file_name"] not in tribfiles:
                tribfiles.append(row["file_name"])
        trib_params["tribfiles"] = ", ".join(tribfiles)

    if missing_stream_km:
        trib_params["tribkm"] = None
    else:
        trib_params["tribkm"] = ", ".join([str(row["stream_km"]) for row in trib_rows])

    result = (trib_rows, trib_params)
    return result


def get_site_files(control_path, control_params, run_type):
    control_path = Path(control_path).expanduser().resolve()
    model_path = control_path.parent
    ext = control_path.suffix.lower()

    met_rows, met_params = _get_met_sites(
        model_path,
        control_params,
        control_path,
        run_type,
        ext,
    )
    trib_rows, trib_params = _get_trib_sites(
        model_path,
        control_params,
        control_path,
        run_type,
        ext,
    )

    result = {
        "met_rows": met_rows,
        "met_params": met_params,
        "trib_rows": trib_rows,
        "trib_params": trib_params,
    }
    return result

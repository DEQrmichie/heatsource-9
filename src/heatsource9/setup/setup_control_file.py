from datetime import datetime
import logging
from pathlib import Path

from heatsource9.io.control_file import import_control_file
from heatsource9.io.input_files import write_input
from heatsource9.io.logging_config import configure_logging
from heatsource9.setup.constants import dtype, sheetnames

logger = logging.getLogger(__name__)


def control_file_dict():
    """
    Returns a dictionary of the control file lines
    with empty values. Dictionary keys correspond to the keys used
    in the control file.

    """
    cf_dict = {
        "usertxt": [2, "Model Description/User Notes", "usertxt", None],
        "name": [3, "Simulation Name", "name", None],
        "inputdir": [4, "Input Directory Path", "inputdir", None],
        "outputdir": [5, "Output Directory Path", "outputdir", None],
        "length": [6, "Stream Length (kilometers)", "length", None],
        "outputkm": [7, "Output Stream Kilometers", "outputkm", None],
        "datastart": [8, "Data Start Date (yyyy-mm-dd)", "datastart", None],
        "modelstart": [9, "Modeling Start Date (yyyy-mm-dd)", "modelstart", None],
        "modelend": [10, "Modeling End Date (yyyy-mm-dd)", "modelend", None],
        "dataend": [11, "Data End Date (yyyy-mm-dd)", "dataend", None],
        "flushdays": [12, "Flush Initial Condition (days)", "flushdays", None],
        "offset": [13, "Time Offset From UTC (hours)", "offset", None],
        "dt": [14, "Model Time Step (minutes)", "dt", None],
        "dx": [15, "Model Distance Step (meters)", "dx", None],
        "outputdt": [16, "Output Time Step (minutes)", "outputdt", None],
        "longsample": [17, "Longitudinal Stream Sample Distance (meters)", "longsample", None],
        "bcfile": [18, "Boundary Condition Input File Name", "bcfile", None],
        "tribsites": [19, "Tributary Inflow Sites", "tribsites", None],
        "accretionfile": [20, "Accretion Input File Name", "accretionfile", None],
        "metsites": [21, "Meteorological Data Sites", "metsites", None],
        "calcevap": [22, "Include Evaporation Losses From Flow (True/False)", "calcevap", None],
        "evapmethod": [23, "Evaporation Method (Mass Transfer/Penman)", "evapmethod", None],
        "wind_a": [24, "Wind Function Coefficient a", "wind_a", None],
        "wind_b": [25, "Wind Function Coefficient b", "wind_b", None],
        "calcalluvium": [26, "Include Deep Alluvium Temperature (True/False)", "calcalluvium", None],
        "alluviumtemp": [27, "Deep Alluvium Temperature (Celsius)", "alluviumtemp", None],
        "morphfile": [28, "Morphology Input Data File Name", "morphfile", None],
        "lcdatafile": [29, "Land Cover Input Data File Name", "lcdatafile", None],
        "lccodefile": [30, "Land Cover Codes Input File Name", "lccodefile", None],
        "trans_count": [31, "Number Of Transects Per Node", "trans_count", None],
        "transsample_count": [32, "Number Of Samples Per Transect", "transsample_count", None],
        "transsample_distance": [33, "Distance Between Transect Samples (meters)", "transsample_distance", None],
        "emergent": [34, "Account For Emergent Veg Shading (True/False)", "emergent", None],
        "lcdatainput": [35, "Land Cover Data Input Type (Codes/Values)", "lcdatainput", None],
        "canopy_data": [36, "Canopy Data Type (LAI/CanopyCover)", "canopy_data", None],
        "lcsampmethod": [37, "Land Cover Sample Method (point/zone)", "lcsampmethod", None],
        "heatsource8": [38, "Use Heat Source 8 Land Cover Methods (True/False)", "heatsource8", None],
    }

    return cf_dict


def write_cf(
    model_dir,
    control_file = None,
    *,
    use_timestamp = False,
    overwrite = False,
    csv_mode = False,
):
    model_path = Path(model_dir)
    if control_file and Path(control_file).is_absolute():
        msg = "control_file must be a file name, not a full path."
        raise ValueError(msg)
    filename = control_file or (
        "HeatSource_Control.csv" if csv_mode else "HeatSource_Control.xlsx"
    )
    if use_timestamp:
        target_name = datetime.now().strftime("%Y-%m-%d_%H%M%S") + "_" + Path(filename).name
    else:
        target_name = Path(filename).name
    target = model_path / target_name

    if target.exists() and not overwrite:
        return target

    rows = [row for _, row in sorted(control_file_dict().items(), key=lambda kv: kv[1][0])]
    write_input(
        path=target,
        headers=["LINE", "PARAMETER", "KEY", "VALUE"],
        rows=rows,
        sheetname=sheetnames["controlfile"],
        csv_mode=csv_mode,
    )
    return target


def _format_control_value(key, value):
    if value is None:
        return ""

    if key not in dtype:
        raise KeyError(f"Unknown control key: {key}")

    expected = dtype[key]
    if expected is bool:
        if isinstance(value, bool):
            return "True" if value else "False"
        text = str(value).strip().lower()
        if text in {"true", "1", "yes", "y"}:
            return "True"
        if text in {"false", "0", "no", "n"}:
            return "False"
        raise TypeError(f"Control key '{key}' expects a True or False value.")
    if expected is int:
        return str(int(float(value)))
    if expected is float:
        return str(float(value))
    if expected == "date":
        text = str(value).strip()
        datetime.strptime(text, "%Y-%m-%d")
        return text
    if expected is str:
        return str(value).strip()
    return str(value).strip()


def setup_cf(
    model_dir,
    control_file = "HeatSource_Control.xlsx",
    *,
    use_timestamp = False,
    overwrite = False,
    strict = True,
    **control_values,
):
    model_path = Path(model_dir).expanduser().resolve()
    configure_logging(model_dir=model_path, overwrite=True)
    try:
        control_path = Path(control_file).expanduser()
        if not control_path.is_absolute():
            control_path = model_path / control_path
        control_path = control_path.resolve()

        ext = control_path.suffix.lower()
        if ext not in [".xlsx", ".csv"]:
            msg = "{0} must be an Excel '.xlsx' or '.csv' file.".format(control_path.name)
            raise ValueError(msg)

        csv_mode = ext == ".csv"

        if not control_values:
            return write_cf(
                model_dir=control_path.parent,
                control_file=control_path.name,
                use_timestamp=use_timestamp,
                overwrite=overwrite,
                csv_mode=csv_mode,
            )

        if not control_path.exists():
            write_cf(
                model_dir=control_path.parent,
                control_file=control_path.name,
                use_timestamp=False,
                overwrite=False,
                csv_mode=csv_mode,
            )
        elif not overwrite:
            return control_path
        else:
            write_cf(
                model_dir=control_path.parent,
                control_file=control_path.name,
                use_timestamp=False,
                overwrite=True,
                csv_mode=csv_mode,
            )

        raw_rows = import_control_file(
            control_path=control_path,
            dtype=dtype,
            control_sheet=sheetnames["controlfile"],
        )["control_rows"]
        base_rows = {k: list(v) for k, v in control_file_dict().items()}
        valid_control_keys = set(base_rows.keys())

        for row in raw_rows:
            key = str(row["key"]).strip()
            if key in base_rows:
                base_rows[key][3] = row["value"]

        for key, value in control_values.items():
            if strict and key not in dtype:
                raise KeyError(f"Unknown dtype for control key: {key}")
            if strict and key not in valid_control_keys and key not in dtype:
                raise KeyError(f"Unknown control key: {key}")
            if key in valid_control_keys:
                base_rows[key][3] = _format_control_value(key, value)

        if use_timestamp:
            output_name = datetime.now().strftime("%Y-%m-%d_%H%M%S") + "_" + control_path.name
        else:
            output_name = control_path.name
        target = model_path / output_name
        if target.exists() and target != control_path and not overwrite:
            return target

        rows = [row for _, row in sorted(base_rows.items(), key=lambda kv: kv[1][0])]
        write_input(
            path=target,
            headers=["LINE", "PARAMETER", "KEY", "VALUE"],
            rows=rows,
            sheetname=sheetnames["controlfile"],
            csv_mode=csv_mode,
        )
        return target
    except Exception as exc:
        msg = str(exc)
        logger.exception(msg)
        raise

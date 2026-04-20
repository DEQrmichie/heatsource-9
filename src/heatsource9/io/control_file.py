import logging
from pathlib import Path
from datetime import datetime


logger = logging.getLogger(__name__)


def validate_control_file(control_file_path):
    """
    Validate the control file path, check for .csv or .xlsx, 
    and return the path string.
    """
    p = Path(control_file_path).expanduser()
    if not p.exists():
        msg = f"Control file not found: {str(p)}"
        raise FileNotFoundError(msg)

    ext = p.suffix.lower()
    if ext not in [".xlsx", ".csv"]:
        msg = "{0} must be an Excel '.xlsx' or '.csv' file.".format(p.name)
        raise ValueError(msg)

    return str(p.resolve())


def cf_path(model_dir, control_file = None):
    """
    Return the full path to the control file.

    If control_file is provided, absolute paths are used directly and
    relative paths are resolved under model_dir.
    If control_file is not provided, this finds HeatSource_Control.* in model_dir.
    If there is not exactly one match, an error is raised.
    The control file must be an Excel '.xlsx' or '.csv' file.
    Can't have both formats in the same directory.
    """
    model_path = Path(model_dir).expanduser().resolve()

    if control_file:
        control_candidate = Path(control_file).expanduser()
        if control_candidate.is_absolute():
            control_path = control_candidate.resolve()
        else:
            control_path = (model_path / control_file).resolve()
    else:
        controls = sorted(model_path.glob("HeatSource_Control.*"))
        if len(controls) == 0:
            msg = (
                "HeatSource_Control file not found. Move the executable or place the control file "
                "in this directory: {0}.".format(model_path)
            )
            raise FileNotFoundError(msg)
        if len(controls) > 1:
            msg = (
                "There is more than one file named 'HeatSource_Control.' in this directory: {0}. "
                "Only one file can exist.".format(model_path)
            )
            raise FileExistsError(msg)
        control_path = controls[0].resolve()

    return Path(validate_control_file(control_path))


def _rename_control_key(key):
    from heatsource9.setup.constants import legacy_cf_format_keys, renamed_cf_keys

    if key in renamed_cf_keys:
        new_key = renamed_cf_keys[key]
        msg = 'Control key "{0}" is deprecated, use "{1}".'.format(key, new_key)
        logger.warning(msg)
        return new_key
    if key in legacy_cf_format_keys:
        new_key = legacy_cf_format_keys[key]
        msg = 'Control key "{0}" is deprecated, use the sites file format instead.'.format(key)
        logger.warning(msg)
        return new_key
    return key


def _to_value_type(key, value, dtype):
    """
    Just a helper to make sure the control file values are 
    the right data type.
    """
    expected = dtype.get(key)
    text = value.strip()
    if text == "":
        return None
    if expected is str:
        return text
    if expected is int:
        return int(float(text))
    if expected is float:
        return float(text)
    if expected is bool:
        title = text.title()
        if title in {"True", "False"}:
            return title == "True"
        msg = f"Control file key '{key}' must be True or False."
        raise TypeError(msg)
    if expected == "date":
        date_text = text.split()[0]
        datetime.strptime(date_text, "%Y-%m-%d")
        return date_text
    if expected == "datetime":
        return text
    return text


def import_control_file(control_path, *, dtype, control_sheet):
    """
    Reads the control file and returns a dictionary with keys:
    model_dir, control_file, control_rows, control_params.
    """
    # this import is here so it's only called when actually reading the control file, 
    # If at the top openpyxl would be imported every time when it's not needed.
    from heatsource9.io.input_files import read_to_list

    cf_list = read_to_list(control_path, skiprows=1, skipcols=0, sheetname=control_sheet)

    rows = []
    params = {}

    for line in cf_list:
        if not line or len(line) < 4:
            continue

        raw_key = str(line[2]).strip()
        key = _rename_control_key(raw_key)
        raw_val = line[3]
        if raw_val is None:
            val = ""
        elif isinstance(raw_val, datetime):
            if key in dtype and dtype[key] == "date":
                val = raw_val.strftime("%Y-%m-%d")
            else:
                val = raw_val.strftime("%Y-%m-%d %H:%M")
        else:
            val = str(raw_val).strip()
            if key in dtype and dtype[key] == "date" and val:
                val = val.split()[0]
                datetime.strptime(val, "%Y-%m-%d")

        rows.append(
            {
                "line": str(line[0]).strip(),
                "parameter": str(line[1]).strip(),
                "key": key,
                "value": val,
            }
        )

        if key:
            if raw_key != key and key in params:
                continue
            params[key] = _to_value_type(key, val, dtype)

    ext = control_path.suffix.lower()
    if ext == ".csv":
        params["metsitesfile"] = params.get("metsitesfile") or "HeatSource_Met_Sites.csv"
        params["tribsitesfile"] = params.get("tribsitesfile") or "HeatSource_Tributary_Sites.csv"
    else:
        params["metsitesfile"] = params.get("metsitesfile") or "HeatSource_Met_Sites.xlsx"
        params["tribsitesfile"] = params.get("tribsitesfile") or "HeatSource_Tributary_Sites.xlsx"

    return {
        "model_dir": str(control_path.parent),
        "control_file": control_path.name,
        "control_rows": rows,
        "control_params": params,
    }

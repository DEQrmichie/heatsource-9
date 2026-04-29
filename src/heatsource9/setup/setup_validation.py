"""Validation functions to support setup and model runs."""
import logging
from decimal import Decimal
from decimal import InvalidOperation

from heatsource9.setup.constants import KM_PRECISION, required_fields, required_files

logger = logging.getLogger(__name__)


def is_blank(value):
    """Return True when the value is None or an empty string."""
    if value is None:
        return True
    if isinstance(value, str) and value.strip() == "":
        return True
    return False


def file_required(run_type, file_key, params = None, source = None):
    """
    Return True if the input file is required for the selected run type.
    """
    params = params or {}
    required = required_files.get(run_type, ())
    if file_key not in required:
        return False

    if file_key == "tribsitesfile":
        return int(params.get("tribsites") or 0) > 0

    if file_key == "metsitesfile":
        return run_type in ("solar", "temperature")

    return True


def conditional_required(field_name, params = None, source = None):
    """
    Return True when a conditional control file key is required
    based on the configuration of other control file keys.
    """
    params = params or {}

    if field_name == "alluviumtemp":
        return params.get("calcalluvium") is True

    if field_name in ("tribfiles", "tribkm"):
        if int(params.get("tribsites") or 0) <= 0:
            return False
        return is_blank(params.get("tribsitesfile"))

    if field_name in ("metfiles", "metkm"):
        run_type = params.get("run_type")
        if run_type in ("solar", "temperature"):
            return is_blank(params.get("metsitesfile"))
        return False

    if field_name == "metheights":
        run_type = params.get("run_type")
        if run_type == "temperature":
            return is_blank(params.get("metsitesfile"))
        return False

    return False


def field_required(run_type, file_key, field_name, params = None, source = None):
    """
    Return True when a field must have a value for this run type and file.
    """
    params = params or {}
    params = dict(params)
    params["run_type"] = run_type

    required = required_fields.get(run_type, {}).get(file_key, ())

    if file_key == "lcdatafile" and run_type in ("solar", "temperature"):
        for p in ("LC", "ELE"):
            if field_name.startswith(p + "_"):
                return True
        return field_name in required

    if field_name not in required:
        return False

    if file_key == "controlfile":
        if field_name in ("tribfiles", "tribkm", "metfiles", "metkm", "metheights", "alluviumtemp"):
            return conditional_required(field_name, params=params, source=source)

    if file_key == "metsitesfile" and field_name == "MET_HEIGHT":
        if run_type == "temperature":
            return True
        if run_type == "solar":
            return False

    return True


def set_default(field_name, run_type = None, file_key = None, params = None, source = None):
    """
    Returns a default value for certain blank control file keys or fields.
    """
    params = params or {}

    if field_name == "outputdt":
        return 60.0

    if field_name == "metheights":
        if run_type == "temperature" and file_key == "controlfile" and is_blank(params.get("metsitesfile")):
            return 2.0

    if file_key == "bcfile" and field_name == "TEMPERATURE":
        return 0.0

    if file_key == "metfiles" and run_type == "solar":
        for prefix in ("WIND_SPEED", "RELATIVE_HUMIDITY", "AIR_TEMPERATURE"):
            if field_name.startswith(prefix):
                return 0.0

    if file_key == "accretionfile" and run_type in ("hydraulics", "temperature"):
        if field_name in ("INFLOW", "TEMPERATURE", "OUTFLOW"):
            return 0.0

    if file_key == "morphfile" and run_type in ("hydraulics", "temperature"):
        if field_name == "HYPORHEIC_PERCENT":
            return 0.0

    return None


def validate_required_field(run_type, file_key, field_name, value, params = None, source = None):
    """
    Validate one field value using required checks and confirmed defaults.

    Returns the existing value when present, raises ValueError for required blanks,
    otherwise returns a default or None.
    """
    params = params or {}
    params = dict(params)
    params["run_type"] = run_type

    if not is_blank(value):
        return value

    default = set_default(
        field_name=field_name,
        run_type=run_type,
        file_key=file_key,
        params=params,
        source=source,
    )
    if default is not None:
        if field_name == "outputdt":
            logger.warning("Control file key 'outputdt' is missing. Defaulting to 60 minutes.")
        return default

    required = field_required(
        run_type=run_type,
        file_key=file_key,
        field_name=field_name,
        params=params,
        source=source,
    )
    if required:
        msg = "Value for required field '{0}' in '{1}' is missing.".format(field_name, file_key)
        raise ValueError(msg)

    return None


def validate_column_values(run_type, file_key, data_dict, params = None, source = None):
    """
    Validate all values in a column based data dictionary.

    Applies validate_required_field to each value and returns the updated dictionary.
    """
    if run_type is None:
        return data_dict
    if not data_dict:
        return data_dict

    params = params or {}
    for field_name, values in list(data_dict.items()):
        data_dict[field_name] = [
            validate_required_field(
                run_type=run_type,
                file_key=file_key,
                field_name=field_name,
                value=value,
                params=params,
                source=source,
            ) for value in values
        ]
    return data_dict


def align_rows_to_kmlist(
    file_key,
    data_dict,
    kmlist,
    file_name,
    longsample,
):
    """
    Align the stream km from the input file to the internal model kmlist.
    This handles unsorted input rows and assigns values to the correct model km.

    Compare each input stream km to the nearest expected model km from kmlist, calculated
    using the longsample and length input into the control file.
    If the difference is greater than 0.001 km (1 meter), it is not a match and an error is raised.
    If the difference is <= 0.001 km, it is treated as a match. This allows small float noise in Excel spreadsheets, for example 3.949999999999 instead of 3.95.
    """
    row_count = len(data_dict["STREAM_KM"])
    step_km = Decimal(str(longsample)) / Decimal("1000")

    expected_index = {
        format(Decimal(str(km)).quantize(KM_PRECISION), "f"): i for i, km in enumerate(kmlist)
    }

    km_aligned = {name: [None for _ in kmlist] for name in data_dict}
    km_found = set()
    km_expected = set()
    km_nomatch = []
    for row_idx in range(row_count):
        raw_km = data_dict["STREAM_KM"][row_idx]
        try:
            raw_km_dec = Decimal(str(raw_km))
        except (InvalidOperation, ValueError, TypeError):
            msg = "Value {0!r} for STREAM_KM in '{1}' is not a valid number.".format(raw_km, file_name)
            raise ValueError(msg)
        snapped_km = (raw_km_dec / step_km).quantize(Decimal("1")) * step_km
        if abs(raw_km_dec - snapped_km) > KM_PRECISION:
            msg = (
                "STREAM_KM {0} in {1} is not a multiple of the longsample ({2} m or {3} km), expected {4}."
            ).format(
                raw_km, file_name, longsample, format(step_km, "f"), format(snapped_km.quantize(Decimal("0.0001")), "f")
            )
            raise ValueError(msg)
        key = format(snapped_km.quantize(KM_PRECISION), "f")
        if key in km_found:
            msg = "Input '{0}' has duplicate STREAM_KM value {1}.".format(file_key, raw_km)
            raise ValueError(msg)
        km_found.add(key)
        if key not in expected_index:
            km_nomatch.append(key)
            continue
        target_index = expected_index[key]
        km_expected.add(target_index)
        for col_name, values in list(data_dict.items()):
            km_aligned[col_name][target_index] = values[row_idx]

    missing = [kmlist[i] for i in range(len(kmlist)) if i not in km_expected]
    if missing:
        first_values = ", ".join([format(km, ".6f") for km in missing[:5]])
        msg = "Input '{0}' is missing rows for model STREAM_KM values: {1}.".format(file_key, first_values)
        raise ValueError(msg)

    if km_nomatch:
        km_nomatch = list(dict.fromkeys(km_nomatch))
        extra_values = ", ".join(km_nomatch)
        msg = (
            "The {0} file named {1} has extra STREAM_KM rows not used by the model. "
            "Extra stream kilometers include: {2}"
        ).format(file_key, file_name, extra_values)
        logger.warning(msg)

    return km_aligned

"""Convert the raw input values into bool, int, float, and string type values."""


def as_bool(value, default = False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text == "true":
        return True
    if text == "false":
        return False
    raise ValueError(f"Invalid boolean value {value!r}; expected 'true' or 'false'.")


def as_int(value, default = 0):
    if value is None or str(value).strip() == "":
        return default
    return int(float(str(value).strip()))


def as_float(value, default = 0.0):
    if value is None or str(value).strip() == "":
        return default
    return float(str(value).strip())


def as_str(value, default = ""):
    if value is None:
        return default
    return str(value).strip()

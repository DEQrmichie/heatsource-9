from datetime import datetime, timedelta, timezone
from math import ceil
from pathlib import Path

from openpyxl.utils import datetime as pyxl_datetime

from heatsource9.io.control_file import cf_path, import_control_file
from heatsource9.io.input_files import write_input
from heatsource9.setup.constants import dtype, sheetnames
from heatsource9.setup.headers import (
    headers_accretion,
    headers_bc,
    headers_inflow,
    headers_lccodes,
    headers_lcdata,
    headers_met,
    headers_morph,
)


def _parse_date_yyyy_mm_dd(value):
    return datetime.strptime(value.strip(), "%Y-%m-%d").replace(tzinfo=timezone.utc)


def _hourly_range(start_date, end_date, as_excel):
    """"
    Hourly range is inclusive through end date + 24h.
    """

    start_dt = _parse_date_yyyy_mm_dd(start_date)
    end_dt = _parse_date_yyyy_mm_dd(end_date) + timedelta(days=1)

    values = []
    current = start_dt
    while current <= end_dt:
        if as_excel:
            values.append(
                pyxl_datetime.to_excel(
                    dt=current,
                    epoch=datetime(1899, 12, 30, 0, 0, tzinfo=timezone.utc),
                )
            )
        else:
            values.append(current.strftime("%Y-%m-%d %H:%M"))
        current += timedelta(hours=1)
    return values


def _compute_stream_kms(length_km, longsample_m):
    """
    Return a list of stream kilometers sorted from upstream to downstream.
    """
    num_nodes = int(ceil(round(length_km * 1000 / longsample_m, 4))) + 1
    kms = [(float(node) * longsample_m) / 1000 for node in range(0, num_nodes)]
    kms.sort(reverse=True)
    return kms


def write_mi(
    model_dir,
    control_file = None,
    *,
    use_timestamp = False,
    overwrite = False,
):
    model_path = Path(model_dir).expanduser().resolve()
    control_path = cf_path(model_path, control_file)

    params = import_control_file(
        control_path=control_path,
        dtype=dtype,
        control_sheet=sheetnames["controlfile"],
    )["control_params"]

    inputdir = Path(str(params.get("inputdir") or (model_path / "Inputs"))).expanduser().resolve()
    outputdir = Path(str(params.get("outputdir") or (model_path / "Output"))).expanduser().resolve()
    for p in (inputdir, outputdir):
        try:
            p.mkdir(parents=True, exist_ok=True)
        except PermissionError as exc:
            raise PermissionError(
                f"Cannot create folder '{p}'. Check write permissions."
            ) from exc

    length_km = params.get("length")
    longsample = params.get("longsample")
    if length_km is None or longsample is None or length_km <= 0 or longsample <= 0:
        raise ValueError("Control file must define 'length' and 'longsample' before model input setup")

    datastart = params.get("datastart")
    dataend = params.get("dataend")
    if not datastart or not dataend:
        raise ValueError("Control file must define 'datastart' and 'dataend' before model input setup")

    csv_mode = control_path.suffix.lower() == ".csv"
    timelist = _hourly_range(datastart, dataend, as_excel=not csv_mode)
    kmlist = _compute_stream_kms(length_km, longsample)

    def out_name(base):
        base_name = str(params.get(base) or (base + ".csv"))
        if use_timestamp:
            return datetime.now().strftime("%Y-%m-%d_%H%M%S") + "_" + base_name
        return base_name

    written = []

    acc = [[None, None, km, None, None, None] for km in kmlist]
    bc = [[t, None, None] for t in timelist]
    morph = [[None, None, km] + [None] * 10 for km in kmlist]
    lc_headers = headers_lcdata(params)
    lcdata = [[None, None, km] + [None] * (len(lc_headers) - 3) for km in kmlist]

    tables = [
        (out_name("accretionfile"), headers_accretion(), acc, sheetnames["accretionfile"]),
        (out_name("bcfile"), headers_bc(), bc, sheetnames["bcfile"]),
        (out_name("morphfile"), headers_morph(), morph, sheetnames["morphfile"]),
        (out_name("lcdatafile"), lc_headers, lcdata, sheetnames["lcdatafile"]),
    ]

    lc_code_headers = headers_lccodes(params)
    if lc_code_headers != [None]:
        tables.append((out_name("lccodefile"), lc_code_headers, [[None]], sheetnames["lccodefile"]))

    met_files = [p.strip() for p in str(params.get("metfiles") or "met.csv").split(",") if p.strip()]
    met_headers = headers_met(params)
    met_rows = [[t] + [None] * (len(met_headers) - 1) for t in timelist]
    for f in met_files:
        if use_timestamp:
            met_name = datetime.now().strftime("%Y-%m-%d_%H%M%S") + "_" + f
        else:
            met_name = f
        tables.append((met_name, met_headers, met_rows, sheetnames["metfiles"]))

    inflow_sites = params.get("inflowsites") or 0
    if inflow_sites > 0:
        inflow_files_value = params.get("inflowinfiles")
        inflow_files = [p.strip() for p in str(inflow_files_value or "").split(",") if p.strip()]
        inflow_headers = headers_inflow(params)
        inflow_rows = [[t] + [None] * (len(inflow_headers) - 1) for t in timelist]
        for f in inflow_files:
            if use_timestamp:
                inflow_name = datetime.now().strftime("%Y-%m-%d_%H%M%S") + "_" + f
            else:
                inflow_name = f
            tables.append((inflow_name, inflow_headers, inflow_rows, sheetnames["inflowinfiles"]))

    for filename, headers, rows, sheet in tables:
        path = inputdir / filename
        if path.exists() and not overwrite:
            written.append(path)
            continue
        write_input(path=path, headers=headers, rows=rows, sheetname=sheet, csv_mode=csv_mode)
        written.append(path)

    return written

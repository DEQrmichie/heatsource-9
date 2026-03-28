

def headers_accretion():
    return ["STREAM_ID", "NODE_ID", "STREAM_KM", "INFLOW", "TEMPERATURE", "OUTFLOW"]


def headers_bc():
    return ["DATETIME", "FLOW", "TEMPERATURE"]


def headers_met(params):
    metsites = params.get("metsites") or 0
    metfiles = [p.strip() for p in str(params.get("metfiles") or "met.csv").split(",") if p.strip()]
    ncols = int(metsites // max(1, len(metfiles)))
    header = ["DATETIME"]
    for n in range(1, ncols + 1):
        header += [
            f"CLOUDINESS{n}",
            f"WIND_SPEED{n}",
            f"RELATIVE_HUMIDITY{n}",
            f"AIR_TEMPERATURE{n}",
        ]
    return header


def headers_inflow(params):
    inflow_sites = params.get("inflowsites") or 0
    inflow_files_value = params.get("inflowinfiles")
    inflow_files = [p.strip() for p in str(inflow_files_value or "").split(",") if p.strip()]
    if inflow_sites <= 0:
        return [None]
    ncols = int(inflow_sites // max(1, len(inflow_files)))
    header = ["DATETIME"]
    for n in range(1, ncols + 1):
        header += [f"FLOW{n}", f"TEMPERATURE{n}"]
    return header


def headers_lccodes(params):
    if (params.get("lcdatainput") or "Codes") != "Codes":
        return [None]
    if (params.get("canopy_data") or "LAI") == "LAI":
        return ["NAME", "CODE", "HEIGHT", "LAI", "k", "OVERHANG", "CANOPY_DEPTH"]
    return ["NAME", "CODE", "HEIGHT", "CANOPY", "OVERHANG", "CANOPY_DEPTH"]


def headers_lcdata(params):
    lcdatainput = params.get("lcdatainput") or "Codes"
    canopy_data = params.get("canopy_data") or "LAI"
    if lcdatainput == "Values":
        prefix = ["HT", "ELE", "LAI", "k", "OH", "CD"] if canopy_data == "LAI" else ["HT", "ELE", "CAN", "OH", "CD"]
    else:
        prefix = ["LC", "ELE"]

    trans_count = params.get("trans_count") or 0
    transsample_count = params.get("transsample_count") or 0
    heatsource8 = params.get("heatsource8") or False

    headers = [
        "STREAM_ID",
        "NODE_ID",
        "STREAM_KM",
        "LONGITUDE",
        "LATITUDE",
        "TOPO_W",
        "TOPO_S",
        "TOPO_E",
    ]

    tran = ["NE", "E", "SE", "S", "SW", "W", "NW"] if heatsource8 else [f"T{x}" for x in range(1, trans_count + 1)]
    zones = list(range(1, transsample_count + 1))

    for p in prefix:
        for ti, t in enumerate(tran):
            for z in zones:
                if p != "ELE" and ti == 0 and z == 1:
                    headers.append(f"{p}_T0_S0")
                headers.append(f"{p}_{t}_S{z}")
    return headers


def headers_morph():
    return [
        "STREAM_ID",
        "NODE_ID",
        "STREAM_KM",
        "ELEVATION",
        "GRADIENT",
        "BOTTOM_WIDTH",
        "CHANNEL_ANGLE_Z",
        "MANNINGS_n",
        "SED_THERMAL_CONDUCTIVITY",
        "SED_THERMAL_DIFFUSIVITY",
        "SED_HYPORHEIC_THICKNESS",
        "HYPORHEIC_PERCENT",
        "POROSITY",
    ]

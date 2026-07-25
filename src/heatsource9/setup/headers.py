

def headers_accretion():
    return ["STREAM_ID", "NODE_ID", "STREAM_KM", "INFLOW", "TEMPERATURE", "OUTFLOW"]


def headers_bc():
    return ["DATETIME", "FLOW", "TEMPERATURE"]


def headers_met_sites():
    return ["COLID", "MET_NAME", "STREAM_KM", "FILE_NAME", "MET_HEIGHT"]


def headers_trib_sites():
    return ["COLID", "TRIB_NAME", "STREAM_KM", "FILE_NAME"]


def headers_met(metsites, metfile_count):
    ncols = int(metsites // metfile_count)
    header = ["DATETIME"]
    for n in range(1, ncols + 1):
        header += [
            f"CLOUDINESS{n}",
            f"WIND_SPEED{n}",
            f"RELATIVE_HUMIDITY{n}",
            f"AIR_TEMPERATURE{n}",
        ]
    return header


def headers_tribfiles(tribsites, tribfile_count):
    ncols = int(tribsites // tribfile_count)
    header = ["DATETIME"]
    for n in range(1, ncols + 1):
        header += [f"FLOW{n}", f"TEMPERATURE{n}"]
    return header


def headers_lccodes(canopy_data):
    if canopy_data == "LAI":
        return ["NAME", "CODE", "HEIGHT", "LAI", "k", "OVERHANG", "CANOPY_DEPTH"]
    return ["NAME", "CODE", "HEIGHT", "CANOPY", "OVERHANG", "CANOPY_DEPTH"]


def headers_lcdata(trans_count, transsample_count, heatsource8):
    prefix = ["LC", "ELE"]

    headers = [
        "STREAM_ID",
        "NODE_ID",
        "STREAM_KM",
        "LONGITUDE",
        "LATITUDE",
        "TOPO_NE",
        "TOPO_E",
        "TOPO_SE",
        "TOPO_S",
        "TOPO_SW",
        "TOPO_W",
        "TOPO_NW",
        "TOPO_N",
    ]

    if heatsource8:
        tran = ["NE", "E", "SE", "S", "SW", "W", "NW"]
    else:
        tran = [f"T{x}" for x in range(1, trans_count + 1)]
    samples = list(range(1, transsample_count + 1))

    for p in prefix:
        for ti, t in enumerate(tran):
            for s in samples:
                if p != "ELE" and ti == 0 and s == 1:
                    headers.append(f"{p}_T0_S0")
                headers.append(f"{p}_{t}_S{s}")
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

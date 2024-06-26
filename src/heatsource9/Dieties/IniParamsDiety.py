"""
IniParams holds regular dictionaries with information
about the model initialization parameters and model inputs.
It is located here so be they can be accessed globally.
"""

# IniParms mostly holds the control file info.
# A couple initial values are suggested to make life easy.
IniParams = {"run_in_python": True,
             "run_type": None,
             "usertxt": None,
             "name": None,
             "inputdir": None,
             "outputdir": None,
             "length": None,
             "outputkm": "all",
             "datastart": None,
             "modelstart": None,
             "modelend": None,
             "dataend": None,
             "flushdays": None,
             "offset": None,
             "dt": 1,
             "dx": 50,
             "longsample": None,
             "bcfile": "bc.csv",
             "inflowsites": None,
             "inflowinfiles": None,
             "inflowkm": None,
             "accretionfile": "accretion.csv",
             "metsites": None,
             "metfiles": "met.csv",
             "metkm": None,
             "calcevap": None,
             "evapmethod": "Mass Transfer",
             "wind_a": 0.00000000151,
             "wind_b": 0.0000000016,
             "calcalluvium": None,
             "alluviumtemp": None,
             "morphfile": "morph.csv",
             "lcdatafile": "lcdata.csv",
             "lccodefile": "lccodes.csv",
             "trans_count": None,
             "transsample_count": None,
             "transsample_distance": None,
             "emergent": "False",
             "lcdatainput": None,
             "canopy_data": None,
             "lcsampmethod": "point",
             "heatsource8": "False", }

# name of excel sheets

sheetnames = {"accretionfile": "Accretion Flow",
              "bcfile": "Boundary Conditions",
              "controlfile": "Control Settings",
              "inflowinfiles": "Tributary Data",
              "metfiles": "Meteorological Data",
              "morphfile": "Morphology Data",
              "lccodefile": "Land Cover Codes",
              "lcdatafile": "Land Cover Data"}

# dype is a dictionary holding the data type for every input. Includes 
# model variable names and input file names. Note the datetime inputs
# are identified as strings here but are converted to integers
# after import. And not to be confusing but the variables identified as 
# bool type here are str in the input files and converted to bool
# upon import. See Inputs class validate() and import_control_file() 
# for details.

dtype = {"run_in_python": bool,
         "usertxt": str,
         "name": str,
         "inputdir": str,
         "outputdir": str,
         "length": float,
         "outputkm": str,
         "datastart": "datetime",
         "modelstart": "datetime",
         "modelend": "datetime",
         "dataend": "datetime",
         "flushdays": int,
         "offset": int,
         "dt": float,
         "dx": int,
         "longsample": int,
         "bcfile": str,
         "inflowsites": int,
         "inflowinfiles": str,
         "inflowkm": str,
         "accretionfile": str,
         "metsites": int,
         "metfiles": str,
         "metkm": str,
         "calcevap": bool,
         "evapmethod": str,
         "wind_a": float,
         "wind_b": float,
         "calcalluvium": bool,
         "alluviumtemp": float,
         "morphfile": str,
         "lcdatafile": str,
         "lccodefile": str,
         "trans_count": int,
         "transsample_count": int,
         "transsample_distance": float,
         "emergent": bool,
         "lcdatainput": str,
         "canopy_data": str,
         "lcsampmethod": str,
         "heatsource8": bool,
         "penman": bool,
         "nodeID": int,
         "streamID": str,
         "km": float,
         "longitude": float,
         "latitude": float,
         "topo": float,
         "lc_code": str,
         "lc_height": float,
         "lc_canopy": float,
         "lc_lai": float,
         "lc_k": float,
         "lc_oh": float,
         "elevation": float,
         "S": float,
         "W_b": float,
         "z": float,
         "n": float,
         "SedThermCond": float,
         "SedThermDiff": float,
         "SedDepth": float,
         "hyp_percent": float,
         "phi": float,
         "Q_cont": float,
         "d_cont": float,
         "Q_in": float,
         "T_in": float,
         "Q_out": float,
         "cloud": float,
         "wind": float,
         "humidity": float,
         "T_air": float,
         "STREAM_ID": str,
         "NODE_ID": int,
         "STREAM_KM": float,
         "INFLOW": float,
         "TEMPERATURE": float,
         "OUTFLOW": float,
         "DATETIME": "datetime",
         "CLOUDINESS": float,
         "WIND_SPEED": float,
         "RELATIVE_HUMIDITY": float,
         "AIR_TEMPERATURE": float,
         "FLOW": float,
         "NAME": str,
         "CODE": str,
         "HEIGHT": float,
         "CANOPY": float,
         "OVERHANG": float,
         "LONGITUDE": float,
         "LATITUDE": float,
         "TOPO_W": float,
         "TOPO_S": float,
         "TOPO_E": float,
         "LC": str,
         "HT": float,
         "ELE": float,
         "CAN": float,
         "LAI": float,
         "k": float,
         "OH": float,
         "ELEVATION": float,
         "GRADIENT": float,
         "BOTTOM_WIDTH": float,
         "CHANNEL_ANGLE_Z": float,
         "MANNINGS_n": float,
         "SED_THERMAL_CONDUCTIVITY": float,
         "SED_THERMAL_DIFFUSIVITY": float,
         "SED_HYPORHEIC_THICKNESS": float,
         "HYPORHEIC_PERCENT": float,
         "POROSITY": float
         }

# head2var is a dictionary crosswalking input headers to model 
# variable names. This is more of a temporary fix because the model 
# variable names and input headers should probably be the same.
head2var = {"STREAM_ID": "streamID",
            "NODE_ID": "nodeID",
            "STREAM_KM": "km",
            "LONGITUDE": "longitude",
            "LATITUDE": "latitude",
            "INFLOW": "Q_in",
            "TEMPERATURE": "T_in",
            "OUTFLOW": "Q_out",
            "ELEVATION": "elevation",
            "GRADIENT": "S",
            "BOTTOM_WIDTH": "W_b",
            "CHANNEL_ANGLE_Z": "z",
            "MANNINGS_n": "n",
            "SED_THERMAL_CONDUCTIVITY": "SedThermCond",
            "SED_THERMAL_DIFFUSIVITY": "SedThermDiff",
            "SED_HYPORHEIC_THICKNESS": "SedDepth",
            "HYPORHEIC_PERCENT": "hyp_percent",
            "POROSITY": "phi",
            "Q_cont": "Q_cont",
            "d_cont": "d_cont"
            }

# iniRange is a dictionary holding min and max values for each input. 
# This is used in the input class by validate()
# to look for potential errors.
iniRange = {"STREAM_KM": [0, 999999],
            "INFLOW": [-3000, 3000],
            "TEMPERATURE": [0, 100],
            "OUTFLOW": [-3000, 3000],
            "CLOUDINESS": [0, 1],
            "WIND_SPEED": [0, 120],
            "RELATIVE_HUMIDITY": [0, 1],
            "AIR_TEMPERATURE": [-90, 60],
            "FLOW": [-3000, 3000],
            "HEIGHT": [0, 2000],
            "CANOPY": [0, 1],
            "LAI": [0, 30],
            "k": [0, 1],
            "OVERHANG": [0, 100],
            "LONGITUDE": [-180, 180],
            "LATITUDE": [-90, 90],
            "TOPO_W": [0, 90],
            "TOPO_S": [0, 90],
            "TOPO_E": [0, 90],
            "HT": [0, 2000],
            "ELE": [-450, 9000],
            "CAN": [0, 1],
            "OH": [0, 100],
            "ELEVATION": [-450, 9000],
            "GRADIENT": [0, 1],
            "BOTTOM_WIDTH": [0, 20000],
            "CHANNEL_ANGLE_Z": [0, 90],
            "MANNINGS_n": [0, 2],
            "SED_THERMAL_CONDUCTIVITY": [0, 40],
            "SED_THERMAL_DIFFUSIVITY": [0, 0.3],
            "SED_HYPORHEIC_THICKNESS": [0, 20],
            "HYPORHEIC_PERCENT": [0, 1],
            "POROSITY": [0, 1]
            }

# Placeholder for a future dictionary with more info about 
# each input variable, name, units, etc.

# keys = ["parameter variable name", "parameter full name", "parameter units", "data type", "minimum value",
# "maximum value"]

# ["DateTime", "Datetime", ]
# ["km", "Stream kilometer", "kilometer", float, 0, 999999]
# ["longitude", "Longitude", "Decimal Degrees - North American Datum 1983", float, -180, 180]
# ["latitude", "Latitude", "Decimal Degrees - North American Datum 1983", float, -90, 90]
# ["elevation", "elevation", "meters", "Morphology, Landcvoer Data", float, -450, 9000]
# ["Gradient", "Gradient", "meters/meters", float, 0.00001, 2]
# ["BottomWidth", "Stream Bottom Width", "meters", float, 0, 3000]
# ["ChannelAngleZ", "Stream Channel Angle Z", "degrees", float, 0, 90]
# ["Mannings_n", "Mannings n coefficient", "unitless", float, 0, 2]
# ["SedThermalConductivity", "Sediment Thermal Conductivity", "Watts/meter/Celsius", float]
# ["SedThermalDiffusivity", "Sediment Thermal Diffusivity", "square cm/second", float]
# ["SedHyporheicThickness", "Sediment Hyporheic Zone Thickness", "meters", float, 0, 5]
# ["%HyporheicExchange", "Percent Hyporheic exchange", "percent", float, 0, 1]
# ["Porosity", "Porosity", "unitless", "float", 0, 1]
# ["Q_cont", ]
# ["d_cont", ]
# ["T_in"]
# ["TopoWest", "Maximum Topographic Shade Angle - West", "degrees", float, 0, 90]
# ["TopoSouth", "Maximum Topographic Shade Angle - South", "degrees", float, 0, 90]
# ["TopoEast", "Maximum Topographic Shade Angle - East", "degrees", float, 0, 90]
# ["Cloudiness", "Percent Cloudiness", "Percent", float, 0, 1]
# ["WindSpeed" "Wind Speed", "meters per second", float, 0, 120]
# ["RelativeHumidity", "Relative Humidity", "unitless", float, 0, 1]
# ["AirTemp", "Air Temperature", "Celsius", float, -22, 266]
# ["Q_in", "In flow", "cubic meters per second", float,]
# ["Temp", "Water Temperature", "Celsius", float]
# ["Q_out", "Out Flow", "cubic meters per second", float]
# ["Name", "Land Cover Name", "unitless", str]
# ["Code", "Land Cover Code", "unitless", str]
# ["Height", "Land Cover Height", "meters", float]
# ["LAI", "Effective Leaf Area Index", float]
# ["k", "k Extinction Coefficient", "unitless", float]
# ["Overhang", "Landcover Overhang", float]

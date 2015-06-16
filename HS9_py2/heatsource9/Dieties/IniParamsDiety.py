""" The IniParams dictionary is simply a regular dictionary that can be 
accessed globally. We keep it here mostly because of legacy architecture from 
heatsource 8. Heatsource 8 used to optimize many of the classes with pysco 
and the primary switches for that used to be located here. Psyco is no longer 
under development and we have moved on to python 2.7 so everything pysco 
related has been removed. The C module switch used to be here also but now
we run all the routines exclusivly in python. Maybe some day we'll update 
the C code. Until then, we'll keep this simple module in place in case we need 
it for the future. Of course, one important thing to remember is to import 
this module or everything breaks.

Good thing that very important caveat is buried deeply in these
notes where no-one will ever read it."""

# Some needed updates to the C module. 
#1. Number of veg zones hard coded as 4
#2. Output from C module does not include solar blocked/passed.

from collections import defaultdict

# We need someting to initialize the dictionary so here it is.
IniParams = {"run_in_python": True,}


# This is a placeholder for a future dictionary that will contain information 
# on each of the model inputs

#keys = ["parameter variable name", "parameter full name", "parameter units", "data type", "minimum value", "maximum value"]

#["DateTime", "Datetime", ]
#["km", "Stream kilometer", "kilometer", float, 0, 999999]
#["Longitude", "Longitude", "Decimal Degrees - North American Datum 1983", float, -180, 180]
#["Latitude", "Latitude", "Decimal Degrees - North American Datum 1983", float, -90, 90]
#["Elevation", "Elevation", "meters", "Morphology, Landcvoer Data", float, -450, 9000]
#["Gradient", "Gradient", "meters/meters", float, 0.00001, 2]
#["BottomWidth", "Stream Bottom Width", "meters", float, 0, 3000]
#["ChannelAngleZ", "Stream Channel Angle Z", "degrees", float, 0, 90]
#["Mannings_n", "Mannings n coeffcient", "unitless", float, 0, 2]
#["SedThermalConductivity", "Sediment Thermal Conductivity", "Watts/meter/Celsius", float]
#["SedThermalDiffusivity", "Sediment Thermal Diffusivity", "square cm/second", float]
#["SedHyporheicThickness", "Sediment Hyporheic Zone Thicknes", "meters", float, 0, 5]
#["%HyporheicExchange", "Percent Hyporheic exchange", "percent", float, 0, 1]
#["Porosity", "Porosity", "unitless", "float", 0, 1]
#["Q_cont", ]
#["d_cont", ]
#["T_in"]
#["TopoWest", "Maximum Topographic Shade Angle - West", "degrees", float, 0, 90]
#["TopoSouth", "Maximum Topographic Shade Angle - South", "degrees", float, 0, 90]
#["TopoEast", "Maximum Topographic Shade Angle - East", "degrees", float, 0, 90]
#["Cloudiness", "Percent Cloudiness", "Percent", float, 0, 1]
#["WindSpeed" "Wind Speed", "meters per second", float, 0, 120]
#["RelativeHumidity", "Relative Humidity", "unitless", float, 0, 1]
#["AirTemp", "Air Temperature", "Celcius", float, -22, 266]
#["Q_in", "In flow", "cubic meters per second", float,]
#["Temp", "Water Temperature", "Celcius", float]
#["Q_out", "Out Flow", "cubic meters per second", float]
#["Name", "Land Cover Name", "unitless", str]
#["Code", "Land Cover Code", "unitless", str]
#["Height", "Land Cover Height", "meters", float]
#["LAI", "Effective Leaf Area Index", float]
#["k", "k Extinction Coeffcient", "unitless", float]
#["Overhang", "Landcover Overhang", float]



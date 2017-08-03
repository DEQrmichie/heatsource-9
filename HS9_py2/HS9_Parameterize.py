"""Example script to demonstrate how to use some of the Input class
methods for finer control during setup and parameterization. The methods
in this script can be used as an alternative to the more generic
HS9_Setup_... scripts.

This script must be located in the same
directory as HeatSource_control.csv. NOTE that executing this script
from Python shell (IDLE) will not identify __file__ correctly and will
result in an error. It must be executed from a command prompt. Your
options are to try to double click on this file and execute it 
using python.exe, or to open a command prompt and execute manually 
by typing:  python -i path/to/this/script/HS9_Run_Parameterize.py
"""

from heatsource9.ModelSetup.Inputs import Inputs
from heatsource9.Dieties.IniParamsDiety import IniParams
from heatsource9.Dieties.IniParamsDiety import dtype
from os.path import abspath
from os.path import dirname
from os.path import exists
from os.path import join
from os.path import realpath


def getScriptPath():
    """Gets the path to the directory where this script is
    being executed from."""
    return abspath(join(dirname(realpath(__file__)), '.'))

model_dir = getScriptPath() + '/'
control_file = 'HeatSource_Control.csv'
nodes_fc = model_dir + r"example_model.gdb\nodes"

if not exists(join(model_dir,control_file)):
    raise Exception("HeatSource_Control.csv not found. \
    Move the executable or place the control file in \
    this directory: {0}.".format(model_dir))

# create an input object
inputs = Inputs(model_dir, control_file)

# Setup control file and parameterize it
inputs.parameterize_cf(overwrite=False,
                       usertxt = "This model is an example model",
                       name = "example model", 
                       inputdir = model_dir + r"inputs/", 
                       outputdir = model_dir + r"outputs/", 
                       length = 1.8, 
                       outputkm = "all", 
                       datastart = "05/06/2003", 
                       modelstart = "07/01/2003", 
                       modelend = "07/14/2003", 
                       dataend = "09/21/2003", 
                       flushdays = 1, 
                       offset = -7, 
                       dt = 1, 
                       dx = 30, 
                       longsample = 50, 
                       bcfile = "bc.csv", 
                       inflowsites = 4, 
                       inflowinfiles = "inflow_01.csv, inflow_02.csv, inflow_03.csv, inflow_04.csv", 
                       inflowkm = "1.65, 1.5, 1.3, 0.85", 
                       accretionfile = "accretion.csv", 
                       metsites = 4, 
                       metfiles = "met_01.csv, met_02.csv, met_03.csv, met_04.csv", 
                       metkm = "1.75, 1.45, 1.10, 0.9", 
                       calcevap = "False", 
                       evapmethod = "Mass Transfer", 
                       wind_a = 1.51E-09, 
                       wind_b = 1.6E-09, 
                       calcalluvium = "True", 
                       alluviumtemp = 12.0, 
                       morphfile = "morphology.csv", 
                       lcdatafile = "lcdata.csv", 
                       lccodefile = "lccodes.csv", 
                       trans_count = 8, 
                       transsample_count = 4, 
                       transsample_distance = 8, 
                       emergent = "True", 
                       lcdatainput = "Codes", 
                       canopy_data = "Canopy", 
                       lcsamplemethod = "point", 
                       heatsource8 = "False")

# imports the control file into input object
inputs.import_control_file()

# write blank inputs
inputs.setup(use_timestamp=False, overwrite=True)

# Parameterize the lcdata and morph inputs directly 
# from nodes feature class. 
# NOTE this method currently requires use of arcpy and an active 
# ArcGIS Desktop license
inputs.parameterize_from_nodes_fc(input_file="lcdatafile", nodes_fc=nodes_fc,
                                  group_val="Example Model", grouping_field="STREAM_ID",
                                  cont_stream_km=False, overwrite=False)

inputs.parameterize_from_nodes_fc(input_file="morphfile", nodes_fc=nodes_fc,
                                  group_val="Example Model", grouping_field="STREAM_ID",
                                  cont_stream_km=False, overwrite=False)

# Parameterize the lccodes input 
lccodes = [('Active River Channel',100,0,0,0), 
           ('Barren - Clearcut',127,0,0,0), 
           ('Brush',128,1,0.4,0), 
           ('Dominate Coniferous',133,32,0.7,1.5), 
           ('Dominate Broadleaf (Riparian)',149,32,0.5,2), 
           ('Dominate Broadleaf (Upland)',150,32,0.5,2), 
           ('Road Unpaved',255,0,0,0)]

inputs.parameterize_lccodes(lccodes)
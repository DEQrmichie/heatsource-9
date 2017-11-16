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
from collections import defaultdict
from os.path import abspath
from os.path import dirname
from os.path import exists
from os.path import isfile
from os.path import join
from os.path import realpath

def parameterize_from_nodes_fc(input_file, nodes_fc,
                               group_val=None,
                               grouping_field="STREAM_ID",
                               cont_stream_km=False, overwrite=False):
    """
    Parameterize the input file using data from the TTools
    node feature class.

    input_file: The input to be parameterized. Can be either:
    "lcdatafile" or "morphfile". Other inputs must use
    parameterize() with the input data supplied as a list.

    nodes_fc: TTools derived point feature class

    group_val: Unique ID to group into a model instance.

    grouping_field: the attribute field in the nodes_fc that
    contains the group_val.

    cont_stream_km: if True overwrites the nodes km so it is
    continuous based on node_dx over all the nodes in each group

    overwrite: if True overwrites an existing file.

    """

    msg = "Parameterize from nodes fc"
    print(msg)

    # try to import arcpy
    try:
        import arcpy
    except:
        msg = "ImportError: No module named arcpy. To use this method ESRI ArcGIS must be installed"
        print(msg)
        SystemExit

    # get the headers
    if input_file == "lcdatafile":
        headers = inputs.headers_lcdata()
    elif input_file == "morphfile":
        headers = inputs.headers_morph()
    else:
        print("input_file = {0}, should be \
        'lcdatafile' or 'morphfile'".format(input_file))
        SystemExit

    # check to see if the file exists before moving on
    if not overwrite and isfile(join(IniParams["inputdir"],
                                     IniParams[input_file])):
        msg = "{0} already exists and will not be overwritten".format(IniParams[input_file])
        print(msg)
        return

    # Get all the column headers in the nodes fc
    nodes_fc_headers = [field.name for field in arcpy.Describe(nodes_fc).fields]

    # find the unique values in the grouping field
    whereclause = """{0} = '{1}'""".format(grouping_field,
                                           group_val)

    # read the node fc into the nodes dict
    nodeDict = read_nodes_fc(nodes_fc, nodes_fc_headers, whereclause)

    # build a list of the data to pass to parameterize()
    outlist = []
    nodeIDs = nodeDict.keys()
    nodeIDs.sort()

    if cont_stream_km:
        kmlist= inputs.stream_kms()

    for i, nodeID in enumerate(nodeIDs):
        row_list = []
        for header in headers:
            if header in nodeDict[nodeID].keys():
                val = nodeDict[nodeID][header]
                if cont_stream_km and header == "STREAM_KM":
                    row_list.append(kmlist[i])
                else:
                    row_list.append(val)
            else:
                # use None when there is no matching key
                row_list.append(None)
        outlist.append(row_list)

    # sort by stream km with headwaters on top
    outlist = sorted(outlist, key=itemgetter(2), reverse=True)

    inputs.write_to_csv(IniParams["inputdir"],
                        IniParams[input_file],
                        outlist, headers)

def read_nodes_fc(nodes_fc, readfields, whereclause):
    """
    Reads an input point file and returns the fields as a
    nested dictionary
    """
    # try to import arcpy
    try:
        import arcpy
    except:
        msg = "ImportError: No module named arcpy. To use this method ESRI ArcGIS must be installed"
        print(msg)
        SystemExit

    nodeDict = nested_dict()
    incursorFields = ["NODE_ID"] + readfields
    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference

    with arcpy.da.SearchCursor(nodes_fc, incursorFields, whereclause, proj) as Inrows:
        for row in Inrows:
            for f in xrange(0,len(readfields)):
                nodeDict[row[0]][readfields[f]] = row[1+f]
    return nodeDict

def nested_dict(): 
    """
    Build a nested dictionary
    """
    return defaultdict(nested_dict)

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
                       canopy_data = "CanopyClosure", 
                       lcsampmethod = "point", 
                       heatsource8 = "False")

# imports the control file into input object
inputs.import_control_file()

# write blank inputs
inputs.setup(use_timestamp=False, overwrite=True)

# Parameterize the lcdata and morph inputs directly 
# from nodes feature class. 
# NOTE this method currently requires use of arcpy and an active 
# ArcGIS Desktop license
parameterize_from_nodes_fc(input_file="lcdatafile", nodes_fc=nodes_fc,
                           group_val="Example Model", grouping_field="STREAM_ID",
                           cont_stream_km=False, overwrite=False)

parameterize_from_nodes_fc(input_file="morphfile", nodes_fc=nodes_fc,
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
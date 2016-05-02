import arcpy
import csv
import platform
from os.path import join
from os.path import exists
from os.path import isfile
from math import ceil
from time import gmtime
from time import strptime
from time import strftime
from collections import defaultdict
from calendar import timegm
from datetime import datetime
from operator import itemgetter 

from ..Dieties.IniParamsDiety import IniParams
from ..Dieties.IniParamsDiety import iniRange
from ..Dieties.IniParamsDiety import dtype
from ..Utils.Printer import Printer as print_console
import logging
logger = logging.getLogger(__name__)

class Inputs(object):
    """
    The Inputs class contains methods to read, parameterize, and write
    Heat Source input files.
    """
    def __init__(self, model_dir, control_file):
        """
        model_dir: path to the directory where the
        control file is located.

        control_file: the control file name. It must be a comma
        delimtted text file.
        
        use_timestamp: if True appends a timestamp to the end of the
        output file.
        
        overwrite: if True overwrites existing files.
        """
        
        # make these local
        self.model_dir = model_dir
        self.control_file = control_file
                    
    def headers(self, input_file="all"):
        """Returns the input file column headers"""
        if input_file == "all":
            return (self.headers_accretion(),
                    self.headers_bc(),
                    self.headers_climate(),
                    self.headers_inflow(),
                    self.headers_lcdata(),
                    self.headers_morph(),
                    self.headers_lccodes())
        
        if input_file == "accretionfile":
            return self.headers_accretion()
        
        if input_file == "bcfile":
            return self.headers_bc()
        
        if input_file == "climatefiles":
            return self.headers_climate()
        
        if input_file == "inflowinfiles":
            return self.headers_inflow()
            
        if input_file == "lcdatafile":
            return self.headers_lcdata()
        
        if input_file == "morphfile":
            return self.headers_morph()
        
        if input_file == "lccodes":
            return self.headers_lccodes()
    
    def headers_accretion(self):
        """Column headers for the accretion input file."""
        return ["STREAM_ID", "NODE_ID","STREAM_KM",
                "INFLOW","TEMPERATURE","OUTFLOW"]
    
    def headers_bc(self):
        """Column headers for the boundary condition file."""
        return ["DATETIME", "INFLOW", "TEMPERATURE"]
    
    def headers_climate(self):
        """Column headers for the climate input file(s)."""
        ncols = int((IniParams["climatesites"] /
                    len(IniParams["climatefiles"].split(","))))
        return ["DATETIME"]+["CLOUDINESS",
                             "WIND_SPEED",
                             "RELATIVE_HUMIDITY",
                             "AIR_TEMPERATURE"] * ncols
    def headers_lcdata(self):
        """
        Generates a list of the landcover data
        file column header names
        """
    
        if IniParams["lcdatainput"] == "Values":
            if IniParams["canopy_data"] == "LAI":
                #Use LAI methods
                type = ["HT","ELE","LAI","k", "OH"]
            else:        
                type = ["HT","ELE","CAN", "OH"]
        else:
            type = ["LC","ELE"]
    
        lcdataheaders =["STREAM_ID", "NODE_ID", "STREAM_KM","LONGITUDE",
                        "LATITUDE","TOPO_W","TOPO_S","TOPO_E"]
        
        if IniParams["heatsource8"]:
            dir = ["NE","E","SE","S","SW","W","NW"]
        else:        
            dir = ["T" + str(x) for x in range(1,IniParams["trans_count"]+ 1)]
    
        zone = range(1,int(IniParams["transsample_count"])+1)
    
        # Concatenate the type, dir, and zone and order in the correct way
        for t in type:
            for d in range(0,len(dir)):
                for z in range(0,len(zone)):
                    if t!="ELE" and d==0 and z==0:
                        # add emergent
                        lcdataheaders.append(t+"_T0_S0")                        
                        lcdataheaders.append(t+"_"+dir[d]+"_S"+str(zone[z]))
                    else:
                        lcdataheaders.append(t+"_"+dir[d]+"_S"+str(zone[z]))
    
        return lcdataheaders     
    
    def headers_lccodes(self):
        """
        Column headers for the land cover codes input file.
        """
        if IniParams["lcdatainput"] == "Codes":
            if IniParams["canopy_data"] == "LAI":
                return ["NAME","CODE","HEIGHT","LAI","k","OVERHANG"]
            else:
                return ["NAME","CODE","HEIGHT","CANOPY","OVERHANG"]
        else:
            return [None]
        
    def headers_cf(self):
        """
        Column headers for the control file.
        """        
        return ["LINE", "PARAMETER", "VALUE"]
    
    def headers_inflow(self):
        """
        Column headers for the tributary inflow input file(s).
        """
        if IniParams["inflowsites"] > 0:
            ncols = int(IniParams["inflowsites"] /
                        len(IniParams["inflowinfiles"].split(",")))
            return ["DATETIME"]+["FLOW","TEMPERATURE"] * ncols
        else:
            return[None]        
        
    def headers_morph(self):
        """
        Column headers for the channel morphology input file.
        """
        return ["STREAM_ID", "NODE_ID","STREAM_KM","ELEVATION",
                "GRADIENT","BOTTOM_WIDTH", "CHANNEL_ANGLE_Z",
                "MANNINGS_n","SED_THERMAL_CONDUCTIVITY",
                "SED_THERMAL_DIFFUSIVITY",
                "SED_HYPORHEIC_THICKNESSS",
                "HYPORHEIC_PERCENT","POROSITY"]
 
    def import_all(self, input_file="all"):
        """Returns all the input data"""
        
        #  read the control file and parameterize IniParams
        self.import_control_file()        
        
        if input_file == "all":
            return (self.import_accretion(),
                    self.import_bc(),
                    self.import_climate(),
                    self.import_inflow(),
                    self.import_lcdata(),
                    self.import_morph(),
                    self.import_lccodes())
        
        if input_file == "accretionfile":
            return self.import_accretion()
        
        if input_file == "bcfile":
            return self.import_bc()
        
        if input_file == "climatefiles":
            return self.import_climate()
        
        if input_file == "inflowinfiles":
            return self.import_inflow()
            
        if input_file == "lcdatafile":
            return self.import_lcdata()
        
        if input_file == "morphfile":
            return self.import_morph()
        
        if input_file == "lccodes":
            return self.import_lccodes()
        
        
    def import_accretion(self):
        """Returns the accretion input data"""
        return self.read_to_dict(IniParams["inputdir"],
                                         IniParams["accretionfile"],
                                         self.headers_accretion())
        
    def import_bc(self, return_list=True, skiprows=1, skipcols=1):
        """Returns the boundary condition data"""
        data = self.read_to_dict(IniParams["inputdir"],
                                IniParams["bcfile"],
                                self.headers_bc())
        
        if return_list:
            return self.dict2list(data, self.headers_bc(),
                           skiprows=1, skipcols=1)
        else:
            return data
    
    def import_climate(self, skiprows=1, skipcols=1):
        """Returns the  climate input data"""
        # split by comma if multiple files
        filenames = IniParams["climatefiles"].split(",")
    
        filenames = [filenames] if not isinstance(filenames, list) else filenames
        
        for i, filename in enumerate(filenames):
            d1 = self.read_to_dict(IniParams["inputdir"],
                                     filename,
                                     self.headers_climate())
            
            d2 = self.dict2list(d1, self.headers_climate(),
                                   skiprows, skipcols)
                
            if i == 0:
                data = d2
            else:
                data = [data[i] + row for i, row in enumerate(d2)]
        
        return data         

    def import_control_file(self):
        """Returns the control file"""
        if not exists(join(self.model_dir, self.control_file)):
            logging.ERROR("HeatSource_Control.csv not \
            found {0}".format(join(self.model_dir,self.control_file)))
            
            raise Exception("HeatSource_Control.csv not found. \
            Move the executable or place the control file in \
            this directory: {0}.".format(self.model_dir))
        
        print_console("Reading control file")
        
        cf_dict = self.control_file_dict()
        cf_list = self.read_to_list(self.model_dir,
                                 self.control_file,
                                 skiprows=1, skipcols=0)
        
        # set up a list to check if a missing value in
        # the contorl file is ok
        if IniParams["run_type"] == 0:
            # For temperature runs None is ok for these inputs
            none_ok = ["usertxt", "name"]
            
        elif IniParams["run_type"] == 1:
            # For solar runs None is ok for these inputs
            none_ok = ["usertxt", "name","flushdays","bcfile",
                       "inflowsites", "inflowinfiles", "inflowkm",
                       "accretionfile", "calcevap", "evapmethod",
                       "wind_a", "wind_b", "calcalluvium", "alluviumtemp"]
                    
        elif IniParams["run_type"] == 2:
            # For hydraulic runs None is ok for these inputs
            none_ok = ["usertxt", "name","lcdatafile", "lccodefile",
                       "trans_count", "transsample_count",
                       "transsample_distance", "emergent",
                       "lcdatainput", "canopy_data", "vegDistMethod",
                       "point"]
        else:
            # This is a setup call, None is ok for all of them
            none_ok = IniParams.keys()
        
            # This is so the iteration happens in descending order
            # so some of the keys are parameterized earlier for the
            # none list. 
            keys = cf_dict.keys()
            keys.sort(reverse=True)
            
            for k in keys:
                v = cf_dict[k]        
                for line in cf_list:
                    if (line[1] == v[1]):
                        # check for missing values
                        if str.strip(line[2]) in ["", None]:
                            if k in none_ok:
                                IniParams[k] = None
                            elif (k == "lccodefile" and
                                IniParams["lcdatainput"] == "Values"):
                                # None ok because lccodes file is not needed
                                IniParams[k] = None
                                
                            elif (k in ["inflowinfiles","inflowkm"] and
                                  IniParams["inflowsites"] == 0):
                                # None ok because there are no inflow sites
                                IniParams[k] = None
                                
                            elif (k == "alluviumtemp" and
                                  IniParams["calcalluvium"] is False):
                                # None ok because not including 
                                # deep alluvium temps
                                IniParams[k] = None
                            else:
                                raise TypeError("Value in control file line {0} is missing".format(line[0]))
                    # now make sure it's the correct data type
                    elif dtype[k] is basestring:
                        IniParams[k] = str.strip(line[2])
                    
                    elif dtype[k] is int:
                        IniParams[k] = int(float(str.strip(line[2])))
                        
                    elif dtype[k] is bool:
                        if str.strip(line[2]).upper() in ["TRUE", "FALSE"]:
                            IniParams[k] = str.strip(line[2]).upper() == "TRUE"
                        else:
                            raise TypeError("Control file line {0} must be TRUE or FALSE".format(line[0]))
                            
                    elif dtype[k] is float:
                        IniParams[k] = float(str.strip(line[2]))
                        
        # Make dates into seconds since UTC epoch
        IniParams["datastart"] = timegm(strptime(IniParams["datastart"] + " 00:00:00" ,"%m/%d/%Y %H:%M:%S"))
        IniParams["dataend"] = timegm(strptime(IniParams["dataend"],"%m/%d/%Y")) + 86400
    
        if IniParams["modelstart"] is None:
            IniParams["modelstart"] = IniParams["datastart"]
        else:
            IniParams["modelstart"] = timegm(strptime(IniParams["modelstart"] + " 00:00:00","%m/%d/%Y %H:%M:%S"))
    
        if IniParams["modelend"] is None:
            IniParams["modelend"] = IniParams["dataend"]
        else:
            IniParams["modelend"] = timegm(strptime(IniParams["modelend"],"%m/%d/%Y")) + 86400
    
        IniParams["flushtimestart"] = IniParams["modelstart"] - IniParams["flushdays"]*86400
                        
        # If the number of transverse samples per direction 
        # is NOT reported, assume 4 (old default)
        if not IniParams["transsample_count"]:
            IniParams["transsample_count"] = 4.0 
    
        # Format for heat source 8 methods same 
        # as 8 directions but no north
        if IniParams["heatsource8"]:
            IniParams["trans_count"] = 7
            IniParams["sample_count"] = int(IniParams["transsample_count"] * 7)
        else:
            # Set the total number landcover sample count (0 = emergent)
            IniParams["sample_count"] = int(IniParams["transsample_count"]
                                            * IniParams["trans_count"])
            
        # Set up our evaporation method
        IniParams["penman"] = False
        if IniParams["calcevap"]:
            IniParams["penman"] = True if IniParams["evapmethod"] == "Penman" else False
    
        # make sure that the timestep divides into 60 minutes, 
        # or we may not land squarely on each hour's starting point.
        if float(60)/IniParams["dt"] - int(float(60)/IniParams["dt"]) > 1e-7:
            raise ValueError("I'm sorry, your timestep ({0}) must evenly divide into 60 minutes.".format(IniParams["dt"]))
        else:
            # make dt measured in seconds
            IniParams["dt"] = IniParams["dt"]*60
    
        # Make sure the output directory ends in a 
        # slash based on system platform
        if (platform.system() == "Windows" and IniParams["outputdir"][-1] != "\\"):
            raise ValueError("Output directory needs to have a backslash at end of the path. ..\\outputfolder\\")
    
        if (platform.system() == "Darwin" and IniParams["outputdir"][-1] != "/"):
            raise ValueError("Output directory needs to have a forward slash at the end of the path. ../outputfolder/")    
    
        # the distance step must be an exact, greater or equal to one, 
        # multiple of the sample rate.
        if (IniParams["dx"]%IniParams["longsample"]
            or IniParams["dx"]<IniParams["longsample"]):
            raise ValueError("Distance step (dx) must be a multiple of the longitudinal stream sample distance")
      

    def control_file_dict(self, **kwargs):
        """
        Returns a dictionary of the control file lines
        with empty values. Dictionary keys coorespond to the
        keys in IniParams.
        """
        
        #TODO RM fix so dict is related to numbers or text symbols
        cf_dict = {"usertxt": [2, "USER NOTES", None],
               "name": [3, "SIMULATION NAME", None],
               "inputdir": [4, "INPUT PATH", None],
               "outputdir": [5, "OUTPUT PATH", None],
               "length": [6, "STREAM LENGTH (KILOMETERS)", None],
               "outputkm": [7, "OUTPUT KILOMETERS", None],
               "datastart": [8, "DATA START DATE (mm/dd/yyyy)", None],
               "modelstart": [9, "MODELING START DATE (mm/dd/yyyy)", None],
               "modelend": [10, "MODELING END DATE (mm/dd/yyyy)", None],
               "dataend": [11, "DATA END DATE (mm/dd/yyyy)", None],
               "flushdays": [12, "FLUSH INITIAL CONDITION (DAYS)", None],
               "offset": [13, "TIME OFFSET FROM UTC (HOURS)", None],
               "dt": [14, "MODEL TIME STEP - DT (MIN)", None],
               "dx": [15, "MODEL DISTANCE STEP - DX (METERS)", None],
               "longsample": [16, "LONGITUDINAL STREAM SAMPLE DISTANCE (METERS)", None],
               "bcfile": [17, "BOUNDARY CONDITION FILE NAME", None],
               "inflowsites": [18, "TRIBUTARY SITES", None],
               "inflowinfiles": [19, "TRIBUTARY INPUT FILE NAMES", None],
               "inflowkm": [20, "TRIBUTARY MODEL KM", None],
               "accretionfile": [21, "ACCRETION INPUT FILE NAME", None],
               "climatesites": [22, "CLIMATE DATA SITES", None],
               "climatefiles": [23, "CLIMATE INPUT FILE NAMES", None],
               "climatekm": [24, "CLIMATE MODEL KM", None],
               "calcevap": [25, "INCLUDE EVAPORATION LOSSES FROM FLOW (TRUE/FALSE)", None],
               "evapmethod": [26, "EVAPORATION METHOD (Mass Transfer/Penman)", None],
               "wind_a": [27, "WIND FUNCTION COEFFICIENT A", None],
               "wind_b": [28, "WIND FUNCTION COEFFICIENT B", None],
               "calcalluvium": [29, "INCLUDE DEEP ALLUVIUM TEMPERATURE (TRUE/FALSE)", None],
               "alluviumtemp": [30, "DEEP ALLUVIUM TEMPERATURE (*C)", None],
               "morphfile": [31, "MORPHOLOGY DATA FILE NAME", None],
               "lcdatafile": [32, "LANDCOVER DATA FILE NAME", None],
               "lccodefile": [33, "LANDCOVER CODES FILE NAME", None],
               "trans_count": [34, "NUMBER OF TRANSECTS PER NODE", None],
               "transsample_count": [35, "NUMBER OF SAMPLES PER TRANSECT", None],
               "transsample_distance": [36, "DISTANCE BETWEEN TRANSESCT SAMPLES (METERS)", None],
               "emergent": [37, "ACCOUNT FOR EMERGENT VEG SHADING (TRUE/FALSE)", None],
               "lcdatainput": [38, "LANDCOVER DATA INPUT TYPE (Codes/Values)", None],
               "canopy_data": [39, "CANOPY DATA TYPE (LAI/CanopyCover)", None],
               "vegDistMethod": [40, "VEGETATION ANGLE CALCULATION METHOD (point/zone)", None],
               "heatsource8": [41, "USE HEAT SOURCE 8 LANDCOVER METHODS (TRUE/FALSE)", None],
               }
        
        return cf_dict    
        
    def import_lccodes(self):
        """Returns the land cover codes data."""
        return self.read_to_dict(IniParams["inputdir"],
                                 IniParams["lccodefile"],
                                 self.headers_lccodes())
    
    def import_lcdata(self, return_list=True, skiprows=1, skipcols=2):
        """Returns the land cover data."""    
        data = self.read_to_dict(IniParams["inputdir"],
                                 IniParams["lcdatafile"],
                                 self.headers_lcdata())
        
        if return_list:
            return self.dict2list(data, self.headers_lcdata(),
                                  skiprows, skipcols)
        else:
            return data        
                                
    def import_elevations(self):
        """Returns the land cover elevation data."""
        return [None]
        
    def import_canopy(self):
        """Returns the land cover canopy data."""
        return [None]
    
    def import_k(self):
        """Returns k, the riparian extintion coeffcient data."""
        return [None]     
    
    def import_inflow(self, return_list=True, skiprows=1, skipcols=1):
        """Returns the tributary inflow data."""
        if IniParams["inflowsites"] > 0:
            # split by comma if multiple files
            filenames = IniParams["inflowinfiles"].split(",")
        
            filenames = [filenames] if not isinstance(filenames, list) else filenames
            
            for i, filename in enumerate(filenames):
                d1 = self.read_to_dict(IniParams["inputdir"],
                                         filename,
                                         self.headers_inflow())
                
                d2 = self.dict2list(d1, self.headers_inflow(),
                                       skiprows, skipcols)
                    
                if i == 0:
                    data = d2
                else:
                    data = [data[i] + row for i, row in enumerate(d2)]
            
            return data
            
        else:
            return []        
        
    def import_morph(self, return_list=False):
        """Returns the channel morphology data"""
        return self.read_to_dict(IniParams["inputdir"],
                                 IniParams["morphfile"],
                                 self.headers_morph())    
            
    def parameterize_from_nodes_fc(self, input_file, nodes_fc,
                                   group_val=None,
                                   grouping_field="STREAM_ID"):
        """
        Paramaterize the input file using data from the TTools
        node feature class.
        
        input_file: The input to be paramatierized. Can be either:
        "lcdatafile" or "morphfile". Other inputs must use
        paramaterize() with the input data supplied as a list.
        
        nodes_fc: TTools derived point feature class
            
        group_val: Unique ID to group into a model instance.
            
        grouping_field: the attribute field in the nodes_fc that
        contains the group_val.
        
        """
        # get the headers
        if input_file == "lcdatafile":
            headers = self.headers_lcdata()
        elif input_file == "morphfile":
            headers = self.headers_morph()
        else:
            logging.ERROR("input_file = {0}, should be \
            'lcdatafile' or 'morphfile'".format(inpu_file))
                
        # Get all the column headers in the nodes fc
        nodes_fc_headers = [field.name for field in arcpy.Describe(nodes_fc).fields]
    
        # find the unique values in the grouping field
        whereclause = """{0} = '{1}'""".format(grouping_field,
                                               group_val)
    
        # read the node fc into the nodes dict
        self.nodeDict = self.read_nodes_fc(nodes_fc, nodes_fc_headers,
                                           whereclause)        
        
        # build a list of the data to pass to paramaterize()
        outlist = []
        nodeIDs = self.nodeDict.keys()
        nodeIDs.sort()                 
        for nodeID in nodeIDs:
            row_list = []
            for header in headers:
                if header in self.nodeDict[nodeID].keys():
                    val = self.nodeDict[nodeID][header]
                    row_list.append(val)
                else:
                    # use None when there is no matching key
                    row_list.append(None)
            outlist.append(row_list)
                
        # sort by stream km with headwaters on top
        outlist = sorted(outlist, key=itemgetter(2), reverse=True)
        
        self.write_to_csv(IniParams["inputdir"],
                          IniParams[input_file],
                          headers, data)
        
    def nested_dict(self): 
        """
        Build a nested dictionary
        """
        return defaultdict(self.nested_dict)
    
    def dict2list(self, data, colnames, skiprows=0, skipcols=0):
        
        d2 = zip(*[[k] + data[k] for k in colnames])
        
        # skiprows
        d3 = d2[-(len(d2)-skiprows):]
        
        # skip cols
        d4 = [line[+skipcols:] for line in d3]      
        
        return [list(i) for i in d4]
    
    def parameterize_cf(self, overwrite=False, **kwargs):
        """
        Writes the control file. Any keyword arguments
        passed will be parameterized into the control file.
        """
        cf_dict = self.control_file_dict()
        
        for k, v in kwargs.items():
            cf_dict[k][2] = v
            
        # sort the list is in the order of the line number
        cf_sorted = sorted(cf_dict.items(), key=itemgetter(1))
        cf_list = [line[1] for line in cf_sorted]        

        if overwrite:
            self.write_to_csv(self.model_dir, self.control_file,
                              cf_list, self.headers_cf())
        else:
            # check to see if the file exists and write if not
            if not isfile(join(self.model_dir, self.control_file)):
                self.write_to_csv(self.model_dir, self.control_file,
                                  cf_list, self.headers_cf())
            else:
                raise Exception("HeatSource_Control.csv already exists.")
        
    def parameterize_lccodes(self, lcdata=None, code_as_ht=False, can=None,
                             lai=None, k=None, oh=None):
        """
        This function identifes unique codes from the land cover
        data file and generates the land cover code file
        based on all the unique codes. If a list of the land cover data
        is not passed as an agrument the file identifed in the cotrol file
        will be read instead. Optional default values can
        be based for the other parameters.
        """
        
        if lcdata is None:
            # read the land cover data
            lcdata = self.import_lcdata(return_list=True,
                                    skiprows=1, skipcols=0)      
        
        codes = set()
        
        for row in lcdata:
            for col in range(8, (IniParams["trans_count"] *
                                 IniParams["transsample_count"]) + 9):
                codes.add(str(row[col]))
        
        codes = list(codes)
        codes.sort()
        
        if IniParams["canopy_data"] == "LAI":
            
            if code_as_ht:
                lccodes = [[None, code, float(code), lai, k, oh] for code in codes]
            else:
                lccodes = [[None, code, None, lai, k, oh] for code in codes]
            
        else:
            if code_as_ht:
                lccodes = [[None, code, float(code), can, oh] for code in codes]
            else:
                lccodes = [[None, code, None, can, oh] for code in codes]
            
        self.write_to_csv(IniParams["inputdir"],
                          IniParams["lccodefile"],
                          lccodes, self.headers_lccodes())    
        
    def read_nodes_fc(self, nodes_fc, readfields, whereclause):
        """
        Reads an input point file and returns the fields as a
        nested dictionary
        """
        nodeDict = self.nested_dict()
        incursorFields = ["NODE_ID"] + readfields
        # Determine input point spatial units
        proj = arcpy.Describe(nodes_fc).spatialReference
                
        with arcpy.da.SearchCursor(nodes_fc, incursorFields, whereclause, proj) as Inrows:
            for row in Inrows:
                for f in xrange(0,len(readfields)):
                    nodeDict[row[0]][readfields[f]] = row[1+f]
        return nodeDict              
    
    def read_to_dict(self, inputdir, filename, colnames):
        """
        Reads a comma delimited text file and returns the data
        as a dictionary with the column header as the key.
        """
        
        # each value in each column is appended to a list
        data=defaultdict(list)
        
        with open(join(inputdir, filename.strip()), "rU") as file_object:
            reader = csv.DictReader(file_object, dialect="excel")
            
            # set the colnames as the dictionary key 
            reader.fieldnames = colnames
            # skip the header row
            reader.next()
            # read a row as {column1: value1, column2: value2,...}
            for row in reader:
                # go over each column name and value
                for k, v in row.items():
                    
                    # if the value is empty '' replace it with a None
                    if v.strip() in ['', None]:
                        v = None
                    # append the value into the appropriate 
                    # list based on column name k
                    #if there is more than one cl
                    data[k].append(v)
                    
        return self.validate(data)    
    
        
    def read_to_list(self, inputdir, filenames, skiprows, skipcols):
        """Reads a comma delimtted text file into a
        list of lists indexed by row number. If there is more than
        one file it appends the data to the row.
        
        adds another index for the number of filenames.
        The return data takes this form:
        data[row][filenameindex][filecolumn]"""
        
        # split by comma if multiple files
        filenames = filenames.split(",")
    
        filenames = [filenames] if not isinstance(filenames, list) else filenames
        i = 1
        for filename in filenames:
            with open(join(inputdir, filename.strip()), "rU") as file_object:
                newfile=[row for row in csv.reader(file_object.read().splitlines(), dialect="excel")]
                numcols = len(newfile[0])
            
            # skip rows    
            newfile = newfile[-(len(newfile)-skiprows):]
            # skip cols
            newfile = [line[+skipcols:] for line in newfile]
            if i == 1:
                data = newfile
            else:
                data = [data[row] + newfile[row] for row in range(0,len(newfile))]
            i = i + 1
        return data
    
    def setup(self, use_timestamp=True, overwrite=True):
        """
        Writes blank input files based on settings in the control file
        """

        print_console("Starting input file setup")
        
        # check if the input/output dirctories exist and 
        # create them if not
        if not exists(IniParams["inputdir"]):
            makedirs(IniParams["inputdir"])
            
        if not exists(IniParams["outputdir"]):
            makedirs(IniParams["outputdir"])         

        now = datetime.now()    
        timestamp = now.strftime("%Y%m%d_%H%M%S")
        timelist= self.datetime_string()	
        kmlist= self.stream_kms()

        # For inflow and climate data there can be a single input 
        # file for each node OR just one input file with multiple 
        # nodes in the file, This creates a list of the trib 
        # and climate file names if there is more than one   
        if IniParams["inflowsites"] > 0:
            tribfiles = IniParams["inflowinfiles"].split(",")
        climatefiles = IniParams["climatefiles"].split(",")	
        
        acclist= [[None, None, km]+[None]*3 for km in kmlist]
        bclist= [[t, None, None] for t in timelist]
        lcdatalist= [[None, None, km]+[None]*
                     (len(self.headers_lcdata())-3) for km in kmlist]
        morphlist= [[None, None, km]+[None]*10 for km in kmlist]
        
        # This writes to csv using the file name from the control 
        # file and adds a timestamp
        
        if use_timestamp:
            acc_file = "input_"+timestamp+"_"+IniParams["accretionfile"]
            bc_file = "input_"+timestamp+"_"+IniParams["bcfile"]
            lcdata_file = "input_"+timestamp+"_"+IniParams["lcdatafile"]
            if IniParams["lcdatainput"] == "Codes":
                lccodes_file = "input_"+timestamp+"_"+IniParams["lccodefile"]
            morph_file = "input_"+timestamp+"_"+IniParams["morphfile"]
        else:
            # Name the files as they appear in the control file
            acc_file = IniParams["accretionfile"]
            bc_file = IniParams["bcfile"]
            lcdata_file = IniParams["lcdatafile"]
            if IniParams["lcdatainput"] == "Codes":
                lccodes_file = IniParams["lccodefile"]
            morph_file = IniParams["morphfile"]
        
        print_console("Writing empty csv files")
        
        if overwrite:
            # overwrite the inputs regardless if they exist or not
            self.write_to_csv(IniParams["inputdir"], acc_file,
                              acclist, self.headers_accretion())
            
            self.write_to_csv(IniParams["inputdir"], bc_file,
                              bclist, self.headers_bc())
            
            self.write_to_csv(IniParams["inputdir"], lcdata_file,
                              lcdatalist, self.headers_lcdata())

            if IniParams["lcdatainput"] == "Codes":
                self.write_to_csv(IniParams["inputdir"], lccodes_file,
                                      [[None]], self.headers_lccodes())
            else:
                print_console("...Landcover input type = Values. land cover codes file not written")
                
            self.write_to_csv(IniParams["inputdir"], morph_file,
                             morphlist, self.headers_morph())
            
        else:
            # check to see if the files exists and write if not
            if not isfile(join(IniParams["inputdir"], acc_file)):
                self.write_to_csv(IniParams["inputdir"], acc_file,
                                  acclist, self.headers_accretion())
                
            if not isfile(join(IniParams["inputdir"], bc_file)):
                self.write_to_csv(IniParams["inputdir"], bc_file,
                                  bclist, self.headers_bc())
                
            if not isfile(join(IniParams["inputdir"], lcdata_file)):
                self.write_to_csv(IniParams["inputdir"], lcdata_file,
                                  lcdatalist, self.headers_lcdata())
                
            if IniParams["lcdatainput"] == "Codes":
                if not isfile(join(IniParams["inputdir"], lccodes_file)):
                    self.write_to_csv(IniParams["inputdir"], lccodes_file,
                                          [[None]], self.headers_lccodes())
            else:
                print_console("...Landcover input type = Values. land cover codes file not written")
                
            if not isfile(join(IniParams["inputdir"], morph_file)):
                self.write_to_csv(IniParams["inputdir"], morph_file,
                                  morphlist, self.headers_morph())

        for file in climatefiles:
            climatelist = [[t] + [None]*4*int((IniParams["climatesites"]/len(climatefiles))) for t in timelist]
            
            if use_timestamp:
                climate_filename = "input_"+timestamp+"_"+file.strip()
            else:
                climate_filename = file.strip()
            
            if overwrite:
                self.write_to_csv(IniParams["inputdir"],
                                  climate_filename,
                                  climatelist, self.headers_climate())
            else:
                if not isfile(join(IniParams["inputdir"], climate_filename)):
                    self.write_to_csv(IniParams["inputdir"],
                                      climate_filename,
                                      climatelist, self.headers_climate())

        if IniParams["inflowsites"] > 0:
            for file in tribfiles:
                inflowlist = [[t] + [None]*2*int((IniParams["inflowsites"]/len(tribfiles))) for t in timelist]
                
                if use_timestamp:
                    trib_filename = "input_"+timestamp+"_"+file.strip()
                else:
                    trib_filename = file.strip()
                    
                if overwrite:
                    self.write_to_csv(IniParams["inputdir"],trib_filename,
                                      inflowlist, self.headers_inflow())
                else:
                    if not isfile(join(IniParams["inputdir"], trib_filename)):
                        self.write_to_csv(IniParams["inputdir"],trib_filename,
                                          inflowlist, self.headers_inflow())
        else:
            print_console("Inflow Sites = 0, inflow file not written")

        print_console("Finished input file setup")
        
    def datetime_string(self):
        """
        Return a date time list in string format (MM/DD/YYYY HH:MM)
        corresponding to the data start and end dates available
        in the control file
        """
        timelist = []
        # hourly timestep
        timelist = range(IniParams["datastart"],IniParams["dataend"]+60,3600) 
        for i, val in enumerate(timelist):
            timelist[i] = strftime("%m/%d/%Y %H:%M",gmtime(val))
        return timelist
    
    def stream_kms(self):
        """Return a list of stream kilometers sorted from
        headwaters to mouth"""
        
        # Since the stream length is formatted as a floating point number 
        # we can sometimes run into trouble when calculating the number 
        # of nodes with dx due to ceil function rounding up on 
        # floating points that are not exact representaitons of the 
        # input value. Therefore we enforce a precision only up to the 
        # ten thousandths decimal place to make sure the number of nodes
        # is correct.        
        num_nodes = int(ceil(round(IniParams["length"]*1000/(IniParams["longsample"]), 4))) +1
        kmlist = []    
        kmlist = [(node * IniParams["longsample"])/1000 for node in range(0,num_nodes)]    
        kmlist.sort(reverse=True)
        return kmlist
    
    def validate(self, data):
        """
        Checks the data type and range.
        If value is blank (None), returns 0.0 for float, 0 for integer,
        and None for string data types.
        """
        data_v = {}
        for key, v in data.iteritems():
                        
            if key not in dtype.keys():
                # This is to find the correct landcover data 
                # key since they are all different
                if "LC" in key:
                    # basestring
                    k = "LC"
                elif any(s in key for s in ["HT","ELE", "LAI", "k", "CAN", "OH"]):
                    # float
                    k = "ELE"
            else:
                k = key
            
            # -- 
                    
            if dtype[k] is float:
                data_v[key] = [float(i) if i is not None else 0.0 for i in v]
            
            elif dtype[k] is int:
                data_v[key] = [int(float(i)) if i is not None else 0 for i in v]
                
            elif dtype[k] is basestring:
                data_v[key] = [str(i) if i is not None else None for i in v]
                
            # --
                
            if (dtype[k] is not basestring and dtype[k] in iniRange.keys()):
                # check the value range
                for val, i in enumerate(v):
                    if not iniRange[k][0] <= val <= iniRange[k][1]:
                        raise ValueError("The value ({0} in {1} is out of range".format(val, k))
        
        return data_v
        
    def write_to_csv(self, outputdir, filenames, outlist, colnames=None):
        """Write the outlist to a comma delimitted text file
        colnames are optional.
        """
        
        # split by comma if multiple files
        filenames = filenames.split(",")        
        
        filenames = [filenames] if not isinstance(filenames, list) else filenames
        
        for filename in filenames:
            with open(join(outputdir, filename), "wb") as file_object:
                writer = csv.writer(file_object,  dialect= "excel")
                if colnames:
                    writer.writerow(colnames)
                writer.writerows(outlist)
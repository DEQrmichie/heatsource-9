from ..Dieties.IniParamsDiety import IniParams
from ..Dieties.IniParamsDiety import iniRange
from ..Dieties.IniParamsDiety import dtype
from ..Utils.Printer import Printer as print_console

import csv
import platform
from os import makedirs
from os.path import exists
from os.path import isfile
from os.path import join
from math import ceil
from time import gmtime
from time import strptime
from time import strftime
from string import digits
from collections import defaultdict
from calendar import timegm
from datetime import datetime
from operator import itemgetter

import logging

# set up logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    filename='heatsource.log',
                    filemode='w')
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

        """

        # make these local
        self.model_dir = model_dir
        self.control_file = control_file

        # For when reading and writing from excel is availiable
        self.read_to_list = self.read_csv_to_list
        self.read_to_dict = self.read_csv_to_dict
        self.write_to_output = self.write_to_csv
        self.setup = self.setup_csv

    def headers(self, input_file="all"):
        """Returns the input file column headers"""
        if input_file == "all":
            return (self.headers_accretion(),
                    self.headers_bc(),
                    self.headers_met(),
                    self.headers_inflow(),
                    self.headers_lcdata(),
                    self.headers_morph(),
                    self.headers_lccodes())

        if input_file == "accretionfile":
            return self.headers_accretion()

        if input_file == "bcfile":
            return self.headers_bc()

        if input_file == "metfiles":
            return self.headers_met()

        if input_file == "inflowinfiles":
            return self.headers_inflow()

        if input_file == "lcdatafile":
            return self.headers_lcdata()

        if input_file == "morphfile":
            return self.headers_morph()

        if input_file == "lccodes":
            return self.headers_lccodes()

    def headers_accretion(self):
        """Returns a list of column headers for
        the accretion input file."""
        return ["STREAM_ID", "NODE_ID", "STREAM_KM",
                "INFLOW", "TEMPERATURE", "OUTFLOW"]

    def headers_bc(self):
        """Returns a list of column headers for
        the boundary condition file."""
        return ["DATETIME", "INFLOW", "TEMPERATURE"]

    def headers_met(self):
        """Returns a list of column headers for
        the met input file(s)."""
        ncols = int((IniParams["metsites"] /
                     len(IniParams["metfiles"].split(","))))
        header = ["DATETIME"]
        for n in range(1, ncols + 1):
            header.append("CLOUDINESS" + str(n))
            header.append("WIND_SPEED" + str(n))
            header.append("RELATIVE_HUMIDITY" + str(n))
            header.append("AIR_TEMPERATURE" + str(n))
        return header

    def headers_lcdata(self):
        """
        Returns a list of column headers for
        the land cover data file.
        """

        if IniParams["lcdatainput"] == "Values":
            if IniParams["canopy_data"] == "LAI":
                # Use LAI methods
                prefix = ["HT", "ELE", "LAI", "k", "OH"]
            else:
                prefix = ["HT", "ELE", "CAN", "OH"]
        else:
            prefix = ["LC", "ELE"]

        lcdataheaders = ["STREAM_ID", "NODE_ID", "STREAM_KM", "LONGITUDE",
                         "LATITUDE", "TOPO_W", "TOPO_S", "TOPO_E"]

        if IniParams["heatsource8"]:
            tran = ["NE", "E", "SE", "S", "SW", "W", "NW"]
        else:
            tran = ["T" + str(x) for x in range(1, IniParams["trans_count"] + 1)]

        zone = range(1, int(IniParams["transsample_count"]) + 1)

        # Concatenate the prefix, transect, and zone and order in the correct way
        for p in prefix:
            for t in range(0, len(tran)):
                for z in range(0, len(zone)):
                    if p != "ELE" and t == 0 and z == 0:
                        # add emergent
                        lcdataheaders.append(p + "_T0_S0")
                        lcdataheaders.append(p + "_" + tran[t] + "_S" + str(zone[z]))
                    else:
                        lcdataheaders.append(p + "_" + tran[t] + "_S" + str(zone[z]))

        return lcdataheaders

    def headers_lccodes(self):
        """
        Returns a list of column headers for
        the land cover codes input file.
        """
        if IniParams["lcdatainput"] == "Codes":
            if IniParams["canopy_data"] == "LAI":
                return ["NAME", "CODE", "HEIGHT", "LAI", "k", "OVERHANG"]
            else:
                return ["NAME", "CODE", "HEIGHT", "CANOPY", "OVERHANG"]
        else:
            return [None]

    def headers_cf(self):
        """
        Returns a list of column headers for
        the control file.
        """
        return ["LINE", "PARAMETER", "KEY", "VALUE"]

    def headers_inflow(self):
        """
        Returns a list of column headers for
        the tributary inflow input file(s).
        """
        if IniParams["inflowsites"] > 0:
            ncols = int(IniParams["inflowsites"] /
                        len(IniParams["inflowinfiles"].split(",")))
            header = ["DATETIME"]
            for n in range(1, ncols + 1):
                header.append("FLOW" + str(n))
                header.append("TEMPERATURE" + str(n))
            return header
        else:
            return [None]

    def headers_morph(self):
        """
        Returns a list of column headers for
        the channel morphology input file.
        """
        return ["STREAM_ID", "NODE_ID", "STREAM_KM", "ELEVATION",
                "GRADIENT", "BOTTOM_WIDTH", "CHANNEL_ANGLE_Z",
                "MANNINGS_n", "SED_THERMAL_CONDUCTIVITY",
                "SED_THERMAL_DIFFUSIVITY",
                "SED_HYPORHEIC_THICKNESSS",
                "HYPORHEIC_PERCENT", "POROSITY"]

    def import_all(self, input_file="all"):
        """Returns all the input data"""

        #  read the control file and parameterize IniParams
        self.import_control_file()

        if input_file == "all":
            return (self.import_accretion(),
                    self.import_bc(),
                    self.import_met(),
                    self.import_inflow(),
                    self.import_lcdata(),
                    self.import_morph(),
                    self.import_lccodes())

        if input_file == "accretionfile":
            return self.import_accretion()

        if input_file == "bcfile":
            return self.import_bc()

        if input_file == "metfiles":
            return self.import_met()

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

    def import_met(self, skiprows=1, skipcols=1):
        """Returns the  met input data"""
        # split by comma if multiple files
        filenames = IniParams["metfiles"].split(",")

        filenames = [filenames] if not isinstance(filenames, list) else filenames

        for i, filename in enumerate(filenames):
            d1 = self.read_to_dict(IniParams["inputdir"],
                                   filename,
                                   self.headers_met())

            d2 = self.dict2list(d1, self.headers_met(),
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
            found {0}".format(join(self.model_dir, self.control_file)))

            raise Exception("HeatSource_Control.csv not found. \
            Move the executable or place the control file in \
            this directory: {0}.".format(self.model_dir))

        msg = "Reading control file"
        logging.info(msg)
        print_console(msg)

        cf_dict = self.control_file_dict()
        cf_list = self.read_to_list(self.model_dir,
                                    self.control_file,
                                    skiprows=1, skipcols=0)

        # set up a list to check if a missing value in
        # the control file is ok
        if IniParams["run_type"] == 0:
            # For temperature runs None is ok for these inputs
            none_ok = ["usertxt", "name"]

        elif IniParams["run_type"] == 1:
            # For solar runs None is ok for these inputs
            none_ok = ["usertxt", "name", "flushdays", "bcfile",
                       "inflowsites", "inflowinfiles", "inflowkm",
                       "accretionfile", "calcevap", "evapmethod",
                       "wind_a", "wind_b", "calcalluvium", "alluviumtemp"]

        elif IniParams["run_type"] == 2:
            # For hydraulic runs None is ok for these inputs
            none_ok = ["usertxt", "name", "lcdatafile", "lccodefile",
                       "trans_count", "transsample_count",
                       "transsample_distance", "emergent",
                       "lcdatainput", "canopy_data", "lcsampmethod",
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
            for line in cf_list:
                if (line[2] == k):
                    # check for missing values
                    if str.strip(line[3]) in ["", None]:
                        if k in none_ok:
                            IniParams[k] = None
                        elif (k == "lccodefile" and
                              IniParams["lcdatainput"] == "Values"):
                            # None ok because lccodes file is not needed
                            IniParams[k] = None

                        elif (k in ["inflowinfiles", "inflowkm"] and
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
                        IniParams[k] = str.strip(line[3])

                    elif dtype[k] is int:
                        IniParams[k] = int(float(str.strip(line[3])))

                    elif dtype[k] is bool:
                        if str.strip(line[3]).upper() in ["TRUE", "FALSE"]:
                            IniParams[k] = str.strip(line[3]).upper() == "TRUE"
                        else:
                            raise TypeError("Control file line {0} must be True or False".format(line[0]))

                    elif dtype[k] is float:
                        IniParams[k] = float(str.strip(line[3]))

        # Make dates into seconds since UTC epoch
        IniParams["datastart"] = timegm(strptime(IniParams["datastart"] + " 00:00:00", "%m/%d/%Y %H:%M:%S"))
        IniParams["dataend"] = timegm(strptime(IniParams["dataend"], "%m/%d/%Y")) + 86400

        if IniParams["modelstart"] is None:
            IniParams["modelstart"] = IniParams["datastart"]
        else:
            IniParams["modelstart"] = timegm(strptime(IniParams["modelstart"] + " 00:00:00", "%m/%d/%Y %H:%M:%S"))

        if IniParams["modelend"] is None:
            IniParams["modelend"] = IniParams["dataend"]
        else:
            IniParams["modelend"] = timegm(strptime(IniParams["modelend"], "%m/%d/%Y")) + 86400

        IniParams["flushtimestart"] = IniParams["modelstart"] - IniParams["flushdays"] * 86400

        # If the number of transverse samples per direction 
        # is NOT reported, assume 4 (old default)
        if not IniParams["transsample_count"]:
            IniParams["transsample_count"] = 4.0

            # Format for heat source 8 methods same
        # as 8 directions but no north
        if IniParams["heatsource8"]:
            IniParams["trans_count"] = 7

        # Set the total number landcover sample count (0 = emergent)
        IniParams["sample_count"] = int(IniParams["transsample_count"]
                                        * IniParams["trans_count"])

        # Set up our evaporation method
        IniParams["penman"] = False
        if IniParams["calcevap"]:
            IniParams["penman"] = True if IniParams["evapmethod"] == "Penman" else False

        # make sure that the timestep divides into 60 minutes, 
        # or we may not land squarely on each hour's starting point.
        if float(60) / IniParams["dt"] - int(float(60) / IniParams["dt"]) > 1e-7:
            raise ValueError(
                "I'm sorry, your timestep ({0}) must evenly divide into 60 minutes.".format(IniParams["dt"]))
        else:
            # make dt measured in seconds
            IniParams["dt"] = IniParams["dt"] * 60

        # Make sure the output directory ends in a 
        # slash based on system platform
        if (platform.system() == "Windows" and IniParams["outputdir"][-1] != "\\"):
            raise ValueError("Output directory needs to have a backslash at end of the path. ..\\outputfolder\\")

        if (platform.system() == "Darwin" and IniParams["outputdir"][-1] != "/"):
            raise ValueError("Output directory needs to have a forward slash at the end of the path. ../outputfolder/")

            # the distance step must be an exact, greater or equal to one,
        # multiple of the sample rate.
        if (IniParams["dx"] % IniParams["longsample"]
                or IniParams["dx"] < IniParams["longsample"]):
            raise ValueError("Distance step (dx) must be a multiple of the longitudinal stream sample distance")

    def control_file_dict(self, **kwargs):
        """
        Returns a dictionary of the control file lines
        with empty values. Dictionary keys correspond to the
        keys in IniParams.
        """

        cf_dict = {"usertxt": [2, "Model Description/User Notes", "usertxt", None],
                   "name": [3, "Simulation Name", "name", None],
                   "inputdir": [4, "Input Directory Path", "inputdir", None],
                   "outputdir": [5, "Output Directory Path", "outputdir", None],
                   "length": [6, "Stream Length (kilometers)", "length", None],
                   "outputkm": [7, "Output Stream Kilometers", "outputkm", None],
                   "datastart": [8, "Data Start Date (mm/dd/yyyy)", "datastart", None],
                   "modelstart": [9, "Modeling Start Date (mm/dd/yyyy)", "modelstart", None],
                   "modelend": [10, "Modeling End Date (mm/dd/yyyy)", "modelend", None],
                   "dataend": [11, "Data End Date (mm/dd/yyyy)", "dataend", None],
                   "flushdays": [12, "Flush Initial Condition (days)", "flushdays", None],
                   "offset": [13, "Time Offset From UTC (hours)", "offset", None],
                   "dt": [14, "Model Time Step (minutes)", "dt", None],
                   "dx": [15, "Model Distance Step (meters)", "dx", None],
                   "longsample": [16, "Longitudinal Stream Sample Distance (meters)", "longsample", None],
                   "bcfile": [17, "Boundary Condition Input File Name", "bcfile", None],
                   "inflowsites": [18, "Tributary Inflow Sites", "inflowsites", None],
                   "inflowinfiles": [19, "Tributary Inflow Input File Name", "inflowinfiles", None],
                   "inflowkm": [20, "Tributary Inflow Model kilometers", "inflowkm", None],
                   "accretionfile": [21, "Accretion Input File Name", "accretionfile", None],
                   "metsites": [22, "Meteorological Data Sites", "metsites", None],
                   "metfiles": [23, "Meteorological Data Input File Name", "metfiles", None],
                   "metkm": [24, "Meteorological Data Model kilometers", "metkm", None],
                   "calcevap": [25, "Include Evaporation Losses From Flow (True/False)", "calcevap", None],
                   "evapmethod": [26, "Evaporation Method (Mass Transfer/Penman)", "evapmethod", None],
                   "wind_a": [27, "Wind Function Coefficient a", "wind_a", None],
                   "wind_b": [28, "Wind Function Coefficient b", "wind_b", None],
                   "calcalluvium": [29, "Include Deep Alluvium Temperature (True/False)", "calcalluvium", None],
                   "alluviumtemp": [30, "Deep Alluvium Temperature (Celsius)", "alluviumtemp", None],
                   "morphfile": [31, "Morphology Input Data File Name", "morphfile", None],
                   "lcdatafile": [32, "Land Cover Input Data File Name", "lcdatafile", None],
                   "lccodefile": [33, "Land Cover Codes Input File Name", "lccodefile", None],
                   "trans_count": [34, "Number Of Transects Per Node", "trans_count", None],
                   "transsample_count": [35, "Number Of Samples Per Transect", "transsample_count", None],
                   "transsample_distance": [36, "Distance Between Transect Samples (meters)", "transsample_distance",
                                            None],
                   "emergent": [37, "Account For Emergent Veg Shading (True/False)", "emergent", None],
                   "lcdatainput": [38, "Land Cover Data Input Type (Codes/Values)", "lcdatainput", None],
                   "canopy_data": [39, "Canopy Data Type (LAI/CanopyCover)", "canopy_data", None],
                   "lcsampmethod": [40, "Land Cover Sample Method (point/zone)", "lcsampmethod", None],
                   "heatsource8": [41, "Use Heat Source 8 Land Cover Methods (True/False)", "heatsource8", None],
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
        """Returns k, the riparian extinction coefficient data."""
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

    def lookup_lccode(self, code, lookup_list):
        """
        Returns a value based on a range provided in the lookup list.
        If the value does not fall within the range None is returned.
        
        lccode: numeric or string value to lookup in the lookup list.
        
        lookup_list: A list of tuples with the form [(min code, max code, return value),].
        The return value is returned if min code <= lccode <= max code.
        If the code is a string the min and max are equal to the code.
        """

        if lookup_list is None or code is None:
            return None

        try:
            float(code)
        except:
            return None

            # sort the list
        lookup_list.sort()

        for row in lookup_list:
            if row[0] <= float(code) <= row[1]:
                return row[2]

        # Not in range, return None
        return None

    def dict2list(self, data, colnames, skiprows=0, skipcols=0):

        d2 = zip(*[[k] + data[k] for k in colnames])

        # skiprows
        d3 = d2[-(len(d2) - skiprows):]

        # skip cols
        d4 = [line[+skipcols:] for line in d3]

        return [list(i) for i in d4]

    def parameterize_cf(self, overwrite=False, use_timestamp=False, **kwargs):
        """
        Writes the control file. Any keyword arguments
        passed will be parameterized into the control file.
        """
        if use_timestamp:
            cf_name = self.timestamp(self.control_file)
        else:
            cf_name = self.control_file

        # check to see if the file exists 
        if not overwrite and isfile(join(self.model_dir, cf_name)):
            msg = "HeatSource_Control.csv already exists. It will not be overwritten"
            logger.warning(msg)
            print_console(msg)
            return

        msg = "Writing control file"
        logger.info(msg)
        print_console(msg)
        cf_dict = self.control_file_dict()

        for k, v in kwargs.items():
            cf_dict[k][3] = v

        # sort the list is in the order of the line number
        cf_sorted = sorted(cf_dict.items(), key=itemgetter(1))
        cf_list = [line[1] for line in cf_sorted]

        self.write_to_output(self.model_dir, cf_name,
                             cf_list, self.headers_cf())

    def timestamp(self, name=""):
        """Returns the name with timestamp prefixed. Name must be a string.

        E.g. 20200314_031415_name
        """
        now = datetime.now()
        timestamp = now.strftime("%Y%m%d_%H%M%S")
        new_name = timestamp + "_" + name

        return new_name

    def parameterize_lccodes(self, lccodes=None, code_as_ht=False,
                             ht_list=None, can_list=None,
                             lai_list=None, k_list=None, oh_list=None,
                             overwrite=False):
        """
        Parameterize and writes the land cover codes file.
        
        lccodes: Optional list of tuples with the landcover code
        information. List takes the form:
        [(landcover name, code, ht, canopy, overhang),] or
        [(landcover name, code, ht, lai, k, overhang),]
        
        If lccodes=None the lcdata file identified in
        the control file will be read instead to identify the unique codes.
        
        The other input parameters including height (ht), canopy (can),
        lai, k, and overhang (oh) are parameterized using a range
        provided in list of tuples for each parameter.
        
        The parameter lookup list must take form
        [(min code, max code, return value),]
        where the return value is returned if min code <= lccode <= max code,
        
        code_as_ht: if True uses the landcover code as the height when
        the code is not found within a range of values passed to ht_list
        or when the ht_list is None.
        
        """
        msg = "Parameterize lccodes"
        logger.info(msg)
        print_console(msg)

        # check to see if the file exists 
        if not overwrite and isfile(join(self.model_dir, IniParams["lccodefile"])):
            msg = "{0} already exists. It will not be overwritten".format(IniParams["lccodefile"])
            logger.warning(msg)
            print_console(msg)
            return

        if lccodes is None:
            # read the land cover data
            lcdata = self.import_lcdata(return_list=True,
                                        skiprows=1, skipcols=0)

            codes = set()

            for row in lcdata:
                for col in range(8, (IniParams["trans_count"] *
                                     IniParams["transsample_count"]) + 9):
                    codes.add(str(row[col]))

            codes = list(codes)
            if code_as_ht:
                codes.sort(key=float)
            else:
                codes.sort()

            lccodes = []

            if IniParams["canopy_data"] == "LAI":
                for code in codes:
                    ht = self.lookup_lccode(code, ht_list)
                    if code_as_ht and ht is None:
                        ht = float(code)
                    lai = self.lookup_lccode(code, lai_list)
                    k = self.lookup_lccode(code, k_list)
                    oh = self.lookup_lccode(code, oh_list)
                    lccodes.append([None, code, ht, lai, k, oh])

            else:
                for code in codes:
                    ht = self.lookup_lccode(code, ht_list)
                    if code_as_ht and ht is None:
                        ht = float(code)
                    can = self.lookup_lccode(code, can_list)
                    oh = self.lookup_lccode(code, oh_list)
                    lccodes.append([None, code, ht, can, oh])

                    # lccodes = [[None, code, float(code), lai, k, oh] for code in codes]
                # else:
                # lccodes = [[None, code, ht, lai, k, oh] for code in codes]
            # else:
            # if code_as_ht:
            # lccodes = [[None, code, float(code), can, oh] for code in codes]
            # else:
            # lccodes = [[None, code, ht, can, oh] for code in codes]

        self.write_to_output(IniParams["inputdir"],
                             IniParams["lccodefile"],
                             lccodes, self.headers_lccodes())

    def read_csv_to_dict(self, inputdir, filename, colnames):
        """
        Reads a comma delimited text file and returns the data
        as a dictionary with the column header as the key.
        """

        # each value in each column is appended to a list
        data = defaultdict(list)

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
                    # if there is more than one cl
                    data[k].append(v)

        return self.validate(data)

    def read_csv_to_list(self, inputdir, filenames, skiprows, skipcols):
        """Reads a comma delimited text file into a
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
                newfile = [row for row in csv.reader(file_object.read().splitlines(), dialect="excel")]
                numcols = len(newfile[0])

            # skip rows    
            newfile = newfile[-(len(newfile) - skiprows):]
            # skip cols
            newfile = [line[+skipcols:] for line in newfile]
            if i == 1:
                data = newfile
            else:
                data = [data[row] + newfile[row] for row in range(0, len(newfile))]
            i = i + 1
        return data

    def setup_csv(self, use_timestamp=True, overwrite=True):
        """
        Writes blank input csv files based on settings in the control file
        """

        msg = "Starting input file setup"
        logging.info(msg)
        print_console(msg)

        # check if the input/output directories exist and 
        # create them if not
        if not exists(IniParams["inputdir"]):
            makedirs(IniParams["inputdir"])

        if not exists(IniParams["outputdir"]):
            makedirs(IniParams["outputdir"])

        now = datetime.now()
        timestamp = self.timestamp()
        timelist = self.datetime_string()
        kmlist = self.stream_kms()

        # For inflow and met data there can be a single input 
        # file for each node OR just one input file with multiple 
        # nodes in the file, This creates a list of the trib 
        # and met file names if there is more than one   
        if IniParams["inflowsites"] > 0:
            tribfiles = IniParams["inflowinfiles"].split(",")
        metfiles = IniParams["metfiles"].split(",")

        acclist = [[None, None, km] + [None] * 3 for km in kmlist]
        bclist = [[t, None, None] for t in timelist]
        lcdatalist = [[None, None, km] + [None] *
                      (len(self.headers_lcdata()) - 3) for km in kmlist]
        morphlist = [[None, None, km] + [None] * 10 for km in kmlist]

        # This writes to csv using the file name from the control 
        # file and adds a timestamp

        if use_timestamp:
            acc_file = timestamp + IniParams["accretionfile"]
            bc_file = timestamp + IniParams["bcfile"]
            lcdata_file = timestamp + IniParams["lcdatafile"]
            if IniParams["lcdatainput"] == "Codes":
                lccodes_file = timestamp + IniParams["lccodefile"]
            morph_file = timestamp + IniParams["morphfile"]
        else:
            # Name the files as they appear in the control file
            acc_file = IniParams["accretionfile"]
            bc_file = IniParams["bcfile"]
            lcdata_file = IniParams["lcdatafile"]
            if IniParams["lcdatainput"] == "Codes":
                lccodes_file = IniParams["lccodefile"]
            morph_file = IniParams["morphfile"]

        msg = "Writing empty csv files"
        logging.info(msg)
        print_console(msg)

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

        for file in metfiles:
            metlist = [[t] + [None] * 4 * int((IniParams["metsites"] / len(metfiles))) for t in timelist]

            if use_timestamp:
                met_filename = timestamp + file.strip()
            else:
                met_filename = file.strip()

            if overwrite:
                self.write_to_csv(IniParams["inputdir"],
                                  met_filename,
                                  metlist, self.headers_met())
            else:
                if not isfile(join(IniParams["inputdir"], met_filename)):
                    self.write_to_csv(IniParams["inputdir"],
                                      met_filename,
                                      metlist, self.headers_met())

        if IniParams["inflowsites"] > 0:
            for file in tribfiles:
                inflowlist = [[t] + [None] * 2 * int((IniParams["inflowsites"] / len(tribfiles))) for t in timelist]

                if use_timestamp:
                    trib_filename = timestamp + file.strip()
                else:
                    trib_filename = file.strip()

                if overwrite:
                    self.write_to_csv(IniParams["inputdir"], trib_filename,
                                      inflowlist, self.headers_inflow())
                else:
                    if not isfile(join(IniParams["inputdir"], trib_filename)):
                        self.write_to_csv(IniParams["inputdir"], trib_filename,
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
        timelist = range(IniParams["datastart"], IniParams["dataend"] + 60, 3600)
        for i, val in enumerate(timelist):
            timelist[i] = strftime("%m/%d/%Y %H:%M", gmtime(val))
        return timelist

    def stream_kms(self):
        """Return a list of stream kilometers sorted from
        headwaters to mouth"""

        # Since the stream length is formatted as a floating point number 
        # we can sometimes run into trouble when calculating the number 
        # of nodes with dx due to ceil function rounding up on 
        # floating points that are not exact representations of the 
        # input value. Therefore we enforce a precision only up to the 
        # ten thousandths decimal place to make sure the number of nodes
        # is correct.        
        num_nodes = int(ceil(round(IniParams["length"] * 1000 / (IniParams["longsample"]), 4))) + 1
        kmlist = []
        kmlist = [(float(node) * IniParams["longsample"]) / 1000 for node in range(0, num_nodes)]
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

            # tranlsate removes numbers from the key, e.g. TEMPERATURE2
            if key.translate(None, digits) not in dtype.keys():
                # This is to find the correct landcover data 
                # key since they are all different
                if "LC" in key:
                    # basestring
                    k = "LC"
                elif any(s in key for s in ["HT", "ELE", "LAI", "k", "CAN", "OH"]):
                    # float
                    k = "ELE"
            else:
                k = key.translate(None, digits)

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
        """Write the outlist to a comma delimited text file
        colnames are optional.
        """

        # split by comma if multiple files
        filenames = filenames.split(",")

        filenames = [filenames] if not isinstance(filenames, list) else filenames

        for filename in filenames:
            with open(join(outputdir, filename), "wb") as file_object:
                writer = csv.writer(file_object, dialect="excel")
                if colnames:
                    writer.writerow(colnames)
                writer.writerows(outlist)

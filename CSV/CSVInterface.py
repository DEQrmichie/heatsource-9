# Heat Source, Copyright (C) 2000-2014, Oregon Department of Environmental Quality

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Main interface class for CSV->Python conversion

CSVInterface provides the single resource for converting the data
in the CSV files into a list of StreamNode classes for use
by the Heat Source model.
"""
# Builtin methods
from __future__ import division, print_function
from itertools import ifilter, izip, chain, repeat, count
from math import ceil, log, degrees, atan
from os.path import exists, join, normpath
from os import unlink
from bisect import bisect
from time import strptime, ctime, strftime, gmtime
from calendar import timegm
from datetime import datetime
import platform
import pandas as pd

# Heat Source Methods
from ..Dieties.IniParamsDiety import IniParams
from ..Stream.StreamNode import StreamNode
from ..Utils.Dictionaries import Interpolator
from ..Utils.easygui import buttonbox

class CSVInterface(object):
    """Reads the Heat Source input CSV files, creates a list of StreamNode instances, and populates those StreamNodes"""

    def __init__(self, inputdir, control_file, log=None, run_type=0):
        self.run_type = run_type
        self.log = log
        self.Reach = {}
        
        self.ReadControlFile(inputdir, control_file)
    
        if run_type == 3:
            # Setup input files
            self.SetupInputFiles(inputdir, control_file)
        else:
            # Get the list of model periods times
            self.flowtimelist = self.GetTimelistUNIX()
            self.continuoustimelist = self.GetTimelistUNIX()
            self.flushtimelist = self.GetTimelistFlushPeriod()         
        
            # Start through the steps of building a reach full of StreamNodes
            self.GetBoundaryConditions()
            self.BuildNodes()
            if IniParams["lidar"]: self.BuildZonesLidar()
            else: self.BuildZonesNormal()
            self.GetTributaryData()
            self.GetClimateData()
            self.SetAtmosphericData()
            self.OrientNodes()    

    def ReadControlFile(self, inputdir, control_file):
        """Reads the initialization parameters from the control file into the Initialization dictionary {IniParams} and does some formatting"""
    
        if exists(join(inputdir,control_file)) == False:
            raise Exception("HeatSource_Control.csv not found. Move the executable or place the control file in this directory: %s." % inputdir)    
        
        print("Reading Control File")
    
        #TODO RM fix so dict is related to numbers or text symbols
        lst = {"usertxt": "# USER TEXT",
                    "name": "# SIMULATION NAME",
                    "length": "# STREAM LENGTH (KILOMETERS)",
                    "outputdir": "# OUTPUT PATH",
                    "inputdir": "# INPUT PATH",
                    "date": "# DATA START DATE (mm/dd/yyyy)",
                    "modelstart": "# MODELING START DATE (mm/dd/yyyy)",
                    "modelend": "# MODELING END DATE (mm/dd/yyyy)",
                    "end": "# DATA END DATE (mm/dd/yyyy)",
                    "flushdays": "# FLUSH INITIAL CONDITION (DAYS)",
                    "offset": "# TIME OFFSET FROM UTC (HOURS)",
                    "dt": "# TIME STEP - DT (MIN)",
                    "dx": "# DISTANCE STEP - DX (METERS)",
                    "longsample": "# LONGITUDINAL SAMPLE RATE (METERS)",
                    "bcfile": "# BOUNDARY CONDITION FILE NAME",
                    "inflowsites": "# TRIBUTARY SITES",
                    "inflowinfiles": "# TRIBUTARY INPUT FILE NAMES",
                    "inflowkm": "# TRIBUTARY MODEL KM",
                    "accretionfile": "# ACCRETION INPUT FILE NAME",
                    "contsites": "# CLIMATE DATA SITES",
                    "contfiles": "# CLIMATE INPUT FILE NAMES",
                    "contkm": "# CLIMATE MODEL KM",
                    "calcevap": "# INCLUDE EVAPORATION LOSSES FROM FLOW (TRUE/FALSE)",
                    "evapmethod": "# EVAPORATION METHOD (Mass Transfer/Penman)",
                    "wind_a": "# WIND FUNCTION COEFFICIENT A",
                    "wind_b": "# WIND FUNCTION COEFFICIENT B",
                    "calcalluvium": "# INCLUDE DEEP ALLUVIUM TEMPERATURE (TRUE/FALSE)",
                    "alluviumtemp": "# DEEP ALLUVIUM TEMPERATURE (*C)",
                    "morphfile": "# MORPHOLOGY DATA FILE NAME",
                    "lcdatafile": "# LANDCOVER DATA FILE NAME",
                    "lccodefile": "# LANDCOVER CODES FILE NAME",
                    "transsample": "# TRANSVERSE SAMPLE RATE (METERS)",
                    "transsample_count": "# NUMBER OF TRANSVERSE SAMPLES IN EACH DIRECTION",
                    "radialsample_count": "# NUMBER OF RADIAL SAMPLES",
                    "emergent": "# ACCOUNT FOR EMERGENT VEG SHADING (TRUE/FALSE)",
                    "lidar": "# LIDAR DATA USED FOR VEG CODES (TRUE/FALSE)",
                    "lcdensity": "# CANOPY COVER VALUE FOR LIDAR DATA",
                    "lcoverhang": "# LANDCOVER STREAM OVERHANG FOR LIDAR DATA (METERS)",
                    "beers_data": "# BEER'S LAW INPUT DATA TYPE (LAI/CanopyCover)",
                    "vegDistMethod": "# VEGETATION ANGLE CALCULATION METHOD (point/zone)",
                    "heatsource8": "# USE HEAT SOURCE 8 LANDCOVER METHODS (TRUE/FALSE)",}
            
        cf = pd.read_csv(join(inputdir,control_file),quotechar='"',quoting=0,header=None,na_values=None)
        # TODO put a checker in here
        #raise Exception("Control file does not have %s" % v) 
        for k,v in lst.iteritems():
            for i in range(0,len(cf)):
                if (cf.ix[i,1] == v):
                    if (self.isNum(cf.ix[i,2]) == True):
                        IniParams[k] = float(str.strip(cf.ix[i,2]))
                    else:
                        IniParams[k] = str.strip(cf.ix[i,2])
        del(cf)
        
        # These might be blank, make them zeros # TODO check this works with new 9.0.0 implementation
        for key in ["inflowsites","flushdays","wind_a","wind_b"]:
            IniParams[key] = 0.0 if not IniParams[key] else IniParams[key]
        
        # Convert string TRUE/FALSE to bool. Note this defaults to false if string is not formated exactly as 'TRUE'
        IniParams["calcalluvium"] = IniParams["calcalluvium"] in ("TRUE")
        IniParams["calcevap"] = IniParams["calcevap"] in ("TRUE")
        IniParams["emergent"] = IniParams["emergent"] in ("TRUE")
        IniParams["lidar"] = IniParams["lidar"] in ("TRUE")
        IniParams["heatsource8"] = IniParams["heatsource8"] in ("TRUE")
        
        # If the number of transverse sample per direction is NOT reported, assume 4 (old default)
        IniParams["transsample_count"] = 4.0 if not IniParams["transsample_count"] else IniParams["transsample_count"]
        
        # If True use heat source 8 default, same as 8 directions but no north
        if IniParams["heatsource8"] == True:
            IniParams["radialsample_count"] = 999 
        else:
            IniParams["radialsample_count"]
        
        # Set the total number landcover sample count (0 = emergent)
        if IniParams["heatsource8"] == True:
            IniParams["sample_count"] = int(IniParams["transsample_count"] * 7)
        else:
            IniParams["sample_count"] = int(IniParams["transsample_count"] * IniParams["radialsample_count"])
        
        # Then make all of these integers because they're used later in for loops
        for key in ["inflowsites","flushdays","contsites", "transsample_count","radialsample_count"]:
            IniParams[key] = int(IniParams[key])
       
        # These need to be strings
        for key in ["inflowkm","contkm"]:
            IniParams[key] = str(IniParams[key])
        
        # Set up our evaporation method
        IniParams["penman"] = False
        if IniParams["calcevap"]:
            IniParams["penman"] = True if IniParams["evapmethod"] == "Penman" else False
           
        # The offset should be negated to work around issues with internal date
        # representation. i.e. Pacific time is -7 from UTC, but the code needs a +7 to work.
        # TODO: This is probably a bug in ChronosDiety, not the time module.
        IniParams["offset"] = -1 * IniParams["offset"]
        
        # Make the dates into datetime instances of the start/stop dates
        IniParams["date"] = timegm(strptime(IniParams["date"] + " 00:00:00" ,"%m/%d/%Y %H:%M:%S"))
        IniParams["end"] = timegm(strptime(IniParams["end"]+ " 23:59:59","%m/%d/%Y %H:%M:%S"))
        
        if IniParams["modelstart"] is None:
            IniParams["modelstart"] = IniParams["date"]
        else:
            IniParams["modelstart"] = timegm(strptime(IniParams["modelstart"] + " 00:00:00","%m/%d/%Y %H:%M:%S"))
        
        if IniParams["modelend"] is None:
            IniParams["modelend"] = IniParams["end"]
        else:
            IniParams["modelend"] = timegm(strptime(IniParams["modelend"] + " 23:59:59","%m/%d/%Y %H:%M:%S"))
            
        IniParams["flushtimestart"] = IniParams["modelstart"] - IniParams["flushdays"]*86400
        
        # make sure alluvium temp is present and a floating point number.
        IniParams["alluviumtemp"] = 0.0 if not IniParams["alluviumtemp"] else float(IniParams["alluviumtemp"])
           
        # make sure that the timestep divides into 60 minutes, or we may not land squarely on each hour's starting point.
        #if 60%IniParams["dt"] > 1e-7:
        if float(60)/IniParams["dt"] - int(float(60)/IniParams["dt"]) > 1e-7:
            raise Exception("I'm sorry, your timestep (%0.2f) must evenly divide into 60 minutes." % IniParams["dt"])
        else:
            IniParams["dt"] = IniParams["dt"]*60 # make dt measured in seconds
         
        # Make sure the output directory ends in a slash based on system platform
        if (platform.system() == 'Windows' and IniParams["outputdir"][-1] != "\\"):
            raise Exception("Output directory needs to have a backslash at end of the path. ..\\outputfolder\\")
        
        if (platform.system() == 'Darwin' and IniParams["outputdir"][-1] != "/"):
            raise Exception("Output directory needs to have a forward slash at the end of the path. ../outputfolder/")    
        
        # Set up the log file in the outputdir
        self.log.SetFile(normpath(join(IniParams["outputdir"],"outfile.log")))
        
        # Make empty Dictionaries for the boundary conditions
        self.Q_bc = Interpolator()
        self.T_bc = Interpolator()
        
        # List of kilometers with continuous data nodes assigned.
        self.ContDataSites = []
        
        # the distance step must be an exact, greater or equal to one, multiple of the sample rate.
        if (IniParams["dx"]%IniParams["longsample"]
            or IniParams["dx"]<IniParams["longsample"]):
            raise Exception("Distance step must be a multiple of the Longitudinal transfer rate")
        # Some convenience variables
        self.dx = IniParams["dx"]
        self.multiple = int(self.dx/IniParams["longsample"]) #We have this many samples per distance step            
    
    def SetupInputFiles(self, inputdir, control_file):
        """Formats and writes blank input files based on settings in the control file"""
        
        print("Starting input file setup")
        
        now = datetime.now()    
        timestamp = now.strftime("%Y%m%d_%H%M%S")
        timelist = self.GetTimeListString()
        kmlist = self.GetStreamKMlist()
        
        # Setup the headers for the Landcover file
        if IniParams["beers_data"] == "LAI":  #Use LAI methods
            type = ['LC','ELE','LAI','k']
            emergentlabel ='LAI_EMERGENT'
        else:        
            type = ['LC','ELE','CCV']
            emergentlabel = 'CCV_EMERGENT'
        
        lcDataColHeaders =['Longitude','Latitude','TopoWest','TopoSouth','TopoEast','EmergentVeg']      
        if IniParams["heatsource8"] == True: # a flag indicating the model should use the heat source 8 methods (same as 8 directions but no north)
            dir = ['NE','E','SE','S','SW','W','NW']
        else:        
            dir = ['D' + str(x) for x in range(1,IniParams["radialsample_count"]+ 1)]
        
        zone = range(1,int(IniParams["transsample_count"])+1)
        
        # Concatenate the type, dir, and zone and order in the correct way
        for t in range(0,len(type)):
            for d in range(0,len(dir)):
                for z in range(0,len(zone)):
                    if t==2 and d==0 and z==0:
                        lcDataColHeaders.append(emergentlabel) # add this when we begin the LAI/Canopy cover series
                        lcDataColHeaders.append(type[t]+'_'+dir[d]+'_'+str(zone[z]))
                    else:
                        lcDataColHeaders.append(type[t]+'_'+dir[d]+'_'+str(zone[z]))
        
        # creates a list of the file names if there is more than one   
        if IniParams["inflowsites"] > 0:
            tribfiles = IniParams["inflowinfiles"].split(",")
        climatefiles = IniParams["contfiles"].split(",")
        
        # This sets the header names and index values (either the stream km or dates))
        # For inflow and climate there can be a single input file for each node OR just one input file with multiple nodes in the file
        # This just repeats the header if needed based on the number of sites and files.
        climatefile = pd.DataFrame(index=timelist, columns=['Cloudiness', 'WindSpeed','RelativeHumidity','AirTemp']*int((IniParams["contsites"]/len(climatefiles))))
        if IniParams["inflowsites"] > 0:
            inflowfile = pd.DataFrame(index=timelist, columns=['Flow','Temp']*int((IniParams["inflowsites"]/len(tribfiles))))
        
        # This sets the header names and index values for the other input files
        bcfile = pd.DataFrame(index=timelist, columns=['Flow', 'Temp'])
        if IniParams["beers_data"] == "LAI":  #Use LAI methods
            lccodes = pd.DataFrame(columns=['Name','Code','Height','LAI','k','Overhang'])
        else:
            lccodes = pd.DataFrame(columns=['Name','Code','Height','CanopyCover','k','Overhang'])
        lcfile = pd.DataFrame(index=kmlist,columns=lcDataColHeaders)
        accfile= pd.DataFrame(index=kmlist,columns=['Inflow','Temp','Outflow'])
        morphfile = pd.DataFrame(index=kmlist,columns=['Elevation','Gradient','BottomWidth','ChannelAngleZ','Mannings_n','SedThermalConductivity','SedThermalDiffusivity','SedHyporheicThickness','%HyporheicExchange','Porosity'])
        
        # This writes to csv using the file name from the control file and adds a timestamp
        print("Writing empty csv files")
        bcfile.to_csv(join(IniParams["inputdir"],'input_'+timestamp+'_'+IniParams["bcfile"]), sep=',', na_rep='NA',index_label='DateTime')
        if IniParams["lidar"] == False:
            lccodes.to_csv(join(IniParams["inputdir"],'input_'+timestamp+'_'+IniParams["lccodefile"]), sep=',', na_rep='NA',index_label='Name')
        else:
            print("...LiDAR = TRUE. land cover codes file not written")
        lcfile.to_csv(join(IniParams["inputdir"],'input_'+timestamp+'_'+IniParams["lcdatafile"]), sep=',', na_rep='NA',index_label='Stream_KM')
        accfile.to_csv(join(IniParams["inputdir"],'input_'+timestamp+'_'+IniParams["accretionfile"]), sep=',', na_rep='NA',index_label='Stream_KM')
        morphfile.to_csv(join(IniParams["inputdir"],'input_'+timestamp+'_'+IniParams["morphfile"]), sep=',', na_rep='NA',index_label='Stream_KM')
        for file in climatefiles:
            climatefile.to_csv(join(IniParams["inputdir"],'input_'+timestamp+'_'+file.strip()), sep=',', na_rep='NA',index_label='DateTime')
        
        if IniParams["inflowsites"] > 0:
            for file in tribfiles:
                inflowfile.to_csv(join(IniParams["inputdir"],'input_'+timestamp+'_'+file.strip()), sep=',', na_rep='NA',index_label='DateTime')
        else:
            print("Inflow Sites = 0 , inflow file not written")
        
        print("Finished input file setup")
    
    def isNum(self,v):
        try:
            inNumberint = int(v)
            isInt = True
        except ValueError:
            isInt = False
            pass
        try:
            inNumberfloat = float(v)
            isFloat = True
        except ValueError:
            isFloat = False
            pass
        if(isInt == True or isFloat == True):
            return True
        else:
            return False

    def OrientNodes(self):
        print("Initializing StreamNodes")
        # Now we manually set each nodes next and previous kilometer values by stepping through the reach
        l = sorted(self.Reach.keys(), reverse=True)
        head = self.Reach[max(l)] # The headwater node
        # Set the previous and next kilometer of each node.
        slope_problems = []
        for i in xrange(len(l)):
            key = l[i] # The current node's key
            # Then, set pointers to the next and previous nodes
            if i == 0: pass
            else: self.Reach[key].prev_km = self.Reach[l[i-1]] # At first node, there's no previous
            try:
                self.Reach[key].next_km = self.Reach[l[i+1]]
            except IndexError:
            ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ## For last node (mouth) we set the downstream node equal to self, this is because
            ## we want to access the node's temp if there's no downstream, and this safes us an
            ## if statement.
                self.Reach[key].next_km = self.Reach[key]
            # Set a headwater node
            self.Reach[key].head = head
            self.Reach[key].Initialize()
            # check for a zero slope. We store all of them before checking so we can print a lengthy error that no-one will ever read.
            if self.Reach[key].S <= 0.0: slope_problems.append(key)
        if self.run_type != 1: # zeros are alright in shade calculations
            if len(slope_problems):
                raise Exception ("The following reaches have zero slope. Kilometers: %s" %",".join(['%0.3f'%i for i in slope_problems]))

    def close(self):
        del self.T_bc, self.Reach

    def CheckEarlyQuit(self): #TODO Need to fix this RM
        """Checks a value to see whether the user wants to stop the model before we completely set everything up"""
        if exists("c:\\quit_heatsource"):
            unlink("c:\\quit_heatsource")
            self.QuitMessage()
        
    def SetAtmosphericData(self):
        """For each node without climate 'continuous' data, use closest (up or downstream) node's data"""
        self.CheckEarlyQuit()
        print("Setting Atmospheric Data")
        sites = self.ContDataSites # Localize the variable for speed
        sites.sort() #Sort the climate site by km. This is necessary for the bisect module
        c = count()
        l = self.Reach.keys()
        # This routine iterates through all nodes and uses bisect to determine which climate site is closest to the node and 
        # initializes that node with the climate data that is closest (up or downstream)
        for km, node in self.Reach.iteritems():
            if km not in sites: # if we are not on the node where the climate data is assigned we have t
                # Kilometer's downstream and upstream
                lower = bisect(sites,km)-1 if bisect(sites,km)-1 > 0 else 0 # zero is the lowest (protect against value of -1)
                # bisect returns the length of a list when the bisecting number is greater than the greatest value.
                # Here we protect by max-ing out at the length of the list.
                upper = min([bisect(sites,km),len(sites)-1])
                # Use the indexes to get the kilometers from the sites list
                down = sites[lower]
                up = sites[upper]
                datasite = self.Reach[up] # Initialize to upstream's continuous data
                if km-down < up-km: # Only if the distance to the downstream node is closer do we use that
                    datasite = self.Reach[down]
                self.Reach[km].ContData = datasite.ContData
                print("Setting Atmospheric Data", c.next()+1, len(l))

    def GetBoundaryConditions(self):
        """Get the boundary conditions"""
        self.CheckEarlyQuit()
        # Get the columns, which is faster than accessing cells
        print("Reading boundary conditions")
        timelist = self.continuoustimelist
        
        # the data block is a tuple of tuples, each corresponding to a timestamp.
        data = pd.read_csv(join(IniParams["inputdir"],IniParams["bcfile"]),quotechar='"',quoting=0, index_col='DateTime')
        data.columns = ['Flow','Temp']
        
        # Check out GetTributaryData() for details on this reformatting of the data
        # for the progress bar
        length = len(data)
        c = count()
        # Now set the discharge and temperature boundary condition dictionaries.

        for i in xrange(len(timelist)):
            time = timelist[i]
            flow, temp = data.Flow.values[i],data.Temp.values[i]
            # Get the flow boundary condition
            if flow == 0 or not flow:
                if self.run_type != 1:
                    raise Exception("Missing flow boundary condition for day %s " % ctime(time))
                else: flow = 0
            self.Q_bc[time] = flow
            # Temperature boundary condition
            t_val = temp if temp is not None else 0.0
            self.T_bc[time] = t_val
            print("Reading boundary conditions",c.next(),length)

        # Next we expand or revise the dictionary to account for the flush period
        # Flush flow: model start value over entire flush period
        for i in xrange(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            self.Q_bc[time] = self.Q_bc[IniParams["modelstart"]]
        # Flush temperature: first 24 hours repeated over flush period
        first_day_time = IniParams["modelstart"]
        second_day = IniParams["modelstart"] + 86400
        for i in xrange(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            self.T_bc[time] = self.T_bc[first_day_time]
            first_day_time += 3600
            if first_day_time >= second_day:
                first_day_time = IniParams["modelstart"]

        self.Q_bc = self.Q_bc.View(IniParams["flushtimestart"], IniParams["modelend"], aft=1)
        self.T_bc = self.T_bc.View(IniParams["flushtimestart"], IniParams["modelend"], aft=1)
    
    def GetLocations(self,ini):
        """Build a list of kilometers corresponding to the tributary inflow or climate data sites"""
        # ini names that are passed: "inflowkm" or "contkm"
        t = ()
        l = self.Reach.keys()
        l.sort()

        if ini == "contkm" or IniParams["inflowsites"] > 0:
            kms = IniParams[ini].split(",") # get a list of sites by km
            kms = tuple([float(line.strip()) for line in kms]) # remove spaces and make float
        else:
            kms= tuple([])    
    
        for site in xrange(0,len(kms)):
            km = kms[site]
            if km is None or not isinstance(km, float):
                # This is a bad dataset if there's no kilometer
                raise Exception("Must have a stream kilometer (e.g. 15.3) for each node in %s page!" % ini)
            key = bisect(l,km)-1
            t += l[key], # Index by kilometer
        return t    
    
    def GetStreamKMlist(self):
        """Build a list of stream kilometers sorted from headwaters to mouth"""
        print("Creating Stream KM list")
        num_nodes = int(ceil(IniParams["length"]*1000/(IniParams["longsample"]))) +1
        kmlist = []    
        kmlist = [(node * IniParams["longsample"])/1000 for node in range(0,num_nodes)]    
        kmlist.sort(reverse=True)
        return kmlist    
            
    def GetTimeListString(self):
        """Build a time list in string format (MM/DD/YYYY HH:MM) corresponding to the data start and end dates available in the control file"""
        print("Creating timelist")
        timelist = []        
        timelist = range(IniParams['date'],IniParams['end']+60,3600) # hourly timestep
        for i in range(0,len(timelist)):
            timelist[i] = strftime("%m/%d/%Y %H:%M",gmtime(timelist[i]))
        return timelist    
               
    def GetTimelistUNIX(self):
        """Build a UNIX time list of floating point time values corresponding to the data start and end dates available in the control file"""        
        timelist = []        
        timelist = range(IniParams['date'],IniParams['end']+60,3600) # hourly timestep
        return tuple(timelist)
                        
    def GetTimelistFlushPeriod(self):
        """Build a UNIX time list that represents the flushing period"""
        #This assumes that data is hourly, not tested with variable input timesteps
        flushtimelist = []
        flushtime = IniParams["flushtimestart"]
        while flushtime < IniParams["modelstart"]:
            flushtimelist += flushtime,
            flushtime += 3600
        return tuple(flushtimelist)

    def GetTributaryData(self):
        """Populate the tributary flow and temperature values for nodes from the Flow Data page"""
        self.CheckEarlyQuit()
        print("Reading inflow data")
        # Get a list of the timestamps that we have data for, and use that to grab the data block
        timelist = self.flowtimelist
        data = pd.DataFrame()
        
        if IniParams["inflowsites"] > 0:
            tribfiles = IniParams["inflowinfiles"].split(",")
            for file in tribfiles:
                newfile = pd.read_csv(join(IniParams["inputdir"],file.strip()),quotechar='"',quoting=0,index_col='DateTime')
                data = pd.concat([data,newfile],axis=1)
               
        # The data is being put into this format
        # | Site 1     | Site 2     | Site 3       | ...
        # [((0.3, 15.7), (0.3, 17.7), (0.02, 18.2)), (, ...))]
        # To facilitate each site having it's own two item tuple.
        # The calls to tuple() just ensure that we are not making lists, which can
        # be changed accidentally. Without them, the line is easier to understand
        data = [tuple(zip(line[0:None:2],line[1:None:2])) for line in data.values]
        
        # Get a tuple of kilometers to use as keys to the location of each tributary 
        kms = self.GetLocations("inflowkm")    
        
        length = len(timelist)
        tm = count() # Which datapoint time are we recording
        nodelist = [] # Quick list of nodes with flow data
        
        if IniParams["inflowsites"] > 0: 
            for time in timelist:
                line = data.pop(0)
                # Error checking?! Naw!!
                c = count()
                for flow, temp in line:
                    i = c.next()
                    node = self.Reach[kms[i]] # Index by kilometer
                    if node not in nodelist or not len(nodelist): nodelist.append(node)
                    if flow is None or (flow > 0 and temp is None):
                        raise Exception("Cannot have a tributary with blank flow or temperature conditions")
                    # Here, we actually set the tribs library, appending to a tuple. Q_ and T_tribs are
                    # tuples of values because we may have more than one input for a given node
                    node.Q_tribs[time] += flow, #Append to tuple
                    node.T_tribs[time] += temp,
                    print("Reading inflow data",tm.next()+1, length)

        # Next we expand or revise the dictionary to account for the flush period
        # Flush flow: model start value over entire flush period
        for i in xrange(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            for node in nodelist:
                node.Q_tribs[time] = node.Q_tribs[IniParams["modelstart"]]
        # Flush temperature: first 24 hours repeated over flush period
        first_day_time = IniParams["modelstart"]
        second_day = IniParams["modelstart"] + 86400
        for i in xrange(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            for node in nodelist:
                node.T_tribs[time] = node.T_tribs[first_day_time]
            first_day_time += 3600
            if first_day_time >= second_day:
                first_day_time = IniParams["modelstart"]

        # Now we strip out the unnecessary values from the dictionaries. This is placed here
        # at the end so we can dispose of it easily if necessary
        for node in nodelist:
            node.Q_tribs = node.Q_tribs.View(IniParams["flushtimestart"], IniParams["modelend"], aft=1)
            node.T_tribs = node.T_tribs.View(IniParams["flushtimestart"], IniParams["modelend"], aft=1)

    def GetClimateData(self):
        """Get data from the input climate data csv file"""
        # This is remarkably similar to GetInflowData. We get a block of data, then set the dictionary of the node
        self.CheckEarlyQuit()
        print("Reading Continuous Data")

        timelist = self.continuoustimelist
        
        climatefiles = IniParams["contfiles"].split(",")
        climatedata = pd.DataFrame()
        for file in climatefiles:
            newfile = pd.read_csv(join(IniParams["inputdir"],file.strip()),quotechar='"',quoting=0, index_col='DateTime')
            climatedata = pd.concat([climatedata,newfile],axis=1)
            
        data = [tuple(zip(line[0:None:4],line[1:None:4],line[2:None:4],line[3:None:4])) for line in climatedata.values]

        # Get a tuple of kilometers to use as keys to the location of each climate node
        kms = self.GetLocations("contkm") 
                
        tm = count() # Which datapoint time are we recording
        length = len(timelist)
        for time in timelist:
            line = data.pop(0)
            c = count()
            for cloud, wind, humid, air in line:
                i = c.next()
                node = self.Reach[kms[i]] # Index by kilometer
                # Append this node to a list of all nodes which have continuous data
                if node.km not in self.ContDataSites:
                    self.ContDataSites.append(node.km)
                # Perform some tests for data accuracy and validity
                if cloud is None: cloud = 0.0
                if wind is None: wind = 0.0
                if cloud < 0 or cloud > 1:
                    if self.run_type == 1: # Alright in shade-a-lator # TODO zeros should not get a pass in solar only runs, fix
                        cloud = 0.0
                    else: raise Exception("Cloudiness (value of '%s' in Continuous Data) must be greater than zero and less than one." % `cloud`) # TODO RM fix this so it gives the km in exception - actualy do for all
                if humid < 0 or humid is None or humid > 1:
                    if self.run_type == 1: # Alright in shade-a-lator
                        humid = 0.0
                    else: raise Exception("Humidity (value of '%s' in Continuous Data) must be greater than zero and less than one." % `humid`)
                if air is None or air < -90 or air > 58:
                    if self.run_type == 1: # Alright in shade-a-lator
                        air = 0.0
                    else: raise Exception("Air temperature input (value of '%s' in Continuous Data) outside of world records, -89 to 58 deg C." % `air`)
                node.ContData[time] = cloud, wind, humid, air
            print("Reading continuous data", tm.next()+1, length)

        # Flush meteorology: first 24 hours repeated over flush period
        first_day_time = IniParams["modelstart"]
        second_day = IniParams["modelstart"] + 86400
        for i in xrange(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            for km in self.ContDataSites:
                node = self.Reach[km]
                node.ContData[time] = node.ContData[first_day_time]
            first_day_time += 3600
            if first_day_time >= second_day:
                first_day_time = IniParams["modelstart"]

        # Now we strip out the climate data outside the model period from the dictionaries. This is placed here
        # at the end so we can dispose of it easily if necessary
        print("Subsetting the Continuous Data to model period")
        tm = count()
        length = len(self.ContDataSites)
        for km in self.ContDataSites:
            node = self.Reach[km]
            node.ContData = node.ContData.View(IniParams["flushtimestart"], IniParams["modelend"], aft=1)
            print("Subsetting ",tm.next()+1, length)
    
    def zipper(self,iterable,mul=2):
        """Zippify list by grouping <mul> consecutive elements together

        Zipper returns a list of lists where the internal lists are groups of <mul> consecutive
        elements from the input list. For example:
        >>> lst = [0,1,2,3,4,5,6,7,8,9]
        >>> zipper(lst)
        [[0],[1,2],[3,4][5,6],[7,8],[9]]
        The first element is a length 1 list because we assume that it is a single element node (headwaters).
        Note that the last element, 9, is alone as well, this method will figure out when there are not
        enough elements to make n equal length lists, and modify itself appropriately so that the remaining list
        will contain all leftover elements. The usefulness of this method is that it will allow us to average over each <mul> consecutive elements
        """
        # From itertools recipes... We use all but the first (boundary node) element
        lst = [i for i in izip(*[chain(iterable[1:], repeat(None, mul-1))]*mul)]
        # Then we tack on the boundary node element
        lst.insert(0,(iterable[0],))
        # Then strip off the None values from the last (if any)
        lst[-1] = tuple(ifilter(lambda x: x is not None,lst[-1]))
        return self.numify(lst)

    def numify(self, lst):
        """Take a list of iterables and remove all values of None or empty strings"""
        # Remove None values at the end of each individual list
        for i in xrange(len(lst)):
            # strip out values of None from the tuple, returning a new tuple
            lst[i] = [x for x in ifilter(lambda x: x is not None, lst[i])]
        # Remove blank strings from within the list
        for l in lst:
            n = []
            for i in xrange(len(l)):
                if l[i] == "": n.append(i)
            n.reverse()
            for i in n: del l[i]
        # Make sure there are no zero length lists because they'll fail if we average
        for i in xrange(len(lst)):
            if len(lst[i]) == 0: lst[i].append(0.0)
        return lst

    def multiplier(self, iterable, predicate=lambda x:x):
        """Return an iterable that was run through the zipper

        Take an iterable and strip the values of None, then send to the zipper
        and apply predicate to each value returned (zipper returns a list)"""
        # This is a way to safely apply a generic lambda function to an iterable.
        # If I were paying attention to design, instead of just hacking away, I would
        # have done this with decorators to modify the function. Now I'm too lazy to
        # re-write it (well, not lazy, but I'm not paid as a programmer, and so I have
        # "better" things to do than optimize our code.)
        # First we strip off the None values.
        stripNone = lambda y: [i for i in ifilter(lambda x: x is not None, y)]
        return [predicate(stripNone(x)) for x in self.zipper(iterable,self.multiple)]

    def zeroOutList(self, lst): # TODO Not used
        """Replace blank values in a list with zeros"""
        test = lambda x: 0.0 if x=="" else x
        return [test(i) for i in lst]

    def GetColumnarData(self):
        """return a dictionary of attributes that are averaged or summed as appropriate"""
        self.CheckEarlyQuit()
        # Pages we grab columns from, and the columns that we grab
        ttools = ["km","Longitude","Latitude"]
        morph = ["Elevation","S","W_b","z","n","SedThermCond","SedThermDiff","SedDepth",
                 "hyp_percent","phi","Q_cont","d_cont"]
        flow = ["Q_in","T_in","Q_out"]
        # Ways that we grab the columns
        sums = ["hyp_percent","Q_in","Q_out"] # These are summed, not averaged
        mins = ["km"]
        aves = ["Longitude","Latitude","Elevation","S","W_b","z","n","SedThermCond",
                "SedThermDiff","SedDepth","phi", "Q_cont","d_cont","T_in"]
        ignore = ["FLIR_Temp","FLIR_Time"] # Ignore in the loop, set them manually

        data = {}
        
        # Get all the "columnar data" from the csv files. Note Stream Km is brought in as an index file and set as data type object.
        # This is not a perfect method but it is done to avoid trailing on floating point values messing with the dictionary.
        # TODO organize the reading of this data into seperate methods so it is easier to find later.
        lcdata = pd.read_csv(join(IniParams["inputdir"],IniParams["lcdatafile"]),quotechar='"',quoting=0)
        morphdata = pd.read_csv(join(IniParams["inputdir"],IniParams["morphfile"]),quotechar='"',quoting=0,index_col='Stream_KM')
        accdata = pd.read_csv(join(IniParams["inputdir"],IniParams["accretionfile"]),quotechar='"',quoting=0,index_col='Stream_KM')
        
        # Add these columns to morph data since they do not exist in the input file.
        morphdata["Q_cont"] = 0
        morphdata["d_cont"] = 0
        
        # Check for missing values # TODO Fix this
        #if pd.isnull(morphdata):
            #raise Exception("Missing data in %s" % IniParams["morphfile"])
        #if pd.isnull(lcdata):
            #raise Exception("Missing data in %s" % IniParams["lcdatafile"])
        
        # With accretion data replace NaN with zeros
        accdata = accdata.fillna(0)
        print("Missing data in %s converted to zeros" % IniParams["accretionfile"])
            
        for i in xrange(len(ttools)):
            data[ttools[i]] = list(lcdata.ix[:,i])
        for i in xrange(len(morph)):
            data[morph[i]] = list(morphdata.ix[:,i])
        for i in xrange(len(flow)):
            data[flow[i]] = list(accdata.ix[:,i])
                    
        #Longitude check
        if max(data[ttools[1]]) > 180:
            long_list = list(data[ttools[1]])
            max1 = max(long_list)
            index1 = long_list.index(max1)
            rkm = data[ttools[0]][index1]
            raise Exception("Longitude must be less than 180 degrees. At rKM = " + str(rkm) + ", longitude = " + str(max1))
        if min(data[ttools[1]]) < -180:
            long_list = list(data[ttools[1]])
            min1 = min(long_list)
            index1 = long_list.index(min1)
            rkm = data[ttools[1]][index1]
            raise Exception("Longitude must be greater than -180 degrees. At rKM = " + str(rkm) + ", longitude = " + str(min1))
        #Latitude check
        if max(data[ttools[2]]) > 90 or min(data[ttools[2]]) < -90:
            raise Exception("Latitude must be greater than -90 and less than 90 degrees")

        # Then sum and average things as appropriate. multiplier() takes a tuple
        # and applies the given lambda function to that tuple.
        for attr in sums:
            data[attr] = self.multiplier(data[attr],lambda x:sum(x))
        for attr in aves:
            data[attr] = self.multiplier(data[attr],lambda x:sum(x)/len(x))
        for attr in mins:
            data[attr] = self.multiplier(data[attr],lambda x:min(x))
        return data

    def BuildNodes(self):
        # This is the worst of the methods but it works. # TODO
        self.CheckEarlyQuit()
        print("Building Stream Nodes")
        Q_mb = 0.0
        # Grab all of the data in a dictionary
        data = self.GetColumnarData()
        #################################
        # Build a boundary node
        node = StreamNode(run_type=self.run_type,Q_mb=Q_mb)
        # Then set the attributes for everything in the dictionary
        for k,v in data.iteritems():
            setattr(node,k,v[0])
        # set the flow and temp boundary conditions for the boundary node
        node.Q_bc = self.Q_bc
        node.T_bc = self.T_bc
        self.InitializeNode(node)
        node.dx = IniParams["longsample"]
        self.Reach[node.km] = node
        ############################################

        #Figure out how many nodes we should have downstream. We use math.ceil() because
        # if we end up with a fraction, that means that there's a node at the end that
        # is not a perfect multiple of the sample distance. We might end up ending at
        # stream kilometer 0.5, for instance, in that case
        vars = (IniParams["length"] * 1000)/IniParams["longsample"]

        num_nodes = int(ceil((vars)/self.multiple))
        for i in range(0, num_nodes):
            node = StreamNode(run_type=self.run_type,Q_mb=Q_mb)
            for k,v in data.iteritems():
                setattr(node,k,v[i+1])# Add one to ignore boundary node
            self.InitializeNode(node)
            self.Reach[node.km] = node
            print("Building Stream Nodes", i+1, num_nodes)
        # Find the mouth node and calculate the actual distance
        mouth = self.Reach[min(self.Reach.keys())]
        mouth_dx = (vars)%self.multiple or 1.0 # number of extra variables if we're not perfectly divisible
        mouth.dx = IniParams["longsample"] * mouth_dx

    def BuildZonesNormal(self):
        """This method builds the sampled vegzones in the case of non-lidar datasets"""
        # Hide your straight razors. This implementation will make you want to use them on your wrists.
        self.CheckEarlyQuit()
        LCcodes = self.GetLandCoverCodes() # Pull the LULC codes
        LCdata = self.GetLandCoverData() # Pull the LULC Data
        
        vheight = []
        vdensity = []
        k = []
        overhang = []
        elevation = []
        average = lambda x:sum(x)/len(x)
        trans_count = IniParams["transsample_count"]
        if IniParams["heatsource8"] == True:
            radial_count = 7 # heat source 8 default
        else:
            radial_count = IniParams["radialsample_count"]

        keys = self.Reach.keys()
        keys.sort(reverse=True) # Downstream sorted list of stream kilometers
        print("Translating LULC Data")
        for i in xrange(6, radial_count*trans_count+7): # For each column of LULC data
            col = list(LCdata.ix[:,i]) # LULC column
            elev = list(LCdata.ix[:,i+radial_count*trans_count]) # Shift by 28 to get elevation column
            # Make a list from the LC codes from the column, then send that to the multiplier
            # with a lambda function that averages them appropriately. Note, we're averaging over
            # the values (e.g. density) not the actual code, which would be meaningless.
            try:
                vheight.append(self.multiplier([LCcodes[x][0] for x in col], average))
                vdensity.append(self.multiplier([LCcodes[x][1] for x in col], average))
                k.append(self.multiplier([LCcodes[x][2] for x in col], average))
                overhang.append(self.multiplier([LCcodes[x][3] for x in col], average))
            except KeyError, (stderr):
                raise Exception("At least one land cover code in %s is blank or not in %s (Code: %s)." % (IniParams["lcdatafile"], IniParams["lccodefile"], stderr.message))
            if i>6:  # There isn't a stream center elevation (that is in the morphology file), so we don't want to read in first elevation value which s actually the last LULC col.
                elevation.append(self.multiplier(elev, average))
            print("Translating LULC Data", i, radial_count*trans_count+7)
            
        for i in xrange(len(keys)):
            node = self.Reach[keys[i]]
            xx = 0
            for dir in xrange(radial_count+1):
                for zone in xrange(trans_count):
                    node.VHeight[dir][zone] = vheight[xx][i]
                    node.VDensity[dir][zone] = vdensity[xx][i]
                    node.k[dir][zone] = k[xx][i]
                    node.Overhang[dir][zone] = overhang[xx][i]
                    xx = xx + 1
                    if dir == 0 and zone == 0: # 0 is emergent, there is only one value at zone = 0
                        break # go to the next dir
                    
        # Average over the topo values
        topo_w = self.multiplier(list(LCdata.TopoWest.values), average)
        topo_s = self.multiplier(list(LCdata.TopoSouth.values), average)
        topo_e = self.multiplier(list(LCdata.TopoEast.values), average)

        # ... and you thought things were crazy earlier! Here is where we build up the
        # values for each node. This is culled from earlier version's VB code and discussions
        # to try to simplify it... yeah, you read that right, simplify it... you should've seen in earlier!
        
        for h in xrange(len(keys)):
            print("Building VegZones", h+1, len(keys))
            node = self.Reach[keys[h]]
            VTS_Total = 0 #View to sky value
            LC_Angle_Max = 0
            # Now we set the topographic elevations in each direction
            node.TopoFactor = (topo_w[h] + topo_s[h] + topo_e[h])/(90*3) # Topography factor Above Stream Surface
            # This is basically a list of directions, each with one of three topographies
            ElevationList = []
            Angle_Incr = 360.0 / radial_count
            DirNumbers = range(1,radial_count+1)
            AngleMid = [x*Angle_Incr for x in DirNumbers]
            for i in xrange(radial_count): # Iterate through each direction
                DirAngle = AngleMid[i]
                if DirAngle < 135:
                    ElevationList.append(topo_e[h])
                elif DirAngle < 225:
                    ElevationList.append(topo_s[h])
                else:
                    ElevationList.append(topo_w[h])
            # Sun comes down and can be full-on, blocked by veg, or blocked by topography. Earlier implementations
            # calculated each case on the fly. Here we chose a somewhat more elegant solution and calculate necessary
            # angles. Basically, there is a minimum angle for which full sun is calculated (top of trees), and the
            # maximum angle at which full shade is calculated (top of topography). Anything in between these is an
            # angle for which sunlight is passing through trees. So, for each direction, we want to calculate these
            # two angles so that late we can test whether we are between them, and only do the shading calculations
            # if that is true.

            for i in xrange(radial_count): # Iterate through each direction
                T_Full = () # lowest angle necessary for full sun
                T_None = () # Highest angle necessary for full shade
                W_Vdens_num = 0.0  #Numerator for the weighted Veg density calculation
                W_Vdens_dem = 0.0  #Denominator for the weighted Veg density calculation
                
                for j in xrange(trans_count): # Iterate through each of the zones
                    Vheight = vheight[i*trans_count+j+1][h]
                    Vdens = vdensity[i*trans_count+j+1][h]
                    Overhang = overhang[i*trans_count+j+1][h]
                    Elev = elevation[i*trans_count+j][h]

                    if not j: # We are at the stream edge, so start over
                        LC_Angle_Max = 0 # New value for each direction
                    else:
                        Overhang = 0 # No overhang away from the stream
                    ##########################################################
                    # Calculate the relative ground elevation. This is the
                    # vertical distance from the stream surface to the land surface
                    SH = Elev - node.Elevation
                    # Then calculate the relative vegetation height
                    VH = Vheight + SH

                    # Calculate the node distance
                    #Adjustment for whether the veg sample represent a zone (see excel interface for explanation)
                    if IniParams["vegDistMethod"] == "zone":
                        adjust = 0.5
                    else:
                        adjust = 0.0
                    
                    # TODO this is a future function where there is a landcover sample at the stream node   
                    if IniParams["heatsource8"] == True:
                        adjust2 = 1
                    else:  # if adjust2 = 0 there is a landcover sample at the stream node 
                        adjust2 = 1
                    
                    LC_Distance = IniParams["transsample"] * (j + adjust2 - adjust)
                    # We shift closer to the stream by the amount of overhang
                    # This is a rather ugly cludge.
                    if not j: LC_Distance -= Overhang
                    if LC_Distance <= 0:
                        LC_Distance = 0.00001
                    # Calculate the minimum sun angle needed for full sun
                    T_Full += degrees(atan(VH/LC_Distance)), # It gets added to a tuple of full sun values
                    
                    # Now get the maximum of bank shade and topographic shade for this direction
                    T_None += degrees(atan(SH/LC_Distance)), # likewise, a tuple of values
                    
                    # Calculate View To Sky
                    veg_angle = degrees(atan(VH/LC_Distance)) - degrees(atan(SH/LC_Distance))
                    if IniParams["beers_data"] == "LAI": #use LAI data
                        LAI_den = Vdens / 12 #  Purpose here is to deveop a density. A LAI of 12 is pretty much closed canopy
                        if LAI_den > 1:
                            LAI_den = 1
                        W_Vdens_num += veg_angle*float(LAI_den)
                    else:
                        W_Vdens_num += veg_angle*float(Vdens)
                    W_Vdens_dem += veg_angle
                    
                    if j == trans_count - 1:
                        if max(T_Full) > 0:   #if bank and/or veg shade is occuring:
                            #Find weighted average the density:
                            #Vdens_mod = (Amount of Veg shade * Veg dens) + (Amount of bank shade * bank dens, i.e. 1) / (Sum of amount of shade)
                            if W_Vdens_dem > 0:
                                Vdens_ave_veg = W_Vdens_num / W_Vdens_dem
                            else:
                                Vdens_ave_veg = 0
                            Vdens_mod = ((max(T_Full)-max(T_None))* Vdens_ave_veg + max(T_None)) / max(T_Full)
                        else:
                            Vdens_mod = 1.0
                        VTS_Total += max(T_Full)*Vdens_mod # Add angle at end of each zone calculation
                node.ShaderList += (max(T_Full), ElevationList[i], max(T_None), T_Full),
            node.ViewToSky = 1 - VTS_Total / (radial_count * 90)


    def BuildZonesLidar(self):
        """Build zones if we are using LiDAR data"""
        #self.CheckEarlyQuit()
        #Tried to keep in the same general form as BuildZonesNormal so blame Metta
        self.CheckEarlyQuit()
        LCdata = self.GetLandCoverData() # Pull the LULC Data
               
        vheight = []
        vdens = []
        k = []
        elevation = []
        average = lambda x:sum(x)/len(x)
        trans_count = IniParams["transsample_count"]
        if IniParams["heatsource8"] == True:
            radial_count = 7 # heat source 8 default
        else:
            radial_count = IniParams["radialsample_count"]

        keys = self.Reach.keys()
        keys.sort(reverse=True) # Downstream sorted list of stream kilometers
        print("Translating LULC Data")
        for i in xrange(6, radial_count*trans_count+7): # For each column of LULC data          
            col = list(LCdata.ix[:,i]) # LULC column
            elev = list(LCdata.ix[:,i+radial_count*trans_count]) # Shift by 7 * "number of trans sample zones" to get elevation column
            if IniParams["heatsource8"] == True:
                dens = list(LCdata.ix[:,i+1+radial_count*trans_count*2])
            else:
                dens = [IniParams["lcdensity"]]*len(col)
            # Make a list from the LC codes from the column, then send that to the multiplier
            # with a lambda function that averages them appropriately. Note, we're averaging over
            # the values (e.g. density) not the actual code, which would be meaningless.
            try:
                vheight.append(self.multiplier([x for x in col], average))
                vdens.append(self.multiplier([x for x in dens], average))
                #k.append(self.multiplier([x for x in col], average)) #TODO
            except KeyError, (stderr):
                raise Exception("Vegetation height/density error" % stderr.message)
            if i>6:  # There isn't a stream center elevation (that is in the morphology file), so we don't want to read in first elevation value which s actually the last LULC col.
                elevation.append(self.multiplier(elev, average))
            print("Reading vegetation heights", i+1, radial_count*trans_count+7)
            
        for i in xrange(len(keys)):
            node = self.Reach[keys[i]]
            xx = 0
            for dir in xrange(radial_count+1):
                for zone in xrange(trans_count):
                    node.VHeight[dir][zone] = vheight[xx][i]
                    node.VDensity[dir][zone] = vdensity[xx][i]
                    node.k[dir][zone] = k[xx][i]
                    node.Overhang[dir][zone] = overhang[xx][i]
                    xx = xx + 1
                    if dir == 0 and zone == 0: # 0 is emergent, there is only one value at zone = 0
                        break # go to the next dir
            
        # Average over the topo values
        topo_w = self.multiplier(list(LCdata.TopoWest.values), average)
        topo_s = self.multiplier(list(LCdata.TopoSouth.values), average)
        topo_e = self.multiplier(list(LCdata.TopoEast.values), average)

        # ... and you thought things were crazy earlier! Here is where we build up the
        # values for each node. This is culled from earlier version's VB code and discussions
        # to try to simplify it... yeah, you read that right, simplify it... you should've seen in earlier!
        
        for h in xrange(len(keys)):
            print("Building VegZones", h+1, len(keys))
            node = self.Reach[keys[h]]
            LC_Angle_Max = 0
            VTS_Total = 0 #View to sky value
            # Now we set the topographic elevations in each direction
            node.TopoFactor = (topo_w[h] + topo_s[h] + topo_e[h])/(90*3) # Topography factor Above Stream Surface
            # This is basically a list of directions, each with one of three topographies
            ElevationList = []
            Angle_Incr = 360.0 / radial_count
            DirNumbers = range(1,radial_count+1)
            AngleMid = [x*Angle_Incr for x in DirNumbers]
            for i in xrange(radial_count): # Iterate through each direction
                DirAngle = AngleMid[i]
                if DirAngle < 135:
                    ElevationList.append(topo_e[h])
                elif DirAngle < 225:
                    ElevationList.append(topo_s[h])
                else:
                    ElevationList.append(topo_w[h])
            # Sun comes down and can be full-on, blocked by veg, or blocked by topography. Earlier implementations
            # calculated each case on the fly. Here we chose a somewhat more elegant solution and calculate necessary
            # angles. Basically, there is a minimum angle for which full sun is calculated (top of trees), and the
            # maximum angle at which full shade is calculated (top of topography). Anything in between these is an
            # angle for which sunlight is passing through trees. So, for each direction, we want to calculate these
            # two angles so that late we can test whether we are between them, and only do the shading calculations
            # if that is true.

            for i in xrange(7): # Iterate through each direction
                T_Full = () # lowest angle necessary for full sun
                T_None = () # Highest angle necessary for full shade
                W_Vdens_num = 0.0  #Numerator for the weighted Veg density calculation
                W_Vdens_dem = 0.0  #Denominator for the weighted Veg density calculation

                for j in xrange(trans_count): # Iterate through each of the zones
                    Vheight = vheight[i*trans_count+j+1][h]
                    if Vheight < 0 or Vheight is None or Vheight > 120:
                        raise Exception("Vegetation height (value of %s in TTools Data) must be greater than zero and less than 120 meters (when LiDAR = True)" % `Vheight`)
                    Vdens = vdens[i*trans_count+j+1][h]
                    Overhang = IniParams["lcoverhang"]
                    Elev = elevation[i*trans_count+j][h]

                    if not j: # We are at the stream edge, so start over
                        LC_Angle_Max = 0 # New value for each direction
                    else:
                        Overhang = 0 # No overhang away from the stream
                    ##########################################################
                    # Calculate the relative ground elevation. This is the
                    # vertical distance from the stream surface to the land surface
                    SH = Elev - node.Elevation
                    # Then calculate the relative vegetation height
                    VH = Vheight + SH

                    # Calculate the node distance.
                    #Different for LiDAR because we assume you are sampling a tree at a specific location
                    #rather than a veg zone which represents the vegetation between two sample points
                    #Adjustment for whether the veg sample represent a zone (see excel interface for explanation)
                    if IniParams["vegDistMethod"] == "zone":
                        adjust = 0.5
                    else:
                        adjust = 0.0
                    LC_Distance = IniParams["transsample"] * (j + 1 - adjust) #This is "+ 1" because j starts at 0
                    # We shift closer to the stream by the amount of overhang
                    # This is a rather ugly cludge.
                    if not j: LC_Distance -= Overhang
                    if LC_Distance <= 0:
                        LC_Distance = 0.00001
                    # Calculate the minimum sun angle needed for full sun
                    T_Full += degrees(atan(VH/LC_Distance)), # It gets added to a tuple of full sun values
                    # Now get the maximum of bank shade and topographic shade for this
                    # direction
                    T_None += degrees(atan(SH/LC_Distance)), # likewise, a tuple of values
                    
                    # Calculate View To Sky
                    veg_angle = degrees(atan(VH/LC_Distance)) - degrees(atan(SH/LC_Distance))
                    if IniParams["beers_data"] == "LAI": #use LAI data
                        LAI_den = Vdens / 12 #  Purpose here is to deveop a density. A LAI of 12 is pretty much closed canopy
                        if LAI_den > 1:
                            LAI_den = 1
                        W_Vdens_num += veg_angle*float(LAI_den)
                    else:
                        W_Vdens_num += veg_angle*float(Vdens)
                    W_Vdens_dem += veg_angle
                    
                    if j == trans_count - 1:
                        if max(T_Full) > 0:   #if bank and/or veg shade is occuring:
                            #Find weighted average the density:
                            #Vdens_mod = (Amount of Veg shade * Veg dens) + (Amount of bank shade * bank dens, i.e. 1) / (Sum of amount of shade)
                            if W_Vdens_dem > 0:
                                Vdens_ave_veg = W_Vdens_num / W_Vdens_dem
                            else:
                                Vdens_ave_veg = 0
                            Vdens_mod = ((max(T_Full)-max(T_None))* Vdens_ave_veg + max(T_None)) / max(T_Full)
                        else:
                            Vdens_mod = 1.0
                        VTS_Total += max(T_Full)*Vdens_mod # Add angle at end of each zone calculation
                node.ShaderList += (max(T_Full), ElevationList[i], max(T_None), T_Full),
            node.ViewToSky = 1 - VTS_Total / (radial_count * 90)

    def GetLandCoverCodes(self):
        """Return the codes from the Land Cover Codes csv input file as a dictionary of dictionaries"""
        self.CheckEarlyQuit()       
        codes,height,dens,k,over = [],[],[],[],[]                
        
        data = pd.read_csv(join(IniParams["inputdir"],IniParams["lccodefile"]),quotechar='"')
        data.columns = ['name', 'code','height','density','k','overhang']
        
        codes = list(data.code.values)
        # make a list of lists with values: [(height[0], dens[0], over[0]), (height[1],...),...]
        vals = [tuple([j for j in i]) for i in zip(data.height.values,data.density.values,data.k.values,data.overhang.values)]
        data = {}

        for i in xrange(len(codes)):
            # Each code is a tuple in the form of (VHeight, VDensity, Overhang)
            data[codes[i]] = vals[i]
            
            if IniParams["beers_data"] == "LAI": #use LAI data
                if vals[i][0] != None and (vals[i][1] < 0):
                    raise Exception("Vegetation Density (value of %s in Land Cover Codes) must be >= 0.0 and <= 1.0" % `vals[i][1]`)           
            else:
                if vals[i][0] != None and (vals[i][1] < 0 or vals[i][1] > 1):
                    raise Exception("Vegetation Density (value of %s in Land Cover Codes) must be >= 0.0 and <= 1.0" % `vals[i][1]`)
        return data
        
    def GetLandCoverData(self):
        """Return all data from the Land Cover Data csv input file as a pandas dataframe"""
        self.CheckEarlyQuit()
        data = pd.read_csv(join(IniParams["inputdir"],IniParams["lcdatafile"]),quotechar='"',quoting=0)
        return data

    def InitializeNode(self, node):
        """Perform some initialization of the StreamNode, and write some values to spreadsheet"""
        # Initialize each nodes tribs dictionary to a tuple
        for time in self.flowtimelist:
            node.Q_tribs[time] = ()
            node.T_tribs[time] = ()
        ##############################################################
        #Now that we have a stream node, we set the node's dx value, because
        # we have most nodes that are long-sample-distance times multiple,
        node.dx = IniParams["dx"] # Nodes distance step.
        node.dt = IniParams["dt"] # Set the node's timestep... this may have to be adjusted to comply with stability
        # Find the earliest temperature boundary condition
        mindate = min(self.T_bc.keys())
        if self.run_type == 2: # Running hydraulics only
            node.T, node.T_prev, node.T_sed = 0.0, 0.0, 0.0
        else:
            if self.T_bc[mindate] is None:
                # Shade-a-lator doesn't need a boundary condition
                if self.run_type == 1: self.T_bc[mindate] = 0.0
                else:  raise Exception("Boundary temperature conditions cannot be blank")
            node.T = self.T_bc[mindate]
            node.T_prev = self.T_bc[mindate]
            node.T_sed = self.T_bc[mindate]
        #we're in shadealator if the runtype is 1. Since much of the heat
        # math is coupled to the shade math, we have to make sure the hydraulic
        # values are not zero or blank because they'll raise ZeroDivisionError
        if self.run_type ==1:
            for attr in ["d_w", "A", "P_w", "W_w", "U", "Disp","Q_prev","Q",
                         "SedThermDiff","SedDepth","SedThermCond"]:
                if (getattr(node, attr) is None) or (getattr(node, attr) == 0):
                    setattr(node, attr, 0.01)
        node.Q_hyp = 0.0 # Assume zero hyporheic flow unless otherwise calculated
        node.E = 0 # Same for evaporation
    
    def QuitMessage(self):
        b = buttonbox("Do you really want to quit Heat Source", "Quit Heat Source", ["Cancel", "Quit"])
        if b == "Quit":
            raise Exception("Model stopped user.")
        else: return

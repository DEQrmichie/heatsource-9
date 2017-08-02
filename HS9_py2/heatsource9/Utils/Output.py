# Heat Source, Copyright (C) 2000-2016, 
# Oregon Department of Environmental Quality

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

from __future__ import division, print_function
from time import ctime
from os.path import join
from copy import deepcopy
from math import ceil
import csv

from ..Dieties.IniParamsDiety import IniParams

class Output(object):
    """Data and file object storage class"""
    def __init__(self, reach, start_time, run_type):
        # Store a sorted list of StreamNodes.
        
        if IniParams["outputkm"].lower() == "all":
            self.nodes = sorted(reach.itervalues(),reverse=True)
        
        else:
            # specfic outputkm
            self.nodes = sorted([reach[km] for km in reach if km in IniParams["outputkm"]],
                                reverse= True)       
        
        # A reference to the model's starting time 
        # (i.e. when spin-up is over)
        self.start_time = start_time

        # run_type is a bit hack-y. If we are running only hydraulics,
        # we fail on division of solar parameters- if running only solar,
        # we fail on hydraulics. This is an easy way to prevent that.
        self.run_type = run_type #0=HS, 1=solar, 2=hydraulics
        # Our first time through, we ignore daily data and don't have 
        # stream geometry calculated, so we have switches for those 
        # (which is a bit dumb)
        self.first_hour = True
        self.first_day = True

        # Filenames and descriptions for each of the output files
        desc = {}
        
        if (run_type in [0, 1] and False):
            # The control file will hold a switch to turn these on
            desc["Heat_DR1"] = "Direct Solar Radiation Flux above Topographic Features (watts/square meter)"
            desc["Heat_DR2"] = "Direct Solar Radiation Flux below Topographic Features (watts/square meter)"
            desc["Heat_DR4"] = "Direct Solar Radiation Flux above the Stream (watts/square meter)"
            desc["Heat_DR5"] = "Direct Solar Radiation Flux Entering Stream (watts/square meter)"
            desc["Heat_DF1"] = "Diffuse Solar Radiation Flux above Topographic Features (watts/square meter)"
            desc["Heat_DF2"] = "Diffuse Solar Radiation Flux below Topographic Features (watts/square meter)"
            desc["Heat_DF4"] = "Diffuse Solar Radiation Flux above the Stream (watts/square meter)"
            desc["Heat_DF5"] = "Diffuse Solar Radiation Flux Entering Stream (watts/square meter)"
            desc["Solar_Azimuth"] = "Solar Azimuth (degrees). Angle of the sun along the horizon with zero degrees corresponding to true north."
            desc["Solar_Elevation"] = "Solar Elevation (degrees). Angle up from the horizon to the sun."            
        
        if run_type in [0, 1]:
            # Solar, Temperature
            desc["Heat_SR1"] = "Solar Radiation Flux above Topographic Features (watts/square meter)"
            desc["Heat_SR2"] = "Solar Radiation Flux below Topographic Features (watts/square meter)"
            desc["Heat_SR3"] = "Solar Radiation Flux below Land Cover (watts/square meter)"
            desc["Heat_SR3b"] = "Solar Radiation Flux blocked by Land Cover (watts/square meter)"
            desc["Heat_SR4"] = "Solar Radiation Flux above the Stream (watts/square meter)"
            desc["Heat_SR5"] = "Solar Radiation Flux Entering Stream (watts/square meter)"
            desc["Shade"] = "Effective Shade"
            desc["VTS"] = "View to Sky. Calculated proportion of the sky hemisphere obscured by land cover."

        if run_type in [0, 2]:
            # Hydro, Temperature
            desc["Hyd_DA"] = "Average Depth (meters)"
            desc["Hyd_DM"] = "Max Depth (meters)"
            desc["Hyd_Flow"] = "Flow Rate (cms)"
            desc["Hyd_Hyp"] = "Hyporheic Exchange (cms)"
            desc["Hyd_Vel"] = "Flow Velocity (meters/second)"
            desc["Hyd_WT"] = "Top Width (meters)"
        if run_type == 0:
            # Temperature
            desc["Heat_SR6"] = "Solar Radiation Flux Received by Stream (watts/square meter)"
            desc["Heat_SR7"] = "Solar Radiation Flux Received by Substrate (watts/square meter)"
            desc["Heat_Cond"] = "Streambed Conduction Flux (watts/square meter)"
            desc["Heat_Long"] = "Longwave Flux (watts/square meter)"
            desc["Heat_Conv"] = "Convection Flux (watts/square meter)"
            desc["Heat_Evap"] = "Evaporation Flux (watts/square meter)"
            desc["Rate_Evap"] = "Evaporation Rate (mm/hour)"
            desc["Temp_H2O"] = "Stream Temperature (*C)"
            desc["Temp_Sed"] = "Sediment Temperature (*C)"
            desc["Hyd_Disp"] = "Hydraulic Dispersion (square meters/second)"

        # Storage dictionary for the data.
        self.data = {}
        for name in desc.keys():
            self.data[name] = {}
        # make a deepcopy of the empty variables dictionary for use later
        self.empty_vars = deepcopy(self.data)
        # Empty dictionary to store file objects
        self.files = {}

        # Here we build up the self.files attribute by cycling through 
        # the filenames and descriptions
        for key in desc.iterkeys():
            # Build the header that will be stamped to each output file
            header = [["File Created:"] + [ctime()]]
            header += [["Heat Source Version:"] + [IniParams["version"]]]
            header += [["Simulation Name:"] + [IniParams["name"]]]
            header += [["User Text:"] + [IniParams["usertxt"]]]
            header += [["Output:"] + [desc[key]]]
            header += [[""]]
            header += [["Datetime"]]
            if key in ["Heat_SR3b"] :
                
                header[6] += ["STREAM_KM"]
                
                if IniParams["heatsource8"]:
                    # Same as 8 directions but no north
                    dir = ["NE","E","SE","S","SW","W","NW"]
                    zone = range(1,int(IniParams["transsample_count"])+1)
                else:
                    dir = ["T" + str(x) for x in range(1,IniParams["trans_count"]+ 1)]
                    zone = range(1,int(IniParams["transsample_count"])+1)
                    
                    # TODO this is a future fuction to have a landcover 
                    # sample at the streamnode
                    #zone = range(0,int(IniParams["transsample_count"]))
                if IniParams["lcdatainput"] == "Values":
                    type = "HT"
                else:
                    type = "LC"
                
                for d in range(0,len(dir)):
                    for z in range(0,len(zone)):
                        header[6] += [type+"_"+ dir[d]+"_S"+str(zone[z])]
                
                header[6] += ["Diffuse_Blocked"]
            else: 
                # Grab a list of all the stream kilometers
                header[6] += [("%0.3f" % x.km) for x in self.nodes]
                           
            # Now create a file object in the dictionary, and write 
            # the header
            self.files[key] = csv.writer(open(join(IniParams["outputdir"], key + ".csv"), "wb"))
            self.files[key].writerows(header)
        
    def __call__(self, time, hour, minute=0, second=0):
        """Call the storage method with a time and an hour"""
        # Ignore this if we're still spinning up of if this is the first
        # hour run (because we don't have channel geometry calculated).
        if time < self.start_time: return
        if self.first_hour:
            self.first_hour = False
            #return
        # Create an Excel-friendly time string
        timestamp = ("%0.7f" % float(time/86400 + 25569))
        # Localize variables to save a bit of time
        nodes = self.nodes
        data = self.data
        # Cycle through each datatype, creating a list of values
        # corresponding to the nodes for this timestamp. Thus, each
        # timestamp conforms to a single line, and each element in the
        # list comprehension conforms to a column. List comprehensions
        # are generally fast (more optimized by the underlying C code)
        # than for loops.

        if (self.run_type < 2 and False):
            # The control file will hold a switch to turn these on.
            data["Solar_Elevation"][timestamp] = [x.head.SolarPos[0] for x in nodes]
            data["Solar_Azimuth"][timestamp] = [x.head.SolarPos[4] for x in nodes]
            data["Heat_DF1"][timestamp] = [x.F_Diffuse[1] for x in nodes]
            data["Heat_DF2"][timestamp] = [x.F_Diffuse[2] for x in nodes]
            data["Heat_DF3"][timestamp] = [x.F_Diffuse[3] for x in nodes]
            data["Heat_DF4"][timestamp] = [x.F_Diffuse[4] for x in nodes]
            data["Heat_DF5"][timestamp] = [x.F_Diffuse[5] for x in nodes]            
            data["Heat_DR1"][timestamp] = [x.F_Direct[1] for x in nodes]
            data["Heat_DR2"][timestamp] = [x.F_Direct[2] for x in nodes]
            data["Heat_DR3"][timestamp] = [x.F_Direct[3] for x in nodes]
            data["Heat_DR4"][timestamp] = [x.F_Direct[4] for x in nodes]
            data["Heat_DR5"][timestamp] = [x.F_Direct[5] for x in nodes]
            
        # Run only with solar
        if self.run_type < 2:         
            data["Heat_SR1"][timestamp] = [x.F_Solar[1] for x in nodes]
            data["Heat_SR2"][timestamp] = [x.F_Solar[2] for x in nodes]
            data["Heat_SR3"][timestamp] = [x.F_Solar[3] for x in nodes]
            data["Heat_SR4"][timestamp] = [x.F_Solar[4] for x in nodes]
            data["Heat_SR5"][timestamp] = [x.F_Solar[5] for x in nodes]
            
        # Run only with hydro
        if self.run_type != 1:
            data["Hyd_DA"][timestamp] = [(x.A / x.W_w) for x in nodes]
            data["Hyd_DM"][timestamp] = [x.d_w for x in nodes]
            data["Hyd_Flow"][timestamp] = [x.Q for x in nodes]
            data["Hyd_Hyp"][timestamp] = [x.Q_hyp for x in nodes]
            data["Hyd_Vel"][timestamp] = [x.U for x in nodes]
            data["Hyd_WT"][timestamp] = [x.W_w for x in nodes]
            
        # Run only with both solar and hydro
        if not self.run_type:
            data["Heat_SR6"][timestamp] = [x.F_Solar[6] for x in nodes]
            data["Heat_SR7"][timestamp] = [x.F_Solar[7] for x in nodes]
            data["Heat_Cond"][timestamp] = [x.F_Conduction for x in nodes]
            data["Heat_Conv"][timestamp] = [x.F_Convection for x in nodes]
            data["Heat_Evap"][timestamp] = [x.F_Evaporation for x in nodes]
            data["Heat_Long"][timestamp] = [x.F_Longwave for x in nodes]
            data["Rate_Evap"][timestamp] = [(x.E / x.dx / x.W_w * 3600 * 1000) for x in nodes] #TODO: Check
            data["Temp_H2O"][timestamp] = [x.T for x in nodes]
            data["Temp_Sed"][timestamp] = [x.T_sed for x in nodes]
            data["Hyd_Disp"][timestamp] = [x.Disp for x in nodes]

        # Zero for an hour means a new day, so we add daily outputs
        # and write to the file. Writing only every day saves us
        # 24xF file accesses where F=len(self.files). Each file access
        # has quite a bit of overhead, so we lump them.
        # Run the daily output on the last hour of the day
        if hour == 23:
            self.write_to_csv(self.run_type < 2, timestamp)
            
        # Run the daily outputs if we are one timestep to
        # the end of the current day. Uncomment when there is an 
        # output every timestep. Run times can get long.
        #if (((23 - hour) * 3600) + ((60-minute) * 60) <=  IniParams["dt"]):
            #self.write_to_csv(self.run_type < 2, timestamp)

    def daily(self, timestamp):
        """Compile and store data that is collected every hour"""
        nodes = self.nodes
        self.data["Shade"][timestamp] = [((x.F_DailySum[1] - x.F_DailySum[4]) / x.F_DailySum[1]) for x in nodes]
        self.data["VTS"][timestamp] = [x.ViewToSky for x in nodes]
        self.data["Heat_SR3b"][timestamp] = []
        # If there's no hour, we're at the beginning of a day, so we 
        # write the values to a file.

    def write_to_csv(self, daily, timestamp):
        if daily:
            # don't call for hydraulics
            self.daily(timestamp)
        # localize the data
        data = self.data
        # Cycle through the file objects
        for name, fileobj in self.files.iteritems():
            # Each time is a single line, so we want to iterate over all 
            # the times stored so far. We can do this because everytime 
            # we store data, we append the time string to self.times
            timelist = sorted(data[name].keys())
            line = []
            if name in ["Heat_SR3b"]:
                timesteps = 86400.0/float(IniParams["dt"])
                for timestamp in timelist:
                    for i, x in enumerate(self.nodes):
                        line += [[timestamp] + [("%0.3f" % x.km)]]
                        if IniParams["heatsource8"] == True:
                            # Seven directions
                            directions = [d for d in range(0,7)]
                        else:
                            directions = [d for d in range(IniParams["trans_count"])]
                        
                        # solar blocked
                        for d in directions:
                            for zone in range(IniParams["transsample_count"]):
                                daily_avg = x.Solar_Blocked[d][zone]  / timesteps
                                line[i] += ["%0.4f" % daily_avg]
                        daily_avg_diffuse = x.Solar_Blocked["diffuse"]  / timesteps
                        line[i] += ["%0.4f" % daily_avg_diffuse]

            else:
                for timestamp in timelist:
                    line += [[timestamp] + [("%0.4f" % x) for x in data[name][timestamp]]]
            
            # finally, write all the lines to the file
            fileobj.writerows(line)
        del data
        # Now empty out the dictionary by simply copying a new one.
        self.data = deepcopy(self.empty_vars)
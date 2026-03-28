"""This module collects model output values and writes them to CSV files."""


from copy import deepcopy
from os.path import join
from time import ctime
import csv

from heatsource9.__version__ import __version__


class OutputWriter:
    """Class to manage output writing during a model run."""

    def __init__(self):
        """Initialize an empty output."""
        self._output= None

    def bind(self, reach, params, run_type):
        """Stores the reach, model parameters, and run type needed for output files"""
        start_time = params["modelstart"]
        stop_time = params["modelend"]
        self._output = _Output(reach=reach, start_time=start_time, stop_time=stop_time, run_type=run_type, params=params)

    def write_step(self, simulation, time, hour, minute, second):
        """Write output values for one model timestep."""
        self._output(time=float(time), hour=int(hour), minute=int(minute), second=int(second))

    def close(self):
        if self._output is not None:
            self._output.close()
            self._output = None


class _Output(object):
    """Class to store output data in memory and write output files."""

    def __init__(self, reach, start_time, stop_time, run_type, params):
        self.params = params

        # Store a sorted list of StreamNodes.
        if str(self.params["outputkm"]).lower() == "all":
            self.nodes = sorted(iter(list(reach.values())), reverse=True)
        else:
            # specific outputkm (already converted to reach keys by ModelSetup)
            self.nodes = sorted([reach[km] for km in reach if km in self.params["outputkm"]], reverse=True)

        # A reference to the model's starting time (i.e. when flushing is over)
        self.start_time = start_time
        self.stop_time = stop_time

        # run_type: "temperature", "solar", or "hydraulics"
        self.run_type = run_type

        # Empty dictionary to hold filenames and descriptions for each of the output files
        desc = {}

        write_extra_solar_outputs = False

        if (run_type in ("temperature", "solar") and write_extra_solar_outputs):
            # These are solar related parameters that the model calculates
            # but are not written to output right now. 
            # Writing these increases run times and they were mostly for testing.
            # Since they are normally not used they are turned off. 
            # To turn on, change write_extra_solar_outputs = True
            # Maybe in the future the control file will hold a switch to turn these on.
            desc["Heat_DR1"] = "Direct Solar Radiation Flux above Topographic Features (watts/square meter)"
            desc["Heat_DR2"] = "Direct Solar Radiation Flux below Topographic Features (watts/square meter)"
            desc["Heat_DR4"] = "Direct Solar Radiation Flux above the Stream (watts/square meter)"
            desc["Heat_DR5"] = "Direct Solar Radiation Flux Entering Stream (watts/square meter)"
            desc["Heat_DF1"] = "Diffuse Solar Radiation Flux above Topographic Features (watts/square meter)"
            desc["Heat_DF2"] = "Diffuse Solar Radiation Flux below Topographic Features (watts/square meter)"
            desc["Heat_DF4"] = "Diffuse Solar Radiation Flux above the Stream (watts/square meter)"
            desc["Heat_DF5"] = "Diffuse Solar Radiation Flux Entering Stream (watts/square meter)"
            desc["Solar_Azimuth"] = (
                "Solar Azimuth (degrees). Angle of the sun along the horizon with zero degrees "
                "corresponding to true north."
            )
            desc["Solar_Elevation"] = "Solar Elevation (degrees). Angle up from the horizon to the sun."

        if run_type in ("temperature", "solar"):
            # Solar, Temperature
            desc["Heat_SR1"] = "Solar Radiation Flux above Topographic Features (watts/square meter)"
            desc["Heat_SR2"] = "Solar Radiation Flux below Topographic Features (watts/square meter)"
            desc["Heat_SR3"] = "Solar Radiation Flux below Land Cover (watts/square meter)"
            desc["Heat_SR3b"] = "Solar Radiation Flux blocked by Land Cover (watts/square meter)"
            desc["Heat_SR4"] = "Solar Radiation Flux above the Stream (watts/square meter)"
            desc["Heat_SR5"] = "Solar Radiation Flux Entering Stream (watts/square meter)"
            desc["Shade"] = "Effective Shade"
            desc["VTS"] = "View to Sky. Calculated proportion of the sky hemisphere obscured by land cover."

        if run_type in ("temperature", "hydraulics"):
            # Hydro, Temperature
            desc["Hyd_DA"] = "Average Depth (meters)"
            desc["Hyd_DM"] = "Max Depth (meters)"
            desc["Hyd_Flow"] = "Flow Rate (cms)"
            desc["Hyd_Hyp"] = "Hyporheic Exchange (cms)"
            desc["Hyd_Vel"] = "Flow Velocity (meters/second)"
            desc["Hyd_WT"] = "Top Width (meters)"

        if run_type == "temperature":
            # Temperature
            desc["Heat_SR6"] = "Solar Radiation Flux Received by Stream (watts/square meter)"
            desc["Heat_SR7"] = "Solar Radiation Flux Received by Substrate (watts/square meter)"
            desc["Heat_Cond"] = "Streambed Conduction Flux (watts/square meter)"
            desc["Heat_Long"] = "Longwave Flux (watts/square meter)"
            desc["Heat_Conv"] = "Convection Flux (watts/square meter)"
            desc["Heat_Evap"] = "Evaporation Flux (watts/square meter)"
            desc["Rate_Evap"] = "Evaporation Rate (mm/hour)"
            desc["Temp_H2O"] = "Stream Temperature (Celsius)"
            desc["Temp_Sed"] = "Sediment Temperature (Celsius)"
            desc["Hyd_Disp"] = "Hydraulic Dispersion (square meters/second)"

        # Empty dictionary to hold the data.
        self.data = {}
        for name in list(desc.keys()):
            self.data[name] = {}

        # make a deepcopy of the empty variables dictionary for use later
        self.empty_vars = deepcopy(self.data)

        # Empty dictionary to store file objects
        self.files = {}

        # Build file objects and write headers
        for key in list(desc.keys()):
            header = [["File Created:"] + [ctime()]]
            header += [["Heat Source Version:"] + [self.params.get("version") or __version__]]
            header += [["Simulation Name:"] + [self.params.get("name", "")]]
            header += [["User Text:"] + [self.params.get("usertxt", "")]]
            header += [["Output:"] + [desc[key]]]
            header += [[""]]
            header += [["Datetime"]]

            if key in ["Heat_SR3b"]:
                header[6] += ["STREAM_KM"]

                if self.params["heatsource8"]:
                    # Same as 8 directions but no north
                    dir = ["NE", "E", "SE", "S", "SW", "W", "NW"]
                    zone = list(range(1, int(self.params["transsample_count"]) + 1))
                else:
                    dir = ["T" + str(x) for x in range(1, self.params["trans_count"] + 1)]
                    zone = list(range(1, int(self.params["transsample_count"]) + 1))

                if self.params["lcdatainput"] == "Values":
                    type_ = "HT"
                else:
                    type_ = "LC"

                for d in range(0, len(dir)):
                    for z in range(0, len(zone)):
                        header[6] += [type_ + "_" + dir[d] + "_S" + str(zone[z])]

                header[6] += ["Diffuse_Blocked"]
            else:
                header[6] += [("%0.3f" % x.km) for x in self.nodes]

            output_path = join(self.params["outputdir"], key + ".csv")
            try:
                self.files[key] = open(output_path, "w", newline="", encoding="utf-8")
            except PermissionError as exc:
                raise PermissionError(
                    f"Cannot write output file '{output_path}'. Check write permissions for this folder."
                ) from exc
            csv.writer(self.files[key]).writerows(header)
            self.files[key].flush()

    def __call__(self, time, hour, minute=0, second=0):
        """
        Store output data for one model time step.

        Skips writing before model start, stores relevant output variables for the run type,
        writes daily output at hour 23, and closes output files at model end.
        """
        # Ignore writing if we're still spinning up.
        if time < self.start_time:
            return

        # Create an Excel friendly time string
        timestamp = ("%0.7f" % float(time / 86400 + 25569))

        # Localize variables to save a bit of time
        nodes = self.nodes
        data = self.data

        write_extra_solar_outputs = False

        if (self.run_type in ("temperature", "solar") and write_extra_solar_outputs ):
            # These are solar related parameters that the model calculates
            # but are not written to output right now. 
            # Writing these increases run times and they were mostly for testing.
            # Since they are normally not used they are turned off. 
            # To turn on change write_extra_solar_outputs = True
            # Maybe in the future the control file will hold a switch to turn these on.
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
        if self.run_type in ("temperature", "solar"):
            data["Heat_SR1"][timestamp] = [x.F_Solar[1] for x in nodes]
            data["Heat_SR2"][timestamp] = [x.F_Solar[2] for x in nodes]
            data["Heat_SR3"][timestamp] = [x.F_Solar[3] for x in nodes]
            data["Heat_SR4"][timestamp] = [x.F_Solar[4] for x in nodes]
            data["Heat_SR5"][timestamp] = [x.F_Solar[5] for x in nodes]

        # Run only with hydro
        if self.run_type != "solar":
            data["Hyd_DA"][timestamp] = [(x.A / x.W_w) for x in nodes]
            data["Hyd_DM"][timestamp] = [x.d_w for x in nodes]
            data["Hyd_Flow"][timestamp] = [x.Q for x in nodes]
            data["Hyd_Hyp"][timestamp] = [x.Q_hyp for x in nodes]
            data["Hyd_Vel"][timestamp] = [x.U for x in nodes]
            data["Hyd_WT"][timestamp] = [x.W_w for x in nodes]

        # Run only with both solar and hydro
        if self.run_type == "temperature":
            data["Heat_SR6"][timestamp] = [x.F_Solar[6] for x in nodes]
            data["Heat_SR7"][timestamp] = [x.F_Solar[7] for x in nodes]
            data["Heat_Cond"][timestamp] = [x.F_Conduction for x in nodes]
            data["Heat_Conv"][timestamp] = [x.F_Convection for x in nodes]
            data["Heat_Evap"][timestamp] = [x.F_Evaporation for x in nodes]
            data["Heat_Long"][timestamp] = [x.F_Longwave for x in nodes]
            data["Rate_Evap"][timestamp] = [(x.E / x.dx / x.W_w * 3600 * 1000) for x in nodes]
            data["Temp_H2O"][timestamp] = [x.T for x in nodes]
            data["Temp_Sed"][timestamp] = [x.T_sed for x in nodes]
            data["Hyd_Disp"][timestamp] = [x.Disp for x in nodes]

        # Run the daily output on the last hour of the day
        if hour == 23:
            self.write_to_csv(self.run_type in ("temperature", "solar"), timestamp)

        if time >= self.stop_time:
            self.close()

    def daily(self, timestamp):
        """Compile and store daily output values from hourly model data."""
        nodes = self.nodes
        self.data["Shade"][timestamp] = [((x.F_DailySum[1] - x.F_DailySum[4]) / x.F_DailySum[1]) for x in nodes]
        self.data["VTS"][timestamp] = [x.ViewToSky for x in nodes]
        self.data["Heat_SR3b"][timestamp] = []
        # If there's no hour, we're at the beginning of a day, so we
        # write the values to a file.

    def write_to_csv(self, daily, timestamp):
        """
        Write model output data to csv files.
        Model output data is written for each configured output file, then the memory
        is reset for the next write cycle.
        """
        if daily:
            # for Shade, VTS, and Heat_SR3b, Not called for hydraulics
            self.daily(timestamp)

        data = self.data

        for name, openfile in list(self.files.items()):
            fileobj = csv.writer(openfile)
            timelist = sorted(data[name].keys())
            line = []

            if name in ["Heat_SR3b"]:
                timesteps = 86400.0 / float(self.params["dt"])
                for timestamp in timelist:
                    for i, x in enumerate(self.nodes):
                        line += [[timestamp] + [("%0.3f" % x.km)]]
                        if self.params["heatsource8"]:
                            # Seven directions
                            directions = [d for d in range(0, 7)]
                        else:
                            directions = [d for d in range(self.params["trans_count"])]

                        # solar blocked
                        for d in directions:
                            for zone in range(int(self.params["transsample_count"])):
                                daily_avg = x.Solar_Blocked[d][zone] / timesteps
                                line[i] += ["%0.4f" % daily_avg]
                        daily_avg_diffuse = x.Solar_Blocked["diffuse"] / timesteps
                        line[i] += ["%0.4f" % daily_avg_diffuse]
            else:
                for timestamp in timelist:
                    line += [[timestamp] + [("%0.4f" % x) for x in data[name][timestamp]]]

            fileobj.writerows(line)
            self.files[name].flush()

        del data
        self.data = deepcopy(self.empty_vars)

    def close(self):
        """Close all open output csv files."""
        for key in list(self.files.keys()):
            self.files[key].close()

"""This module collects model output values and writes them to CSV files."""


from copy import deepcopy
from os.path import join
from time import ctime
import csv
from math import nan

from heatsource9.__version__ import __version__
from heatsource9.domain.clock import excel_time


class OutputWriter:
    """Class to manage output writing during a model run."""

    def __init__(self, reach = None, start_time = None, stop_time = None, run_type = None, params = None):
        """Initialize an empty output or bind immediately if configuration is provided."""
        self._queued_dt_timesteps = 0
        self._write_threshold = None
        self.params = None
        self.nodes = []
        self.start_time = None
        self.stop_time = None
        self.run_type = None
        self.data = {}
        self.empty_vars = {}
        self.files = {}

        if reach is not None:
            self.bind(
                reach = reach,
                params = params,
                run_type = run_type,
                start_time = start_time,
                stop_time = stop_time,
            )

    def bind(self, reach, params, run_type, start_time = None, stop_time = None):
        """Store model output state and create output files for the bound run."""
        self.params = params

        # Store a sorted list of StreamNodes.
        if str(self.params["outputkm"]).lower() == "all":
            self.nodes = sorted(iter(list(reach.values())), reverse=True)
        else:
            # specific outputkm (already converted to reach keys by ModelSetup)
            self.nodes = sorted([reach[km] for km in reach if km in self.params["outputkm"]], reverse=True)

        # A reference to the model's starting time (i.e. when flushing is over)
        self.start_time = params["modelstart"] if start_time is None else start_time
        self.stop_time = params["modelend"] if stop_time is None else stop_time

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
            desc["Heat_SR1"] = "Solar Radiation Flux above Topographic Features (watts/square meter)"
            desc["Heat_SR2"] = "Solar Radiation Flux below Topographic Features (watts/square meter)"
            desc["Heat_SR3"] = "Solar Radiation Flux below Land Cover (watts/square meter)"
            desc["Heat_SR3b"] = "Solar Radiation Flux blocked by Land Cover (watts/square meter)"
            desc["Heat_SR4"] = "Solar Radiation Flux above the Stream (watts/square meter)"
            desc["Heat_SR5"] = "Solar Radiation Flux Entering Stream (watts/square meter)"
            desc["Shade"] = "Effective Shade"
            desc["VTS"] = "View to Sky. Calculated proportion of the sky hemisphere visible from the stream node after obstruction by near stream land cover and the stream bank."

        if run_type in ("temperature", "hydraulics"):
            desc["Hyd_DA"] = "Average Depth (meters)"
            desc["Hyd_DM"] = "Max Depth (meters)"
            desc["Hyd_Flow"] = "Flow Rate (cms)"
            desc["Hyd_Hyp"] = "Hyporheic Exchange (cms)"
            desc["Hyd_Vel"] = "Flow Velocity (meters/second)"
            desc["Hyd_WT"] = "Top Width (meters)"

        if run_type == "temperature":
            desc["Heat_SR6"] = "Solar Radiation Flux Received by Stream (watts/square meter)"
            desc["Heat_SR7"] = "Solar Radiation Flux Received by Substrate (watts/square meter)"
            desc["Heat_Cond"] = "Streambed Conduction Flux (watts/square meter)"
            desc["Heat_Long"] = "Longwave Flux (watts/square meter)"
            desc["Heat_Conv"] = "Convection Flux (watts/square meter)"
            desc["Heat_Evap"] = "Evaporation Flux (watts/square meter)"
            desc["Rate_Evap"] = "Evaporation Rate (mm/hour)"
            desc["Temp_H2O"] = "Stream Temperature (Celsius)"
            desc["Temp_Sed"] = "Sediment Temperature (Celsius)"
            desc["Temp_Hyp"] = "Hyporheic Return Water Temperature (Celsius)"
            desc["Hyd_Disp"] = "Hydraulic Dispersion (square meters/second)"

        # Empty dictionary to hold the data.
        self.data = {}
        for name in list(desc.keys()):
            self.data[name] = {}

        # make a deepcopy of the empty variables dictionary for use later
        self.empty_vars = deepcopy(self.data)

        # Empty dictionary to store file objects
        self._queued_dt_timesteps = 0
        self.files = {}

        # Build file objects and write headers
        for key in list(desc.keys()):
            header = self._headers_outputs(key, desc[key])

            output_path = join(self.params["outputdir"], key + ".csv")
            try:
                self.files[key] = open(output_path, "w", newline="", encoding="utf-8")
            except PermissionError as exc:
                raise PermissionError(
                    f"Cannot write output file '{output_path}'. Check write permissions for this folder."
                ) from exc
            csv.writer(self.files[key]).writerows(header)
            self.files[key].flush()

        self._write_threshold = self.write_threshold()

    def _headers_outputs(self, output_key, output_description):
        """Build csv header rows for one output file."""
        header = [["File Created:"] + [ctime()]]
        header += [["Heat Source Version:"] + [self.params.get("version") or __version__]]
        header += [["Simulation Name:"] + [self.params.get("name", "")]]
        header += [["User Text:"] + [self.params.get("usertxt", "")]]
        header += [["Output:"] + [output_description]]
        header += [[""]]
        header += [["Datetime"]]

        if output_key in ["Heat_SR3b"]:
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

        return header

    def write_threshold(self):
        """
        Return the number of queued non daily output timesteps allowed before
        queued outputs are written to csv.

        This threshold is used to limit problem causing memory growth when
        outputdt is low or the model has a lot of nodes to write to output.
        Either can cause a lot of data in the queue. The threshold is
        calculated from a target queued value
        count of 30,000,000, divided by the estimated number of non daily
        output values stored for one output timestep, which is:

            values per queued timestep = (node count) * (count of non daily output files to write)

        Smaller outputdt values cause output timesteps to be queued more often,
        so the threshold is reached sooner.

        The result is the maximum number of output timesteps that can be queued
        before write_to_csv() is called automatically.

        Daily outputs are Shade, VTS, and Heat_SR3b. Those are not included in
        the count because their size is fixed and the write schedule is at the
        end of each day.
        """
        target_values = 30000000
        nondaily_outputs = [
            name for name in self.data
            if name not in ("Shade", "VTS", "Heat_SR3b")
        ]
        values_per_timestep = len(self.nodes) * len(nondaily_outputs)

        if values_per_timestep <= 0:
            return 10

        return max(10, target_values // values_per_timestep)

    def queue_dt_outputs(self, time):
        """Queue non daily outputs for one model timestep."""
        # Create an Excel friendly time string
        timestamp = excel_time(time)
        nodes = self.nodes
        data = self.data

        write_extra_solar_outputs = False

        if (self.run_type in ("temperature", "solar") and write_extra_solar_outputs):
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

        if self.run_type in ("temperature", "solar"):
            data["Heat_SR1"][timestamp] = [x.F_Solar[1] for x in nodes]
            data["Heat_SR2"][timestamp] = [x.F_Solar[2] for x in nodes]
            data["Heat_SR3"][timestamp] = [x.F_Solar[3] for x in nodes]
            data["Heat_SR4"][timestamp] = [x.F_Solar[4] for x in nodes]
            data["Heat_SR5"][timestamp] = [x.F_Solar[5] for x in nodes]

        if self.run_type in ("temperature", "hydraulics"):
            data["Hyd_DA"][timestamp] = [(x.A / x.Ww) for x in nodes]
            data["Hyd_DM"][timestamp] = [x.Dw for x in nodes]
            data["Hyd_Flow"][timestamp] = [x.Q for x in nodes]
            data["Hyd_Hyp"][timestamp] = [x.Q_hyp for x in nodes]
            data["Hyd_Vel"][timestamp] = [x.U for x in nodes]
            data["Hyd_WT"][timestamp] = [x.Ww for x in nodes]

        if self.run_type == "temperature":
            data["Heat_SR6"][timestamp] = [x.F_Solar[6] for x in nodes]
            data["Heat_SR7"][timestamp] = [x.F_Solar[7] for x in nodes]
            data["Heat_Cond"][timestamp] = [x.F_Conduction for x in nodes]
            data["Heat_Conv"][timestamp] = [x.F_Convection for x in nodes]
            data["Heat_Evap"][timestamp] = [x.F_Evaporation for x in nodes]
            data["Heat_Long"][timestamp] = [x.F_LW for x in nodes]
            data["Rate_Evap"][timestamp] = [(x.Q_evap / x.dx / x.Ww * 3600 * 1000) for x in nodes]
            data["Temp_H2O"][timestamp] = [x.T for x in nodes]
            data["Temp_Sed"][timestamp] = [x.T_sed for x in nodes]
            data["Temp_Hyp"][timestamp] = [x.T_hyp if x.Q_hyp > 0 else nan for x in nodes]
            data["Hyd_Disp"][timestamp] = [x.Disp for x in nodes]

        self._queued_dt_timesteps += 1
        if self._queued_dt_timesteps >= self._write_threshold:
            self.write_to_csv()

    def queue_daily_outputs(self, time):
        """
        Queue daily output rows using the timestamp of the last included timestep.
        Daily outputs are Shade, VTS, and Heat_SR3b. These are not called for hydraulics.
        """

        timestamp = excel_time(time)
        nodes = self.nodes
        self.data["Shade"][timestamp] = [((x.F_DailySum[1] - x.F_DailySum[4]) / x.F_DailySum[1]) for x in nodes]
        self.data["VTS"][timestamp] = [x.ViewToSky for x in nodes]
        self.data["Heat_SR3b"][timestamp] = []

    def write_to_csv(self):
        """
        Write model output data to csv files.
        Model output data is written for each configured output file, then the memory
        is reset for the next write cycle.
        """
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
        self._queued_dt_timesteps = 0
        self.data = deepcopy(self.empty_vars)

    def close(self):
        """Close all open output csv files."""
        for key in list(self.files.keys()):
            self.files[key].close()

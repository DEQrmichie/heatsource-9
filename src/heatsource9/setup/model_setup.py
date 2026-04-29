import logging
from bisect import bisect
from calendar import timegm
from itertools import chain, repeat, count
from math import ceil, degrees, atan, exp
from pathlib import Path
from time import strptime

from heatsource9.io.console import print_console
from heatsource9.io.control_file import import_control_file

from heatsource9.setup.input_setup import InputSetup
from heatsource9.setup.site_setup import get_site_files
from heatsource9.setup.constants import KM_PRECISION, control_keys, dtype, head2var, sheetnames
from heatsource9.setup.setup_validation import align_rows_to_kmlist, validate_required_field
from heatsource9.__version__ import __version__
from heatsource9.domain.clock import Clock, pretty_time
from heatsource9.model.interpolator import Interpolator
from heatsource9.domain.simulation import Simulation
from heatsource9.model.streamnode import StreamNode

logger = logging.getLogger(__name__)


class ModelSetup(object):
    """
    ModelSetup contains methods to build a model and StreamNode instances from
    the input data.
    """

    def __init__(self, control_file_path, run_type):
        """Initialize model setup for the selected run type.
        """
        # run type value
        self.run_type = run_type

        self.params = {"run_type": self.run_type, "version": __version__}

        self.reach = {}
        self.ID2km = {}

        msg = "Starting Model Initialization"
        print_console(msg)

        self.inputs = InputSetup(control_file_path)

        # Boundary condition interpolators
        self.Q_bc = Interpolator()
        self.T_bc = Interpolator()

        # List of kilometers with met data nodes assigned.
        self.metDataSites = []
        self.met_site_rows = []
        self.trib_site_rows = []

        # Convenience variables
        self.dx = None
        self.multiple = None

        # Time lists (set during build)
        self.flowtimelist = ()
        self.continuoustimelist = ()
        self.flushtimelist = ()

    def build(self):
        """Build the model reach and return a Simulation(clock, nodes)."""

        # Load raw control file rows/params
        control = import_control_file(
            control_path=Path(self.inputs.model_dir) / self.inputs.control_file,
            dtype=dtype,
            control_sheet=sheetnames["controlfile"],
        )
        site_data = get_site_files(
            Path(self.inputs.model_dir) / self.inputs.control_file,
            control["control_params"],
            self.run_type,
        )
        merged_params = dict(control["control_params"])
        merged_params.update(site_data["met_params"])
        merged_params.update(site_data["trib_params"])
        self.met_site_rows = site_data["met_rows"]
        self.trib_site_rows = site_data["trib_rows"]
        self._parameterize_control_file(merged_params)

        self.inputs.params = self.params

        self.dx = self.params["dx"]
        self.multiple = int(self.dx / self.params["longsample"])

        # Setup time lists
        self.flowtimelist = self.get_timelist_unix()
        self.continuoustimelist = self.get_timelist_unix()
        self.flushtimelist = self.get_timelist_flush_period()

        # Start through the steps of building a reach full of StreamNodes
        if self.params["run_type"] == "temperature":
            # Temperature
            self.get_boundary_conditions()
            self.build_nodes()
            self.build_zones_w_codes()
            self.get_tributary_data()
            self.get_met_data()
            self.set_atmospheric_data()
            self.orient_nodes()

        elif self.params["run_type"] == "solar":
            # Solar
            self.get_boundary_conditions()
            self.build_nodes()
            self.build_zones_w_codes()
            self.get_met_data()
            self.set_atmospheric_data()
            self.orient_nodes()

        elif self.params["run_type"] == "hydraulics":
            # Hydraulics
            self.get_boundary_conditions()
            self.build_nodes()
            self.get_tributary_data()
            self.orient_nodes()

        # setup output km
        if not self.params["outputkm"] == "all":
            self.params["outputkm"] = self.get_locations("outputkm")

        msg = "Model Initialization Complete"
        print_console(msg)

        clock = Clock(
            start_time=self.params["flushtimestart"],
            end_time=self.params["modelend"],
            timestep_seconds=int(self.params["dt"]),
        )

        kms = sorted(list(self.reach.keys()), reverse=True)
        if not kms:
            raise ValueError("Model must contain at least one node")
        nodes = tuple(self.reach[km] for km in kms)

        return Simulation(clock=clock, nodes=nodes)

    def _parameterize_control_file(self, control_params):
        
        control_params = dict(control_params)

        # Validate control file values with run type requirements and defaults.
        keys = list(control_keys)
        keys.sort(reverse=True)

        for key in keys:
            value = control_params.get(key, None)
            # Control values are typed in import_control_file, then checked for required blanks here.
            self.params[key] = validate_required_field(
                run_type=self.params["run_type"],
                file_key="controlfile",
                field_name=key,
                value=value,
                params=control_params,
                source="controlfile",
            )
            control_params[key] = self.params[key]

        # Make dates into seconds since UTC epoch
        self.params["datastart"] = timegm(strptime(self.params["datastart"] + " 00:00:00", "%Y-%m-%d %H:%M:%S"))
        self.params["dataend"] = timegm(strptime(self.params["dataend"], "%Y-%m-%d")) + 86400

        if self.params["modelstart"] is None:
            self.params["modelstart"] = self.params["datastart"]
        else:
            self.params["modelstart"] = timegm(strptime(self.params["modelstart"] + " 00:00:00", "%Y-%m-%d %H:%M:%S"))

        if self.params["modelend"] is None:
            self.params["modelend"] = self.params["dataend"]
        else:
            self.params["modelend"] = timegm(strptime(self.params["modelend"], "%Y-%m-%d")) + 86400

        if self.params["run_type"] == "solar" and self.params["flushdays"] is None:
            self.params["flushdays"] = 0

        self.params["flushtimestart"] = self.params["modelstart"] - self.params["flushdays"] * 86400

        if self.params["run_type"] in ("solar", "temperature"):
            # If the number of transverse samples per direction is NOT reported, assume 4
            if not self.params["transsample_count"]:
                self.params["transsample_count"] = 4.0

            # Format for heat source 8 methods same as 8 directions but no north
            if self.params["heatsource8"]:
                self.params["trans_count"] = 7

            # Set the total number landcover sample count (0 = emergent)
            self.params["sample_count"] = int(self.params["transsample_count"] * self.params["trans_count"])
        else:
            # Hydraulics runs do not use land cover transect sampling values
            # for model calculations, but StreamNode setup expects integers.
            self.params["trans_count"] = int(self.params.get("trans_count") or 0)
            self.params["transsample_count"] = int(self.params.get("transsample_count") or 0)
            self.params["sample_count"] = 0

        # Set up evaporation method
        if self.params["evapmethod"] == "Penman":
            self.params["penman"] = True
        else:
            self.params["penman"] = False
            
        # convert dt values from minutes to seconds
        self.params["dt"] = int(round(self.params["dt"] * 60))
        self.params["outputdt"] = int(round(self.params["outputdt"] * 60))

        # make sure timestep divides into 60 minutes
        if 3600 % self.params["dt"] != 0:
            raise ValueError(
                "I'm sorry, your timestep ({0}) must evenly divide into 60 minutes.".format(self.params["dt"] / 60)
            )
        if self.params["outputdt"] < self.params["dt"]:
            msg = (
                "Output timestep (outputdt={0}) must be greater than or equal to model timestep (dt={1}).".format(
                    self.params["outputdt"] / 60,
                    self.params["dt"] / 60,
                )
            )
            raise ValueError(msg)
        if self.params["outputdt"] % self.params["dt"] != 0:
            msg = (
                "Output timestep (outputdt={0}) must be an exact multiple of model timestep (dt={1}).".format(
                    self.params["outputdt"] / 60,
                    self.params["dt"] / 60,
                )
            )
            raise ValueError(msg)

        if not isinstance(self.params["longsample"], int) or self.params["longsample"] <= 0:
            msg = "Longitudinal stream sample distance (longsample) must be an integer greater than zero."
            raise ValueError(msg)

        # dx must be a multiple of longsample and >= longsample
        if (self.params["dx"] % self.params["longsample"] or self.params["dx"] < self.params["longsample"]):
            raise ValueError("Distance step (dx) must be a multiple of the longitudinal stream sample distance")

        # Defaults if missing
        if not self.params.get("inputdir"):
            self.params["inputdir"] = self.inputs.model_dir
        if not self.params.get("outputdir"):
            raise ValueError("The control file must define 'outputdir'.")

    def orient_nodes(self):
            # Now we manually set each nodes next and previous 
            # kilometer values by stepping through the reach
            l = sorted(list(self.reach.keys()), reverse=True)
            # The headwater node
            head = self.reach[max(l)]
            # Set the previous and next kilometer of each node.
            slope_problems = []
            for i in range(len(l)):
                key = l[i] # The current node's key
                # Then, set pointers to the next and previous nodes
                if i == 0:
                    pass
                # At first node, there's no previous
                else:
                    self.reach[key].prev_km = self.reach[l[i - 1]]

                try:
                    self.reach[key].next_km = self.reach[l[i + 1]]
                except IndexError:
                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    # For last node (mouth) we set the downstream node equal
                    # to self, this is because we want to access the node's
                    # temp if there's no downstream, and this saves us an
                    # if statement.
                    self.reach[key].next_km = self.reach[key]
                # Set a headwater node
                self.reach[key].head = head
                self.reach[key].initialize()
                # check for a zero slope. We store all of them before 
                # checking so we can print a lengthy error that no-one 
                # will ever read.
                if self.run_type in ("temperature", "hydraulics"):
                    if self.reach[key].S is None or self.reach[key].S <= 0.0:
                        slope_problems.append(key)

            if self.run_type in ("temperature", "hydraulics"):  # zeros are alright in shade calculations
                if len(slope_problems):
                    raise Exception("The following reaches have zero slope. Kilometers: %s" % ",".join(
                        ['%0.3f' % i for i in slope_problems]))

    def set_atmospheric_data(self):
        """For each node without met data, use closest (up or downstream) node's data"""
        print_console("Assigning Meteorological Data to Nodes")

        # Localize the variable for speed
        sites = self.metDataSites

        # Sort the met site by km. This is necessary for 
        # the bisect module
        sites.sort()
        c = count()
        l = list(self.reach.keys())
        # This routine iterates through all nodes and uses bisect to 
        # determine which met site is closest to the node and 
        # initializes that node with the met data that is closest 
        # (up or downstream)
        for km, node in list(self.reach.items()):
            if km not in sites: # if we are not on the node where the 
                # met data is assigned we have t
                # Kilometer's downstream and upstream

                # zero is the lowest (protect against value of -1)
                lower = bisect(sites, km) - 1 if bisect(sites, km) - 1 > 0 else 0
                # bisect returns the length of a list when the bisecting 
                # number is greater than the greatest value.
                # Here we protect by max-ing out at the length of the 
                # list.
                upper = min([bisect(sites, km), len(sites) - 1])
                # Use the indexes to get the kilometers from the
                # sites list
                down = sites[lower]
                up = sites[upper]
                # Initialize to upstream's met data
                datasite = self.reach[up]
                # Only if the distance to the downstream node is closer 
                # do we use that
                if km - down < up - km:
                    datasite = self.reach[down]
                self.reach[km].metData = datasite.metData
                self.reach[km].zm = datasite.zm
            msg = "Assigning Node"
            current = next(c)+1
            print_console(msg, True, current, len(l))

    def get_boundary_conditions(self):
        """Get the boundary conditions"""
        # Get the columns, which is faster than accessing cells
        print_console("Reading boundary conditions")
        timelist = self.continuoustimelist

        # the data block is a tuple of tuples, each corresponding 
        # to a timestamp.
        if self.run_type == "solar":
            # Solar doesn't need a boundary condition
            data = [[0, 0] for i in timelist]
        else:
            data = self.inputs.import_bc()

        # Check out GetTributaryData() for details on this
        # reformatting of the data for the progress bar
        length = len(data)
        c = count()
        # Now set the discharge and temperature 
        # boundary condition dictionaries.

        for i in range(len(timelist)):
            time = timelist[i]
            flow = data[i][0]
            temp = data[i][1]

            # Get the flow boundary condition
            if flow == 0 or not flow:
                if self.run_type != "solar":
                    raise Exception("Missing flow boundary condition at model date/time %s." % pretty_time(time))
                else:
                    flow = 0
            self.Q_bc[time] = flow
            # Temperature boundary condition
            self.T_bc[time] = temp
            msg = "Reading boundary conditions"
            current = next(c)+1
            print_console(msg, True, current, length)

        # Next we expand or revise the dictionary to account for the 
        # flush period
        # Flush flow: model start value over entire flush period
        for i in range(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            self.Q_bc[time] = self.Q_bc[self.params["modelstart"]]
        # Flush temperature: first 24 hours repeated over flush period
        first_day_time = self.params["modelstart"]
        second_day = self.params["modelstart"] + 86400
        for i in range(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            self.T_bc[time] = self.T_bc[first_day_time]
            first_day_time += 3600
            if first_day_time >= second_day:
                first_day_time = self.params["modelstart"]

        self.Q_bc = self.Q_bc.view(self.params["flushtimestart"], self.params["modelend"], aft=1)
        self.T_bc = self.T_bc.view(self.params["flushtimestart"], self.params["modelend"], aft=1)

    def get_locations(self, ini):
        """Build a list of kilometers corresponding to the ini parameter
        that is passed.
    
        ini can equal: "tribkm", "metkm", or "outputkm"
        corresponding to tributary inflow sites, met data sites,
        or the model output kilometers"""

        t = ()
        l = list(self.reach.keys())
        l.sort()

        if (ini == "metkm" or
                ini == "outputkm" or
                self.params["tribsites"] > 0):
            # get a list of sites by km
            kms = self.params[ini].split(",")

            # remove spaces and make float
            kms = tuple([float(line.strip()) for line in kms])
        else:
            kms = tuple([])

        for site in range(0,len(kms)):
            km = kms[site]
            key = bisect(l, km) - 1
            t += l[key],  # Index by kilometer
        return t

    def get_stream_km_list(self):
        """Build a list of stream kilometers sorted from
        headwaters to mouth"""

        # Since the stream length is formatted as a floating point number 
        # we can sometimes run into trouble when calculating the number 
        # of nodes with dx due to ceil function rounding up on 
        # floating points that are not exact representations of the 
        # input value. Therefore we enforce a precision only up to the 
        # ten thousandths decimal place to make sure the number of nodes
        # is correct.        
        num_nodes = int(ceil(round(self.params["length"] * 1000 / (self.params["longsample"]), 4))) + 1
        precision_digits = abs(KM_PRECISION.as_tuple().exponent)
        kmlist = [
            round((node * self.params["longsample"]) / 1000, precision_digits) for node in range(0, num_nodes)
        ]
        kmlist.sort(reverse=True)
        return kmlist

    def get_timelist_unix(self):
        """Build a UNIX time list of floating point time values
        corresponding to the data start and end dates available in
        the control file"""
        timelist = []
        # hourly timestep
        timelist = list(range(self.params["datastart"], self.params["dataend"] + 60, 3600))
        return tuple(timelist)

    def get_timelist_flush_period(self):
        """Build a UNIX time list that represents the flushing period"""
        # This assumes that data is hourly, not tested with
        # variable input timesteps
        flushtimelist = []
        flushtime = self.params["flushtimestart"]
        while flushtime < self.params["modelstart"]:
            flushtimelist += flushtime,
            flushtime += 3600
        return tuple(flushtimelist)

    def get_tributary_data(self):
        """Populate the tributary flow and temperature values for
        nodes from the Flow Data page"""
        print_console("Reading inflow data")
        # Get a list of the timestamps that we have data for, and use 
        # that to grab the data block
        timelist = self.flowtimelist
        data = []
        if self.params["tribsites"] > 0:
            data = self.inputs.import_inflow()

        # The data is being put into this format
        # | Site 1     | Site 2     | Site 3       | ...
        # [((0.3, 15.7), (0.3, 17.7), (0.02, 18.2)), (, ...))]
        # To facilitate each site having it's own two item tuple.
        # The calls to tuple() just ensure that we are not making lists, 
        # which can be changed accidentally. Without them, the line is 
        # easier to understand
        data = [tuple(zip(line[0:None:2], line[1:None:2])) for line in data]

        # Get a tuple of kilometers to use as keys to the 
        # location of each tributary 
        kms = self.get_locations("tribkm")

        length = len(timelist)
        # Which datapoint time are we recording
        tm = count()

        # Quick list of nodes with flow data
        nodelist = []

        if self.params["tribsites"] > 0:
            for time in timelist:
                line = data.pop(0)
                # Error checking?! Naw!!
                c = count()
                for flow, temp in line:
                    i = next(c)
                    node = self.reach[kms[i]] # Index by kilometer
                    if node not in nodelist or not len(nodelist): 
                        nodelist.append(node)
                    if flow is None or (flow > 0 and temp is None):
                        raise Exception("Cannot have a tributary with \
                        blank flow or temperature conditions")
                    # Here, we actually set the tribs library, appending 
                    # to a tuple. Q_ and T_tribs are tuples of values 
                    # because we may have more than one input for a 
                    # given node

                    # Append to tuple
                    node.Q_tribs[time] += flow,
                    node.T_tribs[time] += temp,
                    msg = "Reading inflow data"
                    current = next(tm) + 1
                    print_console(msg, True, current, length * self.params["tribsites"])

        # Next we expand or revise the dictionary to account for the 
        # flush period
        # Flush flow: model start value over entire flush period
        for i in range(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            for node in nodelist:
                node.Q_tribs[time] = node.Q_tribs[self.params["modelstart"]]
        # Flush temperature: first 24 hours repeated over flush period
        first_day_time = self.params["modelstart"]
        second_day = self.params["modelstart"] + 86400
        for i in range(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            for node in nodelist:
                node.T_tribs[time] = node.T_tribs[first_day_time]
            first_day_time += 3600
            if first_day_time >= second_day:
                first_day_time = self.params["modelstart"]

        # Now we strip out the unnecessary values from the dictionaries. 
        # This is placed here at the end so we can dispose of it 
        # easily if necessary
        for node in nodelist:
            node.Q_tribs = node.Q_tribs.view(self.params["flushtimestart"], self.params["modelend"], aft=1)
            node.T_tribs = node.T_tribs.view(self.params["flushtimestart"], self.params["modelend"], aft=1)

    def get_met_data(self):
        """Get data from the input met data file"""
        # This is remarkably similar to GetInflowData. We get a block 
        # of data, then set the dictionary of the node
        print_console("Reading meteorological data")

        timelist = self.continuoustimelist

        metdata = self.inputs.import_met()

        data = [tuple(zip(line[0:None:4], line[1:None:4], line[2:None:4], line[3:None:4])) for line in metdata]

        # Get a tuple of kilometers to use as keys to the location of 
        # each met node
        kms = self.get_locations("metkm")
        if self.met_site_rows:
            metheights = [row["zm"] for row in self.met_site_rows]
        elif self.params.get("metheights"):
            metheights = [float(x.strip()) for x in self.params["metheights"].split(",")]
        else:
            metheights = [2.0 for i in kms]

        tm = count()  # Which datapoint time are we recording
        length = len(timelist)
        for time in timelist:
            line = data.pop(0)
            c = count()
            for cloud, Uzm, humidity, T_air in line:
                i = next(c)
            
                # Index by kilometer
                node = self.reach[kms[i]]
                node.zm = metheights[i]
                # Append this node to a list of all nodes which 
                # have met data
                if node.km not in self.metDataSites:
                    self.metDataSites.append(node.km)
                node.metData[time] = cloud, Uzm, humidity, T_air

            msg = "Reading meteorological data"
            current = next(tm) + 1
            print_console(msg, True, current, length)

        # Flush meteorology: first 24 hours repeated over flush period
        first_day_time = self.params["modelstart"]
        second_day = self.params["modelstart"] + 86400
        for i in range(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            for km in self.metDataSites:
                node = self.reach[km]
                node.metData[time] = node.metData[first_day_time]
            first_day_time += 3600
            if first_day_time >= second_day:
                first_day_time = self.params["modelstart"]

        # Now we strip out the met data outside the model period 
        # from the dictionaries. This is placed here
        # at the end so we can dispose of it easily if necessary

        tm = count()
        length = len(self.metDataSites)
        for km in self.metDataSites:
            node = self.reach[km]
            node.metData = node.metData.view(self.params["flushtimestart"], self.params["modelend"], aft=1)
            msg = "Subsetting met data"
            current = next(tm)+1
            print_console(msg, True, current, length)

    def zipper(self, iterable, mul=2):
        """Zippify list by grouping <mul> consecutive elements together

        Zipper returns a list of lists where the internal lists are
        groups of <mul> consecutive elements from the input list.
    
        For example:
        >>> lst = [0,1,2,3,4,5,6,7,8,9]
        >>> zipper(lst)
        [[0],[1,2],[3,4][5,6],[7,8],[9]]
        The first element is a length 1 list because we assume that it is
        a single element node (headwaters). Note that the last element,
        9, is alone as well, this method will figure out when there are
        not enough elements to make n equal length lists, and modify
        itself appropriately so that the remaining list will contain all
        leftover elements. The usefulness of this method is that it will
        allow us to average over each <mul> consecutive elements
        """
        # From itertools recipes... We use all but the 
        # first (boundary node) element
        lst = [i for i in zip(*[chain(iterable[1:], repeat(None, mul - 1))] * mul)]
        # Then we tack on the boundary node element
        lst.insert(0, (iterable[0],))
        # Then strip off the None values from the last (if any)
        lst[-1] = tuple([x for x in lst[-1] if x is not None])
        return self.numify(lst)

    def numify(self, lst):
        """Take a list of iterables and remove all values of None or empty strings"""
        # Remove None values at the end of each individual list
        for i in range(len(lst)):
            # strip out values of None from the tuple,
            # returning a new tuple
            lst[i] = [x for x in [x for x in lst[i] if x is not None]]
        # Remove blank strings from within the list
        for l in lst:
            n = []
            for i in range(len(l)):
                if l[i] == "":
                    n.append(i)
            n.reverse()
            for i in n:
                del l[i]
        # Make sure there are no zero length lists because they'll
        # fail if we average
        for i in range(len(lst)):
            if len(lst[i]) == 0:
                lst[i].append(0.0)
        return lst

    def multiplier(self, iterable, predicate=lambda x: x):
        """Return an iterable that was run through the zipper

        Take an iterable and strip the values of None, then send to
        the zipper and apply predicate to each value returned
        (zipper returns a list)"""

        # This is a way to safely apply a generic lambda function to an 
        # iterable. If I were paying attention to design, instead of 
        # just hacking away, I would have done this with decorators to 
        # modify the function. Now I'm too lazy to re-write it (well, 
        # not lazy, but I'm not paid as a programmer, and so I have
        # "better" things to do than optimize our code.) First we strip 
        # off the None values.
    
        strip_none = lambda y: [i for i in [x for x in y if x is not None]]
        return [predicate(strip_none(x)) for x in self.zipper(iterable, self.multiple)]

    def get_columnar_data(self):
        """
        Return a dictionary of input attributes that are
        averaged or summed as appropriate
        """
        # columns we grab from the inputs
        lc = ["STREAM_ID", "NODE_ID", "STREAM_KM", "LONGITUDE", "LATITUDE"]

        morph = ["STREAM_ID", "NODE_ID", "STREAM_KM",
                 "ELEVATION", "GRADIENT", "BOTTOM_WIDTH",
                 "CHANNEL_ANGLE_Z", "MANNINGS_n",
                 "SED_THERMAL_CONDUCTIVITY",
                 "SED_THERMAL_DIFFUSIVITY",
                 "SED_HYPORHEIC_THICKNESS", "HYPORHEIC_PERCENT",
                 "POROSITY", "Q_cont", "d_cont"]

        flow = ["INFLOW", "TEMPERATURE", "OUTFLOW"]

        # Operator methods to combine km values (named as model variables)
        # sums = ["Q_hyp_frac", "Q_accr", "Q_with"]
        # mins = ["km"]
        # aves = ["longitude", "latitude", "elevation", "S",
        #         "Wb", "Z", "n",
        #        "Ksed", "Alpha_sed", "Dsed",
        #        "phi", "Q_cont", "d_cont", "T_accr"]

        # Operator methods to combine row values (named as input column names)
        # sums = ["HYPORHEIC_PERCENT", "INFLOW", "OUTFLOW"]
        # mins = ["STREAM_KM"]
        # aves = ["LONGITUDE", "LATITUDE", "ELEVATION", "GRADIENT",
        #         "BOTTOM_WIDTH","CHANNEL_ANGLE_Z","MANNINGS_n",
        #         "SED_THERMAL_CONDUCTIVITY", "SED_THERMAL_DIFFUSIVITY", "SED_HYPORHEIC_THICKNESS",
        #         "POROSITY", "Q_cont","d_cont", "TEMPERATURE"]

        kmlist = self.get_stream_km_list()

        data = {}

        # Read data into a dictionary
        if self.run_type == "temperature":
            lcdata = self.inputs.import_lcdata(return_list=False)
            morphdata = self.inputs.import_morph(return_list=False)
            sums = ["Q_hyp_frac", "Q_accr", "Q_with"]
            lcdata = align_rows_to_kmlist(
                file_key = "lcdatafile",
                data_dict = lcdata,
                kmlist = kmlist,
                file_name = self.params["lcdatafile"],
                longsample = self.params["longsample"],
            )

            mins = ["km"]
            aves = ["longitude", "latitude", "Zs", "S", "Wb", "Z", "n",
                    "Ksed", "Alpha_sed", "Dsed", "Eta",
                    "Q_cont", "d_cont", "T_accr"]

        elif self.run_type == "solar":
            lcdata = self.inputs.import_lcdata(return_list=False)
            morphdata = self.inputs.import_morph(return_list=False)
            lcdata = align_rows_to_kmlist(
                file_key = "lcdatafile",
                data_dict = lcdata,
                kmlist = kmlist,
                file_name = self.params["lcdatafile"],
                longsample = self.params["longsample"],
            )
            sums = []
            mins = ["km"]
            aves = ["longitude", "latitude", "Zs"]

        elif self.run_type == "hydraulics":
            morphdata = self.inputs.import_morph(return_list=False)
            sums = ["Q_hyp_frac", "Q_accr", "Q_with"]
            mins = ["km"]
            aves = ["Zs", "S", "Wb", "Z", "n",
                    "Q_cont", "d_cont"]

        morphdata = align_rows_to_kmlist(
            file_key = "morphfile",
            data_dict = morphdata,
            kmlist = kmlist,
            file_name = self.params["morphfile"],
            longsample = self.params["longsample"],
        )

        # Add these columns to morph data since they do not exist in the input file.
        morphdata["Q_cont"] = [0.0 for km in kmlist]
        morphdata["d_cont"] = [0.0 for km in kmlist]

        if self.run_type in ["temperature", "hydraulics"]:
            if self.params.get("accretionfile"):
                accdata = align_rows_to_kmlist(
                    file_key = "accretionfile",
                    data_dict = self.inputs.import_accretion(),
                    kmlist = kmlist,
                    file_name = self.params["accretionfile"],
                    longsample = self.params["longsample"],
                )
            else:
                accdata = {
                    "INFLOW": [0.0 for km in kmlist],
                    "TEMPERATURE": [0.0 for km in kmlist],
                    "OUTFLOW": [0.0 for km in kmlist],
                }

        node_id_sources = {"morphfile": morphdata["NODE_ID"]}
        if self.run_type in ("temperature", "solar"):
            node_id_sources["lcdatafile"] = lcdata["NODE_ID"]
        if self.run_type in ("temperature", "hydraulics") and self.params.get("accretionfile"):
            node_id_sources["accretionfile"] = accdata["NODE_ID"]

        for i, km in enumerate(kmlist):
            node_ids = {name: values[i] for name, values in node_id_sources.items()}
            if len(set(node_ids.values())) > 1:
                parts = ["{0} NODE_ID={1}".format(name, node_ids[name]) for name in node_ids]
                msg = "NODE_ID mismatch for STREAM_KM {0}: {1}.".format(km, ", ".join(parts))
                raise ValueError(msg)

        # Add the values into the data dictionary but
        # we have to switch the key names because they are not 
        # consistent with model variable names. This needs to be fixed. 
        # TODO
        for k in morph:
            data[head2var[k]] = [i for i in morphdata[k]]

        if self.run_type in ["temperature", "solar"]:
            for k in lc:
                data[head2var[k]] = [i for i in lcdata[k]]

        if self.run_type in ["temperature", "hydraulics"]:
            for k in flow:
                data[head2var[k]] = [i for i in accdata[k]]

        # Then sum and average things as appropriate. 
        # multiplier() takes a tuple and applies the given lambda 
        # function to that tuple.
        for attr in sums:
            data[attr] = self.multiplier(data[attr], lambda x: sum(x))
        for attr in aves:
            data[attr] = self.multiplier(data[attr], lambda x: sum(x) / len(x))
        for attr in mins:
            data[attr] = self.multiplier(data[attr], lambda x: min(x))
        return data

    def build_nodes(self):
        # This is the worst of the methods but it works. # TODO
        print_console("Building Stream Nodes")
        Q_mb = 0.0

        # Grab all of the data in a dictionary
        data = self.get_columnar_data()

        # Build a boundary node
        node = StreamNode(run_type=self.run_type, Q_mb=Q_mb, run_params=self.params)
        # Then set the attributes for everything in the dictionary
        for k, v in list(data.items()):
            setattr(node, k, v[0])
        # set the flow and temp boundary conditions for the boundary node
        node.Q_bc = self.Q_bc
        node.T_bc = self.T_bc
        self.initialize_node(node)
        node.dx = self.params["longsample"]
        self.reach[node.km] = node
        self.ID2km[node.nodeID] = node.km

        # Figure out how many nodes we should have downstream. We use
        # math.ceil() because if we end up with a fraction, that means 
        # that there's a node at the end that is not a perfect multiple 
        # of the sample distance. We might end up ending at stream 
        # kilometer 0.5, for instance, in that case

        vars = (self.params["length"] * 1000) / self.params["longsample"]
        num_nodes = int(ceil(round((vars) / self.multiple, 4)))
        for i in range(0, num_nodes):
            node = StreamNode(run_type=self.run_type, Q_mb=Q_mb, run_params=self.params)
            for k, v in list(data.items()):
                # Add one to ignore boundary node
                setattr(node, k, v[i + 1])
            self.initialize_node(node)
            self.reach[node.km] = node
            self.ID2km[node.nodeID] = node.km
            msg = "Building Stream Nodes"
            print_console(msg, True, i + 1, num_nodes)

        # Find the mouth node and calculate the actual distance
        mouth = self.reach[min(self.reach.keys())]

        # number of extra variables if we're not perfectly divisible
        mouth_dx = (vars) % self.multiple or 1.0
        mouth.dx = self.params["longsample"] * mouth_dx

    def build_zones_w_codes(self):
        """Build zones when the landcover data files contains
        vegetation codes"""

        # Pull the LULC codes
        LCcodes = self.get_lc_codes()

        # Pull the LULC Data
        LCdata = self.inputs.import_lcdata(return_list=True)

        average = lambda x: sum(x) / len(x)
        transsample_count = self.params["transsample_count"]
        radial_count = self.params["trans_count"]
    
        keys = list(self.reach.keys())
    
        # Downstream sorted list of stream kilometers
        keys.sort(reverse=True)

        vheight = []
        vcanopy = []
        overhang = []
        cdepth = []
        elevation = []

        print_console("Translating landcover Data")
        if self.params["canopy_data"] == "LAI":
            # -------------------------------------------------------------
            # using LAI data

            k = []

            # For each column of LULC data
            for i in range(6, radial_count * transsample_count + 7):
            
                # LULC row and index column 
                col = [LCdata[row][i] for row in range(0, len(LCdata))]
                elev = [float(LCdata[row][i + radial_count * transsample_count]) for row in range(0, len(LCdata))]

                # Make a list from the LC codes from the column, then 
                # send that to the multiplier with a lambda function 
                # that averages them appropriately. Note, we're averaging
                # over the values (e.g. density) not the actual code, 
                # which would be meaningless.

                try:
                    vheight.append(self.multiplier([float(LCcodes[x][0])
                                                    for x in col],
                                                   average))
                    vcanopy.append(self.multiplier([float(LCcodes[x][1])
                                                    for x in col],
                                                   average))
                    k.append(self.multiplier([float(LCcodes[x][2])
                                              for x in col],
                                             average))
                    overhang.append(self.multiplier([float(LCcodes[x][3])
                                                     for x in col],
                                                    average))
                    cdepth.append(self.multiplier([float(LCcodes[x][4])
                                                   for x in col],
                                                  average))
                except KeyError as stderr:
                    raise Exception("At least one land cover code in %s is blank or not in %s (Code: %s)." % (
                        self.params["lcdatafile"], self.params["lccodefile"], stderr.message))
                if i > 6:
                    # There isn't a stream center elevation 
                    # (that is in the morphology file), so we don't want 
                    # to read in first elevation value which is actually 
                    # the last LULC col.

                    elevation.append(self.multiplier(elev, average))
                msg = "Translating Land Cover Data"
                print_console(msg, True, i, radial_count * transsample_count + 7)

            for i in range(len(keys)):
                node = self.reach[keys[i]]
                n = 0
                for tran in range(radial_count + 1):
                    for s in range(transsample_count):
                        node.lc_height_top[tran][s] = vheight[n][i]
                        node.lc_lai[tran][s] = vcanopy[n][i]
                        node.lc_canopy_cover[tran][s] = 0.0
                        node.lc_k[tran][s] = k[n][i]
                        node.lc_canopy_depth[tran][s] = cdepth[n][i]
                        node.lc_oh[tran][s] = overhang[n][i]

                        # 0 is emergent, there is only one value at s = 0
                        if tran == 0 and s == 0:
                            # Relative vegetation height is same as veg height
                            node.lc_height_node_top[tran][s] = vheight[n][i]
                            n = n + 1
                            # go to the next tran 
                            break
                        else:
                            # Vegetation height relative to the node, - 1 because there is no emergent elevation
                            node.lc_height_node_top[tran][s] = vheight[n][i] + (elevation[n - 1][i] - node.Zs)
                        n = n + 1
        else:
            # -------------------------------------------------------------
            # using canopy cover data

            # For each column of LULC data
            for i in range(6, radial_count * transsample_count + 7):

                # LULC row and index column 
                col = [LCdata[row][i] for row in range(0, len(LCdata))]
                elev = [float(LCdata[row][i + radial_count * transsample_count])
                        for row in range(0, len(LCdata))]
                # Make a list from the LC codes from the column, then 
                # send that to the multiplier with a lambda function 
                # that averages them appropriately. Note, we're 
                # averaging over the values (e.g. density) not the 
                # actual code, which would be meaningless.

                try:
                    vheight.append(self.multiplier([float(LCcodes[x][0])
                                                    for x in col],
                                                   average))
                    vcanopy.append(self.multiplier([float(LCcodes[x][1])
                                                    for x in col],
                                                   average))
                    overhang.append(self.multiplier([float(LCcodes[x][2])
                                                     for x in col],
                                                    average))
                    cdepth.append(self.multiplier([float(LCcodes[x][3])
                                                   for x in col],
                                                  average))

                except KeyError as stderr:
                    raise Exception("At least one land cover code in %s is blank or not in %s (Code: %s)." % (
                        self.params["lcdatafile"], self.params["lccodefile"], stderr.message))
                if i > 6:

                    # There isn't a stream center elevation
                    # (that is in the morphology file), so we don't want 
                    # to read in first elevation value which is actually 
                    # the last LULC col.

                    elevation.append(self.multiplier(elev, average))

                msg = "Translating Land Cover Data"
                print_console(msg, True, i, radial_count * transsample_count + 7)

            for i in range(len(keys)):
                node = self.reach[keys[i]]
                n = 0
                for tran in range(radial_count + 1):
                    for s in range(transsample_count):
                        node.lc_height_top[tran][s] = vheight[n][i]
                        node.lc_canopy_cover[tran][s] = vcanopy[n][i]
                        node.lc_lai[tran][s] = 0.0
                        node.lc_oh[tran][s] = overhang[n][i]
                        node.lc_canopy_depth[tran][s] = cdepth[n][i]

                        # 0 is emergent, there is only one value at s = 0
                        if tran == 0 and s == 0:
                            # Relative vegetation height is same as veg height
                            node.lc_height_node_top[tran][s] = vheight[n][i]
                            n = n + 1
                            # go to the next tran 
                            break
                        else:
                            # Vegetation height relative to the node, - 1 becaue there is no emergent elevation
                            node.lc_height_node_top[tran][s] = vheight[n][i] + (elevation[n - 1][i] - node.Zs)
                        n = n + 1

        # Average over the topo values
        # topo_ne = self.multiplier([float(LCdata[row][3]) for row in range(0, len(LCdata))], average)
        topo_e = self.multiplier([float(LCdata[row][5])
                                  for row in range(0, len(LCdata))],
                                 average)
        # topo_se = self.multiplier([float(LCdata[row][4]) for row in range(0, len(LCdata))], average)
        topo_s = self.multiplier([float(LCdata[row][4])
                                  for row in range(0, len(LCdata))],
                                 average)
        # topo_sw = self.multiplier([float(LCdata[row][7]) for row in range(0, len(LCdata))], average)
        topo_w = self.multiplier([float(LCdata[row][3])
                                  for row in range(0, len(LCdata))],
                                 average)
        # topo_nw = self.multiplier([float(LCdata[row][9]) for row in range(0, len(LCdata))], average)
        # topo_n = self.multiplier([float(LCdata[row][10]) for row in range(0, len(LCdata))], average)

        # ... and you thought things were crazy earlier! Here is where 
        # we build up the values for each node. This is culled from 
        # heat source version 7 VB code and discussions to try to 
        # simplify it... yeah, you read that right, simplify it... you 
        # should've seen it earlier!

        for h in range(len(keys)):
            msg = "Building land cover zones"
            print_console(msg, True, h + 1, len(keys))
            node = self.reach[keys[h]]
            vts_total = 0  # View to sky value

            # Now we set the topographic elevations in each direction
            # Topography factor Above Stream Surface
            node.TopoFactor = (topo_w[h] + topo_s[h] + topo_e[h]) / (90 * 3)
            # This is basically a list of directions, each 
            # with one of three topographies
            theta_topo_list = []
            angle_incr = 360.0 / radial_count
            dir_numbers = list(range(1, radial_count + 1))
            angle_mid = [x * angle_incr for x in dir_numbers]

            # Iterate through each transect direction
            for i in range(radial_count):
                dir_angle = angle_mid[i]
                if dir_angle < 135:
                    theta_topo_list.append(topo_e[h])
                elif dir_angle < 225:
                    theta_topo_list.append(topo_s[h])
                else:
                    theta_topo_list.append(topo_w[h])

            # Sun comes down and can be full-on, blocked by veg, or 
            # blocked by topography. Earlier implementations calculated 
            # each case on the fly. Here we chose a somewhat more elegant 
            # solution and calculate necessary angles. Basically, there 
            # is a minimum angle for which full sun is calculated 
            # (top of trees), and the maximum angle at which full shade 
            # is calculated (top of topography). Anything in between 
            # these is an angle for which sunlight is passing through 
            # trees. So, for each transect direction, we want to calculate these
            # two angles so that late we can test whether we are between 
            # them, and only do the shading calculations if that is true.

            # Iterate through each transect direction            
            for i in range(radial_count):

                # The minimum sun angle needed for full sun
                theta_full_sun = ()

                # The angle with longest path length in each veg zone
                theta_path = ()

                # Highest angle necessary for full shade in each veg zone
                theta_bank = ()

                # Numerator for the weighted Veg density calculation
                w_vdens_num = 0.0

                # Denominator for the weighted Veg density calculation
                w_vdens_dem = 0.0

                # Iterate through each of the zones
                for s in range(transsample_count):
                    v_height = vheight[i * transsample_count + s + 1][h]
                    v_can = vcanopy[i * transsample_count + s + 1][h]
                    v_overhang = overhang[i * transsample_count + s + 1][h]
                    elev = elevation[i * transsample_count + s][h]

                    # The below code only works when the samples start at stream edge.
                    # When that is implemented this can be turned on.
                    #if s == 0:  # TODO
                        # No overhang away from the stream
                        #v_overhang = 0

                    # Calculate the relative ground elevation. This is 
                    # the vertical distance from the stream surface to 
                    # the land surface
                    SH = elev - node.Zs
                    # Then calculate the relative vegetation height
                    VH = v_height + SH

                    # If lcsampmethod = point we assume you are sampling 
                    # a tree at a specific location rather than a veg 
                    # zone which represents the vegetation between two 
                    # sample points

                    if self.params["lcsampmethod"] == "zone":
                        adj_zone = 0.5
                    else:
                        adj_zone = 0.0

                    lc_distance1 = self.params["transsample_distance"] * (s + 1 - adj_zone)
                    lc_distance2 = self.params["transsample_distance"] * (s + 2 - adj_zone)

                    # We shift closer to the stream by the amount of overhang
                    if not s:
                        lc_distance1 -= v_overhang
                    if lc_distance1 <= 0:
                        lc_distance1 = 0.00001
                    # Calculate the minimum sun angle needed for full sun
                    # It gets added to a tuple of full sun values
                    theta_full_sun += degrees(atan(VH / lc_distance1)),

                    # Calculate angle with longest path length. This is used in the solar flux calcs
                    theta_path += degrees(atan(VH / lc_distance2)),

                    # Now get the maximum of bank shade and topographic 
                    # shade for this transect direction.
                    # likewise, a tuple of values
                    theta_bank += degrees(atan(SH / lc_distance1)),

                    # Calculate View To Sky
                    veg_angle = degrees(atan(VH / lc_distance1)) - degrees(atan(SH / lc_distance1))
                    if self.params["canopy_data"] == "LAI":
                        # use LAI data
                        Vk = k[i * transsample_count + s + 1][h]

                        # calculate canopy cover from LAI and k
                        LAI_can = 1 - exp(-Vk * v_can)
                        if LAI_can > 1:
                            LAI_can = 1
                        w_vdens_num += veg_angle * float(LAI_can)
                    else:
                        w_vdens_num += veg_angle * float(v_can)
                    w_vdens_dem += veg_angle

                    if s == transsample_count - 1:
                        if max(theta_full_sun) > 0:
                            # if bank and/or veg shade is occurring:
                            # Find weighted average the density:
                            # vdens_mod = (Amount of Veg shade * Veg dens) +
                            # (Amount of bank shade * bank dens, i.e. 1) / 
                            # (Sum of amount of shade)
                            if w_vdens_dem > 0:
                                vdens_ave_veg = w_vdens_num / w_vdens_dem
                            else:
                                vdens_ave_veg = 0
                            vdens_mod = ((max(theta_full_sun) - max(theta_bank)) * vdens_ave_veg + max(theta_bank)) / max(theta_full_sun)
                        else:
                            vdens_mod = 1.0
                        vts_total += max(theta_full_sun) * vdens_mod  # Add angle at end of each zone calculation
                theta_full_sun_max = max(theta_full_sun)
                theta_bank_max = max(theta_bank)
                node.ShaderList += (theta_full_sun_max, theta_topo_list[i], theta_bank_max, theta_full_sun, theta_path),
            node.ViewToSky = 1 - vts_total / (radial_count * 90)

    def get_lc_codes(self):
        """Return the codes from the Land Cover Codes input file
        as a dictionary of dictionaries"""

        data = self.inputs.import_lccodes()
        if self.params["canopy_data"] == "LAI":  # using LAI data

            # make a list of lists with values: [(height[0], lai[0], k[0], over[0], canopy_depth[0]), (height[1],...),...]
            vals = [tuple([float(j) for j in i]) for i in zip(data["HEIGHT"], data["LAI"], data["k"], data["OVERHANG"], data["CANOPY_DEPTH"])]
            codes = list(data["CODE"])  # CHECK
            data = {}

            for i, code in enumerate(codes):
                # Each code is a tuple in the form of (lc_height_top, lc_lai, lc_k, lc_oh, lc_canopy_depth)
                data[code] = vals[i]

        else:
            # make a list of lists with values: [(height[0], canopy[0], over[0], canopy_depth[0]), (height[1],...),...]
            vals = [tuple([float(j) for j in i]) for i in zip(data["HEIGHT"], data["CANOPY"], data["OVERHANG"], data["CANOPY_DEPTH"])]
            codes = list(data["CODE"])  # CHECK
            data = {}

            for i, code in enumerate(codes):
                # Each code is a tuple in the form of (lc_height_top, lc_canopy_cover, lc_oh, lc_canopy_depth)
                data[code] = vals[i]
                if vals[i][0] is not None and (vals[i][1] < 0 or vals[i][1] > 1):
                    raise ValueError(
                        "Canopy value of {0} in Land Cover Codes) must be >= 0.0 and <= 1.0".format(vals[i][1]))
        return data

    def initialize_node(self, node):
        """Perform some initialization of the StreamNode,
        and write some values to spreadsheet"""
        # Initialize each nodes tribs dictionary to a tuple
        for time in self.flowtimelist:
            node.Q_tribs[time] = ()
            node.T_tribs[time] = ()
        ##############################################################
        # Now that we have a stream node, we set the node's dx value, because
        # we have most nodes that are long-sample-distance times multiple,
        node.dx = self.params["dx"]  # Nodes distance step.
        node.dt = self.params["dt"]  # Set the node's timestep... this may have to be adjusted to comply with stability
        # Find the earliest temperature boundary condition
        mindate = min(self.T_bc.keys())
        if self.run_type == "hydraulics":  # Running hydraulics only
            node.T, node.T_prev, node.T_sed, node.T_hyp = 0.0, 0.0, 0.0, 0.0
        else:
            node.T = self.T_bc[mindate]
            node.T_prev = self.T_bc[mindate]
            node.T_sed = self.T_bc[mindate]
            node.T_hyp = self.T_bc[mindate]
        # we're in shadealator if the run type is "solar". Since much of the heat
        # math is coupled to the shade math, we have to make sure the hydraulic
        # values are not zero or blank because they'll raise ZeroDivisionError
        if self.run_type == "solar":
            for attr in ["Dw", "A", "Pw", "Ww", "U", "Disp", "Q_prev", "Q",
                         "Alpha_sed", "Dsed", "Ksed"]:
                if (getattr(node, attr) is None) or (getattr(node, attr) == 0):
                    setattr(node, attr, 0.01)
        if getattr(node, "F_DailySum", None) is None:
            node.F_DailySum = [0] * 5
        if getattr(node, "Solar_Blocked", None) is None:
            node.Solar_Blocked = {
                i: [0] * int(self.params["transsample_count"])
                for i in range(int(self.params["trans_count"]))
            }
            node.Solar_Blocked["diffuse"] = 0
        node.Q_hyp = 0.0  # Assume zero hyporheic flow unless otherwise calculated
        node.Q_evap = 0  # Same for evaporation

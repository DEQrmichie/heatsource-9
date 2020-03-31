# Heat Source, Copyright (C) 2000-2019, 
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

from __future__ import division

# Heat Source Methods
from ..Dieties.IniParamsDiety import IniParams
from ..Dieties.IniParamsDiety import head2var
from .Inputs import Inputs
from ..Stream.StreamNode import StreamNode
from ..Utils.Dictionaries import Interpolator
from ..Utils.Printer import Printer as print_console

# Builtin methods
from itertools import ifilter, izip, chain, repeat, count
from math import ceil, log, degrees, atan
from bisect import bisect
from time import ctime
import logging

logger = logging.getLogger(__name__)

class ModelSetup(object):
    """
    ModelSetup contains methods to build a model and
    StreamNode instances from the input data.
    """

    def __init__(self, model_dir, control_file, run_type=0):
        self.run_type = run_type
        self.reach = {}
        self.ID2km = {}

        msg = "Starting Model Initialization"
        logger.info(msg)
        print_console(msg)

        # create an input object to mange the data
        self.inputs = Inputs(model_dir, control_file)

        # read control file and parameterize IniParams
        self.inputs.import_control_file()

        # Make empty Dictionaries for the boundary conditions

        self.Q_bc = Interpolator()
        self.T_bc = Interpolator()

        # List of kilometers with met data nodes assigned.
        self.metDataSites = []

        # Some convenience variables
        self.dx = IniParams["dx"]

        # We have this many samples per distance step
        self.multiple = int(self.dx / IniParams["longsample"])

        # Setup for a model run
        # Get the list of model periods times
        self.flowtimelist = self.get_timelist_unix()
        self.continuoustimelist = self.get_timelist_unix()
        self.flushtimelist = self.get_timelist_flush_period()

        # Start through the steps of building a reach 
        # full of StreamNodes
        self.get_boundary_conditions()
        self.build_nodes()
        if IniParams["lcdatainput"] == "Values":
            self.build_zones_w_values()
        else:
            self.build_zones_w_codes()
        self.get_tributary_data()
        self.get_met_data()
        self.set_atmospheric_data()
        self.orient_nodes()

        # setup output km
        if not IniParams["outputkm"] == "all":
            IniParams["outputkm"] = self.get_locations("outputkm")

        msg = "Model Initialization Complete"
        logger.info(msg)
        print_console(msg)

    def orient_nodes(self):
        # Now we manually set each nodes next and previous 
        # kilometer values by stepping through the reach
        l = sorted(self.reach.keys(), reverse=True)
        # The headwater node
        head = self.reach[max(l)]
        # Set the previous and next kilometer of each node.
        slope_problems = []
        for i in xrange(len(l)):
            key = l[i]  # The current node's key
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
            if self.reach[key].S <= 0.0:
                slope_problems.append(key)

        if self.run_type != 1:  # zeros are alright in shade calculations
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
        l = self.reach.keys()
        # This routine iterates through all nodes and uses bisect to 
        # determine which met site is closest to the node and 
        # initializes that node with the met data that is closest 
        # (up or downstream)
        for km, node in self.reach.iteritems():
            if km not in sites:  # if we are not on the node where the
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
            msg = "Assigning Node"
            current = c.next() + 1
            logger.debug('{0} {1} {2}'.format(msg, True, current, len(l)))
            print_console(msg, True, current, len(l))

    def get_boundary_conditions(self):
        """Get the boundary conditions"""
        # Get the columns, which is faster than accessing cells
        print_console("Reading boundary conditions")
        timelist = self.continuoustimelist

        # the data block is a tuple of tuples, each corresponding 
        # to a timestamp.      
        data = self.inputs.import_bc()

        # Check out GetTributaryData() for details on this
        # reformatting of the data for the progress bar
        length = len(data)
        c = count()
        # Now set the discharge and temperature 
        # boundary condition dictionaries.

        for i in xrange(len(timelist)):
            time = timelist[i]
            flow = data[i][0]
            temp = data[i][1]

            # Get the flow boundary condition
            if flow == 0 or not flow:
                if self.run_type != 1:
                    raise Exception("Missing flow boundary condition for day %s " % ctime(time))
                else:
                    flow = 0
            self.Q_bc[time] = flow
            # Temperature boundary condition
            t_val = temp if temp is not None else 0.0
            self.T_bc[time] = t_val
            msg = "Reading boundary conditions"
            current = c.next() + 1
            logger.debug('{0} {1} {2}'.format(msg, current, length))
            print_console(msg, True, current, length)

        # Next we expand or revise the dictionary to account for the 
        # flush period
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

        self.Q_bc = self.Q_bc.view(IniParams["flushtimestart"], IniParams["modelend"], aft=1)
        self.T_bc = self.T_bc.view(IniParams["flushtimestart"], IniParams["modelend"], aft=1)

    def get_locations(self, ini):
        """Build a list of kilometers corresponding to the ini parameter
        that is passed.
        
        ini can equal: "inflowkm", "metkm", or "outputkm"
        corresponding to tributary inflow sites, met data sites,
        or the model output kilometers"""

        t = ()
        l = self.reach.keys()
        l.sort()

        if (ini == "metkm" or
                ini == "outputkm" or
                IniParams["inflowsites"] > 0):
            # get a list of sites by km
            kms = IniParams[ini].split(",")

            # remove spaces and make float
            kms = tuple([float(line.strip()) for line in kms])
        else:
            kms = tuple([])

        for site in xrange(0, len(kms)):
            km = kms[site]
            if km is None or not isinstance(km, float):
                # This is a bad dataset if there's no kilometer
                raise Exception("Must have a stream kilometer (e.g. 15.3) for each node in %s page!" % ini)
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
        num_nodes = int(ceil(round(IniParams["length"] * 1000 / (IniParams["longsample"]), 4))) + 1
        kmlist = []
        kmlist = [(node * IniParams["longsample"]) / 1000 for node in range(0, num_nodes)]
        kmlist.sort(reverse=True)
        return kmlist

    def get_timelist_unix(self):
        """Build a UNIX time list of floating point time values
        corresponding to the data start and end dates available in
        the control file"""
        timelist = []
        # hourly timestep
        timelist = range(IniParams["datastart"], IniParams["dataend"] + 60, 3600)
        return tuple(timelist)

    def get_timelist_flush_period(self):
        """Build a UNIX time list that represents the flushing period"""
        # This assumes that data is hourly, not tested with
        # variable input timesteps
        flushtimelist = []
        flushtime = IniParams["flushtimestart"]
        while flushtime < IniParams["modelstart"]:
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
        if IniParams["inflowsites"] > 0:
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
        kms = self.get_locations("inflowkm")

        length = len(timelist)
        # Which datapoint time are we recording
        tm = count()

        # Quick list of nodes with flow data
        nodelist = []

        if IniParams["inflowsites"] > 0:
            for time in timelist:
                line = data.pop(0)
                # Error checking?! Naw!!
                c = count()
                for flow, temp in line:
                    i = c.next()
                    node = self.reach[kms[i]]  # Index by kilometer
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
                    current = tm.next() + 1
                    logger.debug('{0} {1} {2}'.format(msg, current, length * IniParams["inflowsites"]))
                    print_console(msg, True, current, length * IniParams["inflowsites"])

        # Next we expand or revise the dictionary to account for the 
        # flush period
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

        # Now we strip out the unnecessary values from the dictionaries. 
        # This is placed here at the end so we can dispose of it 
        # easily if necessary
        for node in nodelist:
            node.Q_tribs = node.Q_tribs.view(IniParams["flushtimestart"], IniParams["modelend"], aft=1)
            node.T_tribs = node.T_tribs.view(IniParams["flushtimestart"], IniParams["modelend"], aft=1)

    def get_met_data(self):
        """Get data from the input met data csv file"""
        # This is remarkably similar to GetInflowData. We get a block 
        # of data, then set the dictionary of the node
        print_console("Reading meteorological data")

        timelist = self.continuoustimelist

        metdata = self.inputs.import_met()

        data = [tuple(zip(line[0:None:4], line[1:None:4], line[2:None:4], line[3:None:4])) for line in metdata]

        # Get a tuple of kilometers to use as keys to the location of 
        # each met node
        kms = self.get_locations("metkm")

        tm = count()  # Which datapoint time are we recording
        length = len(timelist)
        for time in timelist:
            line = data.pop(0)
            c = count()
            for cloud, wind, humidity, T_air in line:
                i = c.next()

                # Index by kilometer
                node = self.reach[kms[i]]
                # Append this node to a list of all nodes which 
                # have met data
                if node.km not in self.metDataSites:
                    self.metDataSites.append(node.km)
                # Perform some tests for data accuracy and validity
                if cloud is None:
                    cloud = 0.0
                if wind is None:
                    wind = 0.0
                if cloud < 0 or cloud > 1:
                    # Alright in shade-a-lator 
                    # # TODO zeros should not get a passed in 
                    # solar only runs, fix
                    if self.run_type == 1:
                        cloud = 0.0
                    else:
                        raise Exception(
                            "Cloudiness (value of '%s' in Meteorological Data) must be greater than zero and less "
                            "than one." % cloud)
                        # TODO RM fix this so it gives the km in exception - do for all
                if humidity < 0 or humidity is None or humidity > 1:
                    if self.run_type == 1:  # Alright in shade-a-lator
                        humidity = 0.0
                    else:
                        raise Exception(
                            "Humidity (value of '%s' in Meteorological Data) must be greater than zero and less than "
                            "one." % humidity)
                if T_air is None or T_air < -90 or T_air > 58:
                    if self.run_type == 1:  # Alright in shade-a-lator
                        T_air = 0.0
                    else:
                        raise Exception(
                            "Air temperature input (value of '%s' in Meteorological Data) outside of world records, "
                            "-89 to 58 deg C." % T_air)
                node.metData[time] = cloud, wind, humidity, T_air

            msg = "Reading meteorological data"
            current = tm.next() + 1
            logger.debug('{0} {1} {2}'.format(msg, current, length))
            print_console(msg, True, current, length)

        # Flush meteorology: first 24 hours repeated over flush period
        first_day_time = IniParams["modelstart"]
        second_day = IniParams["modelstart"] + 86400
        for i in xrange(len(self.flushtimelist)):
            time = self.flushtimelist[i]
            for km in self.metDataSites:
                node = self.reach[km]
                node.metData[time] = node.metData[first_day_time]
            first_day_time += 3600
            if first_day_time >= second_day:
                first_day_time = IniParams["modelstart"]

        # Now we strip out the met data outside the model period 
        # from the dictionaries. This is placed here
        # at the end so we can dispose of it easily if necessary

        tm = count()
        length = len(self.metDataSites)
        for km in self.metDataSites:
            node = self.reach[km]
            node.metData = node.metData.view(IniParams["flushtimestart"], IniParams["modelend"], aft=1)
            msg = "Subsetting met data"
            current = tm.next() + 1
            logger.info('{0} {1} {2}'.format(msg, current, length))
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
        lst = [i for i in izip(*[chain(iterable[1:], repeat(None, mul - 1))] * mul)]
        # Then we tack on the boundary node element
        lst.insert(0, (iterable[0],))
        # Then strip off the None values from the last (if any)
        lst[-1] = tuple(ifilter(lambda x: x is not None, lst[-1]))
        return self.numify(lst)

    def numify(self, lst):
        """Take a list of iterables and remove all values of None or empty strings"""
        # Remove None values at the end of each individual list
        for i in xrange(len(lst)):
            # strip out values of None from the tuple, 
            # returning a new tuple
            lst[i] = [x for x in ifilter(lambda x: x is not None, lst[i])]
        # Remove blank strings from within the list
        for l in lst:
            n = []
            for i in xrange(len(l)):
                if l[i] == "":
                    n.append(i)
            n.reverse()
            for i in n:
                del l[i]
        # Make sure there are no zero length lists because they'll 
        # fail if we average
        for i in xrange(len(lst)):
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

        strip_none = lambda y: [i for i in ifilter(lambda x: x is not None, y)]
        return [predicate(strip_none(x)) for x in self.zipper(iterable, self.multiple)]

    def get_columnar_data(self):
        """
        Return a dictionary of input attributes that are
        averaged or summed as appropriate
        """
        # columns from we grab from the inputs
        lc = ["STREAM_ID", "NODE_ID", "STREAM_KM", "LONGITUDE", "LATITUDE"]

        morph = ["ELEVATION", "GRADIENT", "BOTTOM_WIDTH",
                 "CHANNEL_ANGLE_Z", "MANNINGS_n",
                 "SED_THERMAL_CONDUCTIVITY",
                 "SED_THERMAL_DIFFUSIVITY",
                 "SED_HYPORHEIC_THICKNESSS", "HYPORHEIC_PERCENT",
                 "POROSITY", "Q_cont", "d_cont"]

        flow = ["INFLOW", "TEMPERATURE", "OUTFLOW"]

        # Ways that we grab the columns (named as model variables)
        # These are summed, not averaged
        sums = ["hyp_percent", "Q_in", "Q_out"]
        mins = ["km"]
        aves = ["longitude", "latitude", "elevation", "S", "W_b", "z", "n",
                "SedThermCond", "SedThermDiff", "SedDepth", "phi",
                "Q_cont", "d_cont", "T_in"]

        # sums = ["HYPORHEIC_PERCENT","INFLOW","OUTFLOW"]
        # mins = ["STREAM_KM"]
        # aves = ["LONGITUDE","LATITUDE","ELEVATION","GRADIENT",
        # "BOTTOM_WIDTH","CHANNEL_ANGLE_Z","MANNINGS_n",
        # "SED_THERMAL_CONDUCTIVITY",
        # "SED_THERMAL_DIFFUSIVITY", "SED_HYPORHEIC_THICKNESSS",
        # "POROSITY", "Q_cont","d_cont", "TEMPERATURE"]

        data = {}

        # Read data into a dictionary
        lcdata = self.inputs.import_lcdata(return_list=False)
        morphdata = self.inputs.import_morph(return_list=False)
        accdata = self.inputs.import_accretion()

        # Add these columns to morph data since they do not 
        # exist in the input file. 
        # TODO
        kmlist = self.get_stream_km_list()
        morphdata["Q_cont"] = [0.0 for km in kmlist]
        morphdata["d_cont"] = [0.0 for km in kmlist]

        # Add the values into the data dictionary but
        # we have to switch the key names because they are not 
        # consistent with model variable names. This needs to be fixed. 
        # TODO
        for k in lc:
            data[head2var[k]] = [i for i in lcdata[k]]

        for k in morph:
            data[head2var[k]] = [i for i in morphdata[k]]

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
        node = StreamNode(run_type=self.run_type, Q_mb=Q_mb)
        # Then set the attributes for everything in the dictionary
        for k, v in data.iteritems():
            setattr(node, k, v[0])
        # set the flow and temp boundary conditions for the boundary node
        node.Q_bc = self.Q_bc
        node.T_bc = self.T_bc
        self.initialize_node(node)
        node.dx = IniParams["longsample"]
        self.reach[node.km] = node
        self.ID2km[node.nodeID] = node.km

        # Figure out how many nodes we should have downstream. We use
        # math.ceil() because if we end up with a fraction, that means 
        # that there's a node at the end that is not a perfect multiple 
        # of the sample distance. We might end up ending at stream 
        # kilometer 0.5, for instance, in that case

        vars = (IniParams["length"] * 1000) / IniParams["longsample"]
        num_nodes = int(ceil(round((vars) / self.multiple, 4)))
        for i in range(0, num_nodes):
            node = StreamNode(run_type=self.run_type, Q_mb=Q_mb)
            for k, v in data.iteritems():
                setattr(node, k, v[i + 1])  # Add one to ignore boundary node
            self.initialize_node(node)
            self.reach[node.km] = node
            self.ID2km[node.nodeID] = node.km
            msg = "Building Stream Nodes"
            logger.debug('{0} {1} {2}'.format(msg, i + 1, num_nodes))
            print_console(msg, True, i + 1, num_nodes)

        # Find the mouth node and calculate the actual distance
        mouth = self.reach[min(self.reach.keys())]

        # number of extra variables if we're not perfectly divisible
        mouth_dx = (vars) % self.multiple or 1.0
        mouth.dx = IniParams["longsample"] * mouth_dx

    def build_zones_w_codes(self):
        """Build zones when the landcover data files contains
        vegetation codes"""

        # Pull the LULC codes
        LCcodes = self.get_lc_codes()

        # Pull the LULC Data
        LCdata = self.inputs.import_lcdata(return_list=True)

        average = lambda x: sum(x) / len(x)
        transsample_count = IniParams["transsample_count"]
        radial_count = IniParams["trans_count"]

        keys = self.reach.keys()

        # Downstream sorted list of stream kilometers
        keys.sort(reverse=True)

        vheight = []
        vcanopy = []
        overhang = []
        elevation = []

        print_console("Translating landcover Data")
        if IniParams["canopy_data"] == "LAI":
            # -------------------------------------------------------------
            # using LAI data

            k = []

            # For each column of LULC data
            for i in xrange(6, radial_count * transsample_count + 7):

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
                except KeyError, stderr:
                    raise Exception("At least one land cover code in %s is blank or not in %s (Code: %s)." % (
                        IniParams["lcdatafile"], IniParams["lccodefile"], stderr.message))
                if i > 6:
                    # There isn't a stream center elevation 
                    # (that is in the morphology file), so we don't want 
                    # to read in first elevation value which is actually 
                    # the last LULC col.

                    elevation.append(self.multiplier(elev, average))
                msg = "Translating Land Cover Data"
                logger.info('{0} {1} {2}'.format(msg, i, radial_count * transsample_count + 7))
                print_console(msg, True, i, radial_count * transsample_count + 7)

            for i in xrange(len(keys)):
                node = self.reach[keys[i]]
                n = 0
                for tran in xrange(radial_count + 1):
                    for s in xrange(transsample_count):
                        node.lc_height[tran][s] = vheight[n][i]
                        node.lc_canopy[tran][s] = vcanopy[n][i]
                        node.lc_k[tran][s] = k[n][i]
                        node.lc_oh[tran][s] = overhang[n][i]

                        # 0 is emergent, there is only one value at s = 0
                        if tran == 0 and s == 0:
                            # Relative vegetation height is same as veg height
                            node.lc_height_rel[tran][s] = vheight[n][i]
                            n = n + 1
                            # go to the next tran 
                            break
                        else:
                            # Vegetation height relative to the node, - 1 because there is no emergent elevation
                            node.lc_height_rel[tran][s] = vheight[n][i] + (elevation[n - 1][i] - node.elevation)
                        n = n + 1
        else:
            # -------------------------------------------------------------
            # using canopy cover data

            # For each column of LULC data
            for i in xrange(6, radial_count * transsample_count + 7):

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
                except KeyError, stderr:
                    raise Exception("At least one land cover code in %s is blank or not in %s (Code: %s)." % (
                        IniParams["lcdatafile"], IniParams["lccodefile"], stderr.message))
                if i > 6:
                    # There isn't a stream center elevation
                    # (that is in the morphology file), so we don't want 
                    # to read in first elevation value which is actually 
                    # the last LULC col.

                    elevation.append(self.multiplier(elev, average))

                msg = "Translating Land Cover Data"
                logger.info('{0} {1} {2}'.format(msg, i, radial_count * transsample_count + 7))
                print_console(msg, True, i, radial_count * transsample_count + 7)

            for i in xrange(len(keys)):
                node = self.reach[keys[i]]
                n = 0
                for tran in xrange(radial_count + 1):
                    for s in xrange(transsample_count):
                        node.lc_height[tran][s] = vheight[n][i]
                        node.lc_canopy[tran][s] = vcanopy[n][i]
                        node.lc_oh[tran][s] = overhang[n][i]

                        # 0 is emergent, there is only one value at s = 0
                        if tran == 0 and s == 0:
                            # Relative vegetation height is same as veg height
                            node.lc_height_rel[tran][s] = vheight[n][i]
                            n = n + 1
                            # go to the next tran 
                            break
                        else:
                            # Vegetation height relative to the node, - 1 becaue there is no emergent elevation
                            node.lc_height_rel[tran][s] = vheight[n][i] + (elevation[n - 1][i] - node.elevation)
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

        for h in xrange(len(keys)):
            msg = "Building land cover zones"
            logger.info('{0} {1} {2}'.format(msg, h + 1, len(keys)))
            print_console(msg, True, h + 1, len(keys))
            node = self.reach[keys[h]]
            vts_total = 0  # View to sky value
            # Now we set the topographic elevations in each direction
            # Topography factor Above Stream Surface
            node.TopoFactor = (topo_w[h] + topo_s[h] + topo_e[h]) / (90 * 3)
            # This is basically a list of directions, each 
            # with one of three topographies
            elevation_list = []
            angle_incr = 360.0 / radial_count
            dir_numbers = range(1, radial_count + 1)
            angle_mid = [x * angle_incr for x in dir_numbers]

            # Iterate through each transect direction
            for i in xrange(radial_count):
                dir_angle = angle_mid[i]
                if dir_angle < 135:
                    elevation_list.append(topo_e[h])
                elif dir_angle < 225:
                    elevation_list.append(topo_s[h])
                else:
                    elevation_list.append(topo_w[h])

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
            for i in xrange(radial_count):

                # The minimum sun angle needed for full sun
                t_full = ()

                # The angle with longest path length in each veg zone
                t_path = ()

                # Highest angle necessary for full shade in each veg zone
                t_none = ()

                # Numerator for the weighted Veg density calculation
                w_vdens_num = 0.0

                # Denominator for the weighted Veg density calculation
                w_vdens_dem = 0.0

                # Iterate through each of the zones
                for s in xrange(transsample_count):
                    v_height = vheight[i * transsample_count + s + 1][h]
                    v_can = vcanopy[i * transsample_count + s + 1][h]
                    v_overhang = overhang[i * transsample_count + s + 1][h]
                    elev = elevation[i * transsample_count + s][h]

                    # The below code only works when the samples start at stream edge.
                    # When that is implemented this can be turned on.
                    #if s == 0: # TODO
                        # No overhang away from the stream
                        #v_overhang = 0

                    # Calculate the relative ground elevation. This is 
                    # the vertical distance from the stream surface to 
                    # the land surface
                    SH = elev - node.elevation
                    # Then calculate the relative vegetation height
                    VH = v_height + SH

                    # If lcsampmethod = point we assume you are sampling 
                    # a tree at a specific location rather than a veg 
                    # zone which represents the vegetation between two 
                    # sample points

                    if IniParams["lcsampmethod"] == "zone":
                        adjust = 0.5
                    else:
                        adjust = 0.0

                    lc_distance1 = IniParams["transsample_distance"] * (s + 1 - adjust)
                    lc_distance2 = IniParams["transsample_distance"] * (s + 2 - adjust)

                    # We shift closer to the stream by the amount of overhang
                    if not s:
                        lc_distance1 -= v_overhang
                    if lc_distance1 <= 0:
                        lc_distance1 = 0.00001
                    # Calculate the minimum sun angle needed for full sun
                    # It gets added to a tuple of full sun values
                    t_full += degrees(atan(VH / lc_distance1)),

                    # Calculate angle with longest path length. This is used in the solar flux calcs
                    t_path += degrees(atan(VH / lc_distance2)),

                    # Now get the maximum of bank shade and topographic 
                    # shade for this transect direction.
                    # likewise, a tuple of values
                    t_none += degrees(atan(SH / lc_distance1)),

                    # Calculate View To Sky
                    veg_angle = degrees(atan(VH / lc_distance1)) - degrees(atan(SH / lc_distance1))
                    if IniParams["canopy_data"] == "LAI":
                        # use LAI data
                        Vk = k[i * transsample_count + s + 1][h]

                        # Purpose here is to calculate a LAI where 
                        # gap fraction = 0.1% (basically zero)
                        LAI_den = v_can / -log(0.001) / Vk
                        if LAI_den > 1:
                            LAI_den = 1
                        w_vdens_num += veg_angle * float(LAI_den)
                    else:
                        w_vdens_num += veg_angle * float(v_can)
                    w_vdens_dem += veg_angle

                    if s == transsample_count - 1:
                        if max(t_full) > 0:
                            # if bank and/or veg shade is occurring:
                            # Find weighted average the density:
                            # vdens_mod = (Amount of Veg shade * Veg dens) +
                            # (Amount of bank shade * bank dens, i.e. 1) / 
                            # (Sum of amount of shade)
                            if w_vdens_dem > 0:
                                vdens_ave_veg = w_vdens_num / w_vdens_dem
                            else:
                                vdens_ave_veg = 0
                            vdens_mod = ((max(t_full) - max(t_none)) * vdens_ave_veg + max(t_none)) / max(t_full)
                        else:
                            vdens_mod = 1.0
                        vts_total += max(t_full) * vdens_mod  # Add angle at end of each zone calculation
                node.ShaderList += (max(t_full), elevation_list[i], max(t_none), t_full, t_path),
            node.ViewToSky = 1 - vts_total / (radial_count * 90)

    def build_zones_w_values(self):
        """Build zones when the landcover data files contains explicit
        vegetation data instead of codes"""

        # Pull the LULC Data
        LCdata = self.inputs.import_lcdata(return_list=True)

        average = lambda x: sum(x) / len(x)
        transsample_count = IniParams["transsample_count"]
        radial_count = IniParams["trans_count"]
        shiftcol = radial_count * transsample_count  # Shift to get to each data type column

        keys = self.reach.keys()
        keys.sort(reverse=True)  # Downstream sorted list of stream kilometers

        vheight = []
        vcanopy = []
        overhang = []
        elevation = []

        print_console("Translating Land Cover Data")
        if IniParams["canopy_data"] == "LAI":
            # -------------------------------------------------------------
            # using LAI data

            k = []

            for i in xrange(6, shiftcol + 7):  # For each column of LULC data
                heightcol = [float(LCdata[row][i]) for row in range(0, len(LCdata))]
                elevcol = [float(LCdata[row][i + 1 + shiftcol]) for row in range(0, len(LCdata))]
                laicol = [float(LCdata[row][i + 1 + (shiftcol * 2)]) for row in range(0, len(LCdata))]
                kcol = [float(LCdata[row][i + 2 + (shiftcol * 3)]) for row in range(0, len(LCdata))]
                ohcol = [float(LCdata[row][i + 3 + (shiftcol * 4)]) for row in range(0, len(LCdata))]

                # Make a list from the LC codes from the column, then send 
                # that to the multiplier with a lambda function that averages 
                # them appropriately. Note, we're averaging over the values 
                # (e.g. density) not the actual code, which would be meaningless.
                try:
                    vheight.append(self.multiplier([float(x) for x in heightcol], average))
                    vcanopy.append(self.multiplier([float(x) for x in laicol], average))
                    k.append(self.multiplier([float(x) for x in kcol], average))
                    overhang.append(self.multiplier([float(x) for x in ohcol], average))

                except KeyError, stderr:
                    raise Exception("Vegetation height/density error" % stderr.message)
                if i > 6:
                    # There isn't a stream center elevation (that is in 
                    # the morphology file), so we don't want to read in first 
                    # elevation value which s actually the last LULC col.
                    elevation.append(self.multiplier(elevcol, average))
                msg = "Reading vegetation heights"
                logger.debug('{0} {1} {2}'.format(msg, i + 1, shiftcol + 7))
                print_console(msg, True, i + 1, shiftcol + 7)

            for i in xrange(len(keys)):
                node = self.reach[keys[i]]
                n = 0
                for tran in xrange(radial_count + 1):
                    for s in xrange(transsample_count):
                        node.lc_height[tran][s] = vheight[n][i]
                        node.lc_canopy[tran][s] = vcanopy[n][i]
                        node.lc_k[tran][s] = k[n][i]
                        node.lc_oh[tran][s] = overhang[n][i]

                        # 0 is emergent, there is only one value at s = 0
                        if tran == 0 and s == 0:
                            # Relative vegetation height is same as veg height
                            node.lc_height_rel[tran][s] = vheight[n][i]
                            n = n + 1
                            # go to the next tran 
                            break
                        else:
                            # Vegetation height relative to the node, - 1 becaue there is no emergent elevation
                            node.lc_height_rel[tran][s] = vheight[n][i] + (elevation[n - 1][i] - node.elevation)
                        n = n + 1

        else:
            # -------------------------------------------------------------
            # using canopy cover data

            for i in xrange(6, shiftcol + 7):  # For each column of LULC data
                heightcol = [float(LCdata[row][i]) for row in range(0, len(LCdata))]
                elevcol = [float(LCdata[row][i + 1 + shiftcol]) for row in range(0, len(LCdata))]
                dencol = [float(LCdata[row][i + 1 + (shiftcol * 2)]) for row in range(0, len(LCdata))]
                ohcol = [float(LCdata[row][i + 2 + (shiftcol * 3)]) for row in range(0, len(LCdata))]

                # Make a list from the LC codes from the column, then s
                # end that to the multiplier with a lambda function that
                # averages them appropriately. Note, we're averaging over
                # the values (e.g. density) not the actual code, which 
                # would be meaningless.
                try:
                    vheight.append(self.multiplier([float(x) for x in heightcol], average))
                    vcanopy.append(self.multiplier([float(x) for x in dencol], average))
                    overhang.append(self.multiplier([float(x) for x in ohcol], average))

                except KeyError, stderr:
                    raise Exception("Vegetation height/density error" % stderr.message)
                if i > 6:
                    # There isn't a stream center elevation (that is in 
                    # the morphology file), so we don't want to read 
                    # in first elevation value which s actually the 
                    # last LULC col.
                    elevation.append(self.multiplier(elevcol, average))

                msg = "Reading vegetation heights"
                logger.debug('{0} {1} {2}'.format(msg, i + 1, shiftcol + 7))
                print_console(msg, True, i + 1, shiftcol + 7)

            for i in xrange(len(keys)):
                node = self.reach[keys[i]]
                n = 0
                for tran in xrange(radial_count + 1):
                    for s in xrange(transsample_count):
                        node.lc_height[tran][s] = vheight[n][i]
                        node.lc_canopy[tran][s] = vcanopy[n][i]
                        node.lc_oh[tran][s] = overhang[n][i]

                        # 0 is emergent, there is only one value at s = 0
                        if tran == 0 and s == 0:
                            # Relative vegetation height is same as veg height
                            node.lc_height_rel[tran][s] = vheight[n][i]
                            n = n + 1
                            # go to the next tran 
                            break
                        else:
                            # Vegetation height relative to the node, - 1 becaue there is no emergent elevation
                            node.lc_height_rel[tran][s] = vheight[n][i] + (elevation[n - 1][i] - node.elevation)
                        n = n + 1

        # Average over the topo values
        # topo_ne = self.multiplier([float(LCdata[row][3]) for row in range(0, len(LCdata))], average)
        topo_e = self.multiplier([float(LCdata[row][5]) for row in range(0, len(LCdata))], average)
        # topo_se = self.multiplier([float(LCdata[row][4]) for row in range(0, len(LCdata))], average)
        topo_s = self.multiplier([float(LCdata[row][4]) for row in range(0, len(LCdata))], average)
        # topo_sw = self.multiplier([float(LCdata[row][7]) for row in range(0, len(LCdata))], average)
        topo_w = self.multiplier([float(LCdata[row][3]) for row in range(0, len(LCdata))], average)
        # topo_nw = self.multiplier([float(LCdata[row][9]) for row in range(0, len(LCdata))], average)
        # topo_n = self.multiplier([float(LCdata[row][10]) for row in range(0, len(LCdata))], average)        

        # ... and you thought things were crazy earlier! Here is where 
        # we build up the values for each node. This is culled from heat 
        # source version 7 VB code and discussions to try to simplify 
        # it... yeah, you read that right, simplify it... you should've 
        # seen it earlier!

        for h in xrange(len(keys)):
            msg = "Building VegZones"
            logger.info('{0} {1} {2}'.format(msg, h + 1, len(keys)))
            print_console(msg, True, h + 1, len(keys))

            node = self.reach[keys[h]]
            vts_total = 0  # View to sky value
            # Now we set the topographic elevations in each direction
            # Topography factor Above Stream Surface
            node.topo_factor = (topo_w[h] + topo_s[h] + topo_e[h]) / (90 * 3)
            # This is basically a list of directions, each with one 
            # of three topographies
            elevation_list = []
            angle_incr = 360.0 / radial_count
            dir_numbers = range(1, radial_count + 1)
            angle_mid = [x * angle_incr for x in dir_numbers]
            for i in xrange(radial_count):  # Iterate through each transect direction
                dir_angle = angle_mid[i]
                if dir_angle < 135:
                    elevation_list.append(topo_e[h])
                elif dir_angle < 225:
                    elevation_list.append(topo_s[h])
                else:
                    elevation_list.append(topo_w[h])
            # Sun comes down and can be full-on, blocked by veg, or blocked by topography. Earlier implementations
            # calculated each case on the fly. Here we chose a somewhat more elegant solution and calculate necessary
            # angles. Basically, there is a minimum angle for which full sun is calculated (top of trees), and the
            # maximum angle at which full shade is calculated (top of topography). Anything in between these is an
            # angle for which sunlight is passing through trees. So, for each direction, we want to calculate these
            # two angles so that late we can test whether we are between them, and only do the shading calculations
            # if that is true.

            for i in xrange(radial_count):  # Iterate through each transect direction
                # The minimum sun angle needed for full sun
                t_full = ()

                # The angle with longest path length in each veg zone
                t_path = ()

                # Highest angle necessary for full shade in each veg zone
                t_none = ()

                # Numerator for the weighted Veg density calculation
                w_vdens_num = 0.0

                # Denominator for the weighted Veg density calculation
                w_vdens_dem = 0.0

                for s in xrange(transsample_count):  # Iterate through each of the zones
                    v_height = vheight[i * transsample_count + s + 1][h]
                    if v_height < 0 or v_height is None or v_height > 120:
                        raise Exception(
                            "Vegetation height (value of %s in Land Cover Data) must be greater than zero and less "
                            "than 120 meters" % v_height)
                    v_can = vcanopy[i * transsample_count + s + 1][h]
                    v_overhang = overhang[i * transsample_count + s + 1][h]
                    elev = elevation[i * transsample_count + s][h]

                    # The below code only works when the samples start at stream edge.
                    # When that is implemented this can be turned on.
                    #if s == 0: # TODO
                        # No overhang away from the stream
                        #v_overhang = 0

                    # Calculate the relative ground elevation. This is the
                    # vertical distance from the stream surface to the land surface
                    SH = elev - node.elevation
                    # Then calculate the relative vegetation height
                    VH = v_height + SH

                    # Calculate the distance to the node from the 
                    # current landcover sample location.
                    # If lcsampmethod = point we assume you are sampling 
                    # a tree at a specific location rather than within a 
                    # zone which represents the vegetation between two 
                    # sample points
                    if IniParams["lcsampmethod"] == "zone":
                        adjust = 0.5
                    else:
                        adjust = 0.0
                    lc_distance1 = IniParams["transsample_distance"] * (
                            s + 1 - adjust)  # This is "+ 1" because s starts at 0
                    lc_distance2 = IniParams["transsample_distance"] * (
                            s + 2 - adjust)  # This is "+ 2" because we want to get to the farthest end of the zone
                    # We shift closer to the stream by the amount of overhang
                    # This is a rather ugly cludge.
                    if not s:
                        lc_distance1 -= v_overhang
                    if lc_distance1 <= 0:
                        lc_distance1 = 0.00001
                    # Calculate the minimum sun angle needed for full sun
                    t_full += degrees(atan(VH / lc_distance1)),  # It gets added to a tuple of full sun values

                    # Calculate angle with longest path length. This is used in the solar flux calcs
                    t_path += degrees(atan(VH / lc_distance2)),

                    # Now get the maximum of bank shade and topographic shade for this
                    # transect direction
                    t_none += degrees(atan(SH / lc_distance1)),  # likewise, a tuple of values

                    # Calculate View To Sky
                    veg_angle = degrees(atan(VH / lc_distance1)) - degrees(atan(SH / lc_distance1))
                    if IniParams["canopy_data"] == "LAI":
                        # use LAI data

                        Vk = k[i * transsample_count + s + 1][h]
                        # use LAI data
                        # Purpose here is to calculate a LAI where 
                        # gap fraction = 0.01% (basically zero)
                        LAI_den = v_can / -log(0.001) / Vk

                        if LAI_den > 1:
                            LAI_den = 1
                        w_vdens_num += veg_angle * float(LAI_den)
                    else:
                        w_vdens_num += veg_angle * float(v_can)
                    w_vdens_dem += veg_angle

                    if s == transsample_count - 1:
                        if max(t_full) > 0:
                            # if bank and/or veg shade is occurring:
                            # Find weighted average the density:
                            # vdens_mod = (Amount of Veg shade * Veg dens) +
                            # (Amount of bank shade * bank dens, i.e. 1) / 
                            # (Sum of amount of shade)
                            if w_vdens_dem > 0:
                                vdens_ave_veg = w_vdens_num / w_vdens_dem
                            else:
                                vdens_ave_veg = 0
                            vdens_mod = ((max(t_full) - max(t_none)) * vdens_ave_veg + max(t_none)) / max(t_full)
                        else:
                            vdens_mod = 1.0
                        vts_total += max(t_full) * vdens_mod  # Add angle at end of each zone calculation
                node.ShaderList += (max(t_full), elevation_list[i], max(t_none), t_full, t_path),
            node.ViewToSky = 1 - vts_total / (radial_count * 90)

    def get_lc_codes(self):
        """Return the codes from the Land Cover Codes csv input file
        as a dictionary of dictionaries"""

        data = self.inputs.import_lccodes()
        if IniParams["canopy_data"] == "LAI":  # using LAI data

            # make a list of lists with values: [(height[0], lai[0], k[0], over[0]), (height[1],...),...]
            vals = [tuple([float(j) for j in i]) for i in zip(data["HEIGHT"], data["LAI"], data["k"], data["OVERHANG"])]
            codes = list(data["CODE"])  # CHECK
            data = {}

            for i, code in enumerate(codes):
                # Each code is a tuple in the form of (lc_height, lc_canopy, lc_K, lc_oh)
                data[code] = vals[i]

        else:
            # make a list of lists with values: [(height[0], canopy[0], over[0]), (height[1],...),...]
            vals = [tuple([float(j) for j in i]) for i in zip(data["HEIGHT"], data["CANOPY"], data["OVERHANG"])]
            codes = list(data["CODE"])  # CHECK
            data = {}

            for i, code in enumerate(codes):
                # Each code is a tuple in the form of (lc_height, lc_canopy, lc_oh)
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
        node.dx = IniParams["dx"]  # Nodes distance step.
        node.dt = IniParams["dt"]  # Set the node's timestep... this may have to be adjusted to comply with stability
        # Find the earliest temperature boundary condition
        mindate = min(self.T_bc.keys())
        if self.run_type == 2:  # Running hydraulics only
            node.T, node.T_prev, node.T_sed = 0.0, 0.0, 0.0
        else:
            if self.T_bc[mindate] is None:
                # Shade-a-lator doesn't need a boundary condition
                if self.run_type == 1:
                    self.T_bc[mindate] = 0.0
                else:
                    raise Exception("Boundary temperature conditions cannot be blank")
            node.T = self.T_bc[mindate]
            node.T_prev = self.T_bc[mindate]
            node.T_sed = self.T_bc[mindate]
        # we're in shadealator if the runtype is 1. Since much of the heat
        # math is coupled to the shade math, we have to make sure the hydraulic
        # values are not zero or blank because they'll raise ZeroDivisionError
        if self.run_type == 1:
            for attr in ["d_w", "A", "P_w", "W_w", "U", "Disp", "Q_prev", "Q",
                         "SedThermDiff", "SedDepth", "SedThermCond"]:
                if (getattr(node, attr) is None) or (getattr(node, attr) == 0):
                    setattr(node, attr, 0.01)
        node.Q_hyp = 0.0  # Assume zero hyporheic flow unless otherwise calculated
        node.E = 0  # Same for evaporation

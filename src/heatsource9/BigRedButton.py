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

"""BigRedButton holds the main Heat Source model controls.
The ModelControl class loads and controls the model run.
A model instance object is created using the ModelSetup class.
"""
# Heat Source modules
from heatsource9.Dieties.IniParamsDiety import IniParams
from heatsource9.ModelSetup.ModelSetup import ModelSetup
from heatsource9.Dieties.ChronosDiety import Chronos
from heatsource9.Utils.Printer import Printer as print_console
from heatsource9.Utils.Output import Output as O
from heatsource9.__version__ import version_string
from heatsource9.ModelSetup.Inputs import Inputs

# Built-in modules
from builtins import next
from builtins import range
from builtins import object
import logging
from itertools import count
import traceback
import os
from glob import glob
import argparse
from time import time as Time

# set up logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    filename='heatsource.log',
                    filemode='w')
logger = logging.getLogger(__name__)


class ModelControl(object):
    """Main model control class for Heat Source.

    While it works, this class is basically a one-off hack
    that was designed as an interim solution to model control.
    As it is, it's a fairly dumb method which simply grabs
    a list of StreamNodes from the ModelSetup class,
    and loops through them using the iterator inside
    ChronosDiety to advance the time. Ideally, this would be
    a more sophisticated system that could hold and control
    multiple stream reaches, modeling upstream reaches
    and using the results as tributary inputs to the down-
    stream reaches. For that, we'd likely need a separate
    Reach class. Since this was essentially an interim
    solution to the problem, don't hesitate to improve it.
    """

    def __init__(self, model_dir, control_file, run_type=0):
        """
        model_dir is the path to the directory where the
        control file is located.
        
        control_file is the control file name. It must be a comma
        delimited text file.
        
        run_type is one of 0, 1, or 2 for Heat Source (Temperature),
        Solar only, hydraulics only respectively.
        """

        # Add run type and heat source model version into IniParams
        IniParams["version"] = version_string
        IniParams["run_type"] = run_type

        msg = ("\n")
        msg += "Heat Source Version:  {0}\n".format(IniParams["version"])
        msg += ("\n\n")
        msg += ("Oregon Department of Environmental Quality\n\n"
                "This program comes with ABSOLUTELY NO WARRANTY. "
                "Appropriate model \n"
                "use and application are the sole responsibility "
                "of the user. \n"
                "This is free software, and you are welcome to "
                "redistribute it under \n"
                "certain conditions described in the License.\n\n")

        print_console(msg)
        logger.info(msg)

        # Create a ModelSetup instance.
        self.HS = ModelSetup(model_dir, control_file, run_type)

        # This is the list of StreamNode instances- we sort it in reverse
        # order because we number stream kilometer from the mouth to the
        # headwater, but we want to run the model from headwater to mouth.
        self.reachlist = sorted(iter(list(self.HS.reach.values())), reverse=True)

        # This if statement prevents us from having to test every 
        # timestep. We just call self.run_all(), which is a classmethod 
        # pointing to the correct method.
        if run_type == 0:
            self.run_all = self.run_hs
        elif run_type == 1:
            self.run_all = self.run_sh
        elif run_type == 2:
            self.run_all = self.run_hy
        else:
            logger.error("Bad run_type: {0}. Must be 0, 1, or 2.\
            Something wrong with the executable".format(run_type))
            raise Exception("Bad run_type: {0}. Must be 0, 1, or 2. \
            Something wrong with the executable".format(run_type))

        # Create a Chronos iterator that controls all model time.
        Chronos.clock(start=IniParams["modelstart"],
                      stop=IniParams["modelend"],
                      dt=IniParams["dt"],
                      spin=IniParams["flushdays"],
                      offset=IniParams["offset"])

        # This is the output class, which is essentially just a list
        # of file objects and an append method which writes to them
        # every so often.
        self.Output = O(self.HS.reach, IniParams["modelstart"], run_type)

    def run(self):
        """Run the model one time

        Use the Chronos instance and list of StreamNodes to cycle
        through each timestep and spacestep, calling the appropriate
        StreamNode functions to calculate heat and hydraulics."""

        msg = "Starting Simulation:  {0}".format(IniParams["name"])
        logger.info(msg)
        print_console(msg)

        # Current time of the Chronos clock (i.e. this timestep)
        time = Chronos.TheTime

        # Stop time for Chronos
        stop = Chronos.stop

        # Start time for Chronos, model start, not flush/spin start.
        start = Chronos.start

        # flush in seconds
        flush = start - (IniParams["flushdays"] * 86400)
        # Number of timesteps is based on the division of the timesteps 
        # into the hour. In other words 1 day with a 1 minute dt is 
        # 1440 timesteps, while a 3 minute dt is only 480 timesteps. 
        # Thus, we define the timesteps by dividing dt (now in seconds) 
        # by 3600
        timesteps = (stop - flush) / IniParams["dt"]

        # Counter iterator for counting current timesteps passed
        cnt = count()

        # Volume of water flowing out of mouth (for simple mass balance)
        out = 0

        # Current computer time- for estimating total model runtime
        time1 = Time()

        # Localize run_type for a bit more speed
        ################################################################
        # So, it's simple and stupid. We basically just cycle through the
        # time until we get to the model stop time, and each cycle, we 
        # call a method that cycles through the StreamNodes, and 
        # calculates them in order. A smarter way would be to thread 
        # this so we can start calculating the second timestep (using 
        # another CPU or core) while the first one is still unfinished.
        while time <= stop:
            year, month, day, hour, minute, second, JD, offset, JDC = Chronos.time_tuple()
            # zero hour+minute+second means first timestep of new day
            if not (hour + minute + second):
                # zero out the daily flux sum at this point.
                for nd in self.reachlist:
                    nd.F_DailySum = [0] * 5
                    nd.Solar_Blocked = {}

                    # Number of radial sample directions
                    for i in range(IniParams["trans_count"]):
                        # A spot for each land cover sample
                        nd.Solar_Blocked[i] = [0] * IniParams["transsample_count"]
                    nd.Solar_Blocked["diffuse"] = 0

            self.run_all(time, hour, minute, second, JD, JDC)

            # If minute and second are both zero, we are at the top of 
            # the hour. 
            if (minute == 0 and second == 0):
                ts = next(cnt) # Number of actual timesteps per tick
                
                # Number of timesteps in one hour
                hr = 60 / (IniParams["dt"] / 60)
                # This writes a line to the status bar.
                msg = "Timesteps:"
                logger.info('{0} {1} {2}'.format(msg, (ts) * hr, timesteps))
                print_console(msg, True, (ts) * hr, timesteps)

                # Call the Output class to update the textfiles. We call
                # this every hour and store the data, then we write to 
                # file every day. Limiting disk access saves us 
                # considerable time.
                self.Output(time, hour)

            # ---- 
            # Uncomment to output every timestep and 
            # comment section above
            # ts = cnt.next()
            # msg = "Timesteps:"
            # logger.info('{0} {1} {2}'.format(msg, ts, timesteps))
            # print_console("Timesteps:", True, ts, timesteps)
            # self.Output(time, hour, minute, second)
            # ---- 

            # We've made it through the entire stream without an error, 
            # so we update our mass balance
            # by adding the discharge of the mouth...
            out += self.reachlist[-1].Q
            # and tell Chronos that we're moving time forward.
            time = Chronos(tick=True)

        # So, here we are at the end of a model run. First we 
        # calculate how long all of this took
        total_time = (Time() - time1) / 60
        # Calculate the mass balance inflow
        balances = [x.Q_mass for x in self.reachlist]
        total_inflow = sum(balances)
        # Calculate how long the model took to run each timestep for 
        # each spacestep. This is how long (on average) it took to 
        # calculate a single StreamNode's information one time. 
        # Ideally, for performance and impatience reasons, we want this
        # to be somewhere around or less than 1 microsecond.
        microseconds = (total_time / timesteps / len(self.reachlist)) * 1e6

        message = "Simulation Complete"

        message += "\n\nFinished in {0:0.1f} minutes ".format(total_time)
        message += "(spent {0:0.3f} microseconds ".format(microseconds)
        message += "in each stream node).\n"
        message += "Water Balance: {0:0.3f}/{1:0.3f}\n".format(total_inflow, out)
        message += "Simulation: {0}\n".format(IniParams["name"])
        message += "Outputs: {0}\n\n".format(IniParams["outputdir"])

        logger.info(message)
        print_console(message)
        # raw_input('Press <ENTER> to close this console')

    #############################################################
    # Three different versions of the run() routine, depending on the 
    # run_type.

    def run_hs(self, time, H, M, S, JD, JDC):
        """Call both hydraulic and solar routines for each StreamNode"""
        [x.CalcDischarge(time) for x in self.reachlist]
        [x.CalcHeat(time, H, M, S, JD, JDC) for x in self.reachlist]
        [x.maccormick2(time) for x in self.reachlist]

    def run_hy(self, time, H, M, S, JD, JDC):
        """Call hydraulic routines for each StreamNode"""
        [x.CalcDischarge(time) for x in self.reachlist]

    def run_sh(self, time, H, M, S, JD, JDC):
        """Call solar routines for each StreamNode"""
        [x.CalcHeat(time, H, M, S, JD, JDC, True) for x in self.reachlist]

def run(model_dir, control_file):
    """Run full temperature model"""
    try:
        HSP = ModelControl(model_dir, control_file, 0)
        HSP.run()
        del HSP
    except:
        msg = "Error: {0}".format(traceback.format_exc())
        logging.error(msg)
        print_console(msg)


def run_solar(model_dir, control_file):
    """Run solar routines only"""
    try:
        HSP = ModelControl(model_dir, control_file, 1)
        HSP.run()
    except:
        msg = "Error: {0}".format(traceback.format_exc())
        logging.error(msg)
        print_console(msg)


def run_hydraulics(model_dir, control_file):
    """Run hydraulics only"""
    try:
        HSP = ModelControl(model_dir, control_file, 2)
        HSP.run()
    except:
        msg = "Error: {0}".format(traceback.format_exc())
        logging.error(msg)
        print_console(msg)


def setup_cf(model_dir, control_file, use_timestamp=False, overwrite=False, **kwargs):
    """Write a blank control file or use **kwargs to parameterize it."""
    try:
        # create an input object
        inputs = Inputs(model_dir, control_file)
        
        # Write a blank control file
        inputs.parameterize_cf(use_timestamp=use_timestamp, overwrite=overwrite, **kwargs)
        
    except:
        msg = "Error: {0}".format(traceback.format_exc())
        logging.error(msg)
        print_console(msg)


def setup_mi(model_dir, control_file, use_timestamp=False, overwrite=False):
    """Write blank input files. Control file must already be parameterized."""
    try:
        # create an input object
        inputs = Inputs(model_dir, control_file)

        # Control file must already be parameterized
        inputs.import_control_file()
        
        # Write blank input files,
        inputs.setup(use_timestamp=use_timestamp, overwrite=overwrite)

    except:
        msg = "Error: {0}".format(traceback.format_exc())
        logging.error(msg)
        print_console(msg)


def hs():
    """Entry point to run Heat Source from command line."""

    pkg_dir = os.path.abspath(os.path.dirname(__file__))

    parser = argparse.ArgumentParser(prog='hs', prefix_chars='-',
                                     usage='''hs <command> [options]''')

    # general options
    parser.add_argument('-v', action='version', version='hs {0} installed at {1}.'.format(version_string, pkg_dir),
                        help='The heat source version and install directory.')
    parser.add_argument('-md', '--model-dir', nargs='?', default=os.getcwd() + '/',
                        help='Path to the model directory. If not used default is current working directory.')

    subparsers = parser.add_subparsers(title='options', dest='command')

    # run command
    run_parser = subparsers.add_parser('run', help='Command to run a model with arguments -t | -s | -hy')
    run_parser.add_argument('-t', '--temperature', action='store_true', help='Runs a temperature model.')
    run_parser.add_argument('-s', '--solar', action='store_true', help='Runs solar routines only.')
    run_parser.add_argument('-hy', '--hydraulics', action='store_true', help='Runs hydraulics only.')

    # setup command
    setup_parser = subparsers.add_parser('setup', help='Command to setup a model with arguments -cf | -mi')
    setup_parser.add_argument('-cf', '--control-file', action='store_true', help='Writes a blank control file.')
    setup_parser.add_argument('-mi', '--model-inputs', action='store_true',
                              help='Write blank input files. Control file must already be parameterized.')
    setup_parser.add_argument('-csv', '--csv-mode', action='store_true',
                              help='Use -csv to write a csv (Unicode UTF-8) formatted control file. If this option is '
                                   'not included an Excel file (.xlsx) is written.')
    setup_parser.add_argument('-t', '--timestamp', action='store_true',
                              help='Use -t to add a timestamp to the file name.')
    setup_parser.add_argument('-o', '--overwrite', action='store_true',
                              help='Use -o to overwrite any existing file.')

    arg = parser.parse_args()

    if arg.command == 'run':

        full_path = glob(os.path.join(arg.model_dir, "HeatSource_Control.*"))

        # checks to make sure HeatSource_Control.xlsx or HeatSource_Control.csv exists
        if len(full_path) == 0:
            raise Exception("HeatSource_Control file not found in {0}.".format(arg.model_dir))

        if len(full_path) > 1:
            raise Exception("There is more than one file named 'HeatSource_Control.' \
        in this directory: {0}. Only one file can exist.".format(arg.model_dir))

        control_file = os.path.split(full_path[0])[1]
        control_ext = os.path.splitext(control_file)[1]

        if control_ext not in [".xlsx", ".csv"]:
            raise Exception("{0} must be an Excel '.xlsx' or '.csv' file.".format(control_file))
            
        if arg.temperature:
            run(arg.model_dir, control_file)

        elif arg.solar:
            run_solar(arg.model_dir, control_file)
                
        elif arg.run.hydraulics:
            run_hydraulics(arg.model_dir, control_file)

    if arg.command == 'setup':

        if arg.csv_mode:
            control_file = 'HeatSource_Control.csv'
        else:
            control_file = 'HeatSource_Control.xlsx'

        if arg.control_file:
            # Write a blank control file
            setup_cf(arg.model_dir, control_file, use_timestamp=arg.timestamp, overwrite=arg.overwrite)

        if arg.model_inputs:
            # Write blank input files,
            setup_mi(arg.model_dir, control_file,
                     use_timestamp=arg.timestamp, overwrite=arg.overwrite)


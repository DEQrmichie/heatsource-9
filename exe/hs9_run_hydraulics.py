"""This script imports the heatsource module and executes the 
hydraulic model routines. This script must be located in the same
directory as HeatSource_control.csv. NOTE that executing this script
from Python shell (IDLE) will not identify __file__ correctly and will
result in an error. It must be executed from a command prompt. Your
options are to try to double click on this file and execute it 
using python.exe, or to open a command prompt and execute manually 
by typing: python -i path/to/this/script/hs9_run_hydraulics.py 

Command line:
> hs run -hy

usage: hs <command> [options]

optional arguments:
  -h, --help         show this help message and exit
  -t, --temperature  Runs a temperature model.
  -s, --solar        Runs solar routines only.
  -hy, --hydraulics  Runs hydraulics only.

"""

from heatsource9 import BigRedButton
from os.path import abspath
from os.path import dirname
from os.path import exists
from os.path import join
from os.path import realpath


def getScriptPath():
    """Gets the path to the directory where the script is being executed from."""
    return abspath(join(dirname(realpath(__file__)), '.'))


model_dir = join(getScriptPath(), '')
control_file = 'HeatSource_Control.csv'

if not exists(join(model_dir, control_file)):
    raise Exception("HeatSource_Control.csv not found. \
    Move the executable or place the control file in \
    this directory: {0}.".format(model_dir))

# Run Heat Source Hydraulics only, run_type = 2
BigRedButton.run_hydraulics(model_dir, control_file)

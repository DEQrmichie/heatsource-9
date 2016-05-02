"""DO NOT DELETE THIS FILE. This script imports the heatsource module
and executes the model routines. This script must be located in the same
directory as HeatSource_control.csv. NOTE that executing this script
from Python shell (IDLE) will not identify __file__ correctly and will
result in an error. It must be executed from a command prompt. Your
options are to double click on this file and execute using python.exe,
use the batch command files (which point to these files), or
open terminal and execture manually by typing:
python -i path/to/this/script/HS9_Run_Temperature.py"""

from heatsource9 import BigRedButton
from os.path import abspath
from os.path import dirname
from os.path import exists
from os.path import join
from os.path import realpath

def getScriptPath():
    """Gets the path to the directory where the script is being executed from."""
    return abspath(join(dirname(realpath(__file__)), '.'))

model_dir = getScriptPath() + '/'
control_file = 'HeatSource_Control.csv'

if not exists(join(model_dir,control_file)):
    raise Exception("HeatSource_Control.csv not found. \
    Move the executable or place the control file in \
    this directory: {0}.".format(model_dir))

# Run Heat Source Temperature, run_type = 0
BigRedButton.RunHS(model_dir, control_file)

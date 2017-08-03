"""This script writes a control file with empty values. 
This script must be located in the same directory 
as HeatSource_control.csv. NOTE that executing this script
from Python shell (IDLE) will not identify __file__ correctly and will
result in an error. It must be executed from a command prompt. Your
options are to try to double click on this file and execute it 
using python.exe, or to open a command prompt and execute manually 
by typing: python -i path/to/this/script/HS9_Setup_Control_File.py"""

from heatsource9.ModelSetup.Inputs import Inputs
from os.path import abspath
from os.path import dirname
from os.path import join
from os.path import realpath

def getScriptPath():
    """Gets the path to the directory where the script is being executed from."""
    return abspath(join(dirname(realpath(__file__)), '.'))

model_dir = getScriptPath() + '/'
control_file = 'HeatSource_Control.csv'

# create an input object
inputs = Inputs(model_dir, control_file)

# Write a blank control file
inputs.parameterize_cf(overwrite=False)



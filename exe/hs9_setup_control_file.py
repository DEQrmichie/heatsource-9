"""This script writes a control file with empty values. 
This script must be located in the same directory 
as HeatSource_control.csv. NOTE that executing this script
from Python shell (IDLE) will not identify __file__ correctly and will
result in an error. It must be executed from a command prompt. Your
options are to try to double click on this file and execute it 
using python.exe, or to open a command prompt and execute manually 
by typing: python -i path/to/this/script/hs9_setup_control_file.py

Command line:
> hs setup -cf

usage: hs <command> [options]

optional arguments:
  -h, --help           show this help message and exit
  -cf, --control-file  Writes a blank control file.
  -mi, --model-inputs  Write blank input files. Control file must already be
                       parameterized.
  -t, --timestamp      Use -t to add a timestamp to the file name.
  -o, --overwrite      Use -o to overwrite any existing file.
  
"""

from heatsource9 import BigRedButton
from os.path import abspath
from os.path import dirname
from os.path import join
from os.path import realpath


def getScriptPath():
    """Gets the path to the directory where the script is being executed from."""
    return abspath(join(dirname(realpath(__file__)), '.'))


model_dir = join(getScriptPath(), '')
control_file = 'HeatSource_Control.csv'

BigRedButton.setup_cf(model_dir, control_file,
                      use_timestamp=False, overwrite=False)

#!/usr/bin/python3

"""This script is used for model setup and writes empty input files
specific to the parametrization settings in the control file.
It is used by pyinstaller to build the
standalone executable. This script can also be used
directly. It must be located in the same directory as
HeatSource_control.csv. Run by opening a windows command
prompt and type:

cd path/to/this/script/
py -m hs9_setup_model_inputs

Or create a .bat file that contains the heat source command:
hs setup -mi

Command line:
> hs setup -mi

usage: hs <command> [options]

optional arguments:
  -h, --help           show this help message and exit
  -cf, --control-file  Writes a blank control file.
  -mi, --model-inputs  Write blank input files. Control file must already be
                       parameterized.
  -t, --timestamp      Use -t to add a timestamp to the file name.
  -o, --overwrite      Use -o to overwrite any existing file.
  -csv, --csv-mode     Use -csv to write a csv (Unicode UTF-8) formatted control file instead of .xlsx. Default is .xlsx.

"""

from heatsource9 import BigRedButton
import sys
from os.path import abspath
from os.path import dirname
from os.path import join
from os.path import realpath
from os.path import split
from os.path import splitext
from glob import glob

if getattr(sys, 'frozen', False):
    # path to the directory where the exe is being executed from
    application_path = dirname(sys.executable)
else:
    # path to the directory where the script is being executed from
    application_path = abspath(join(dirname(realpath(__file__)), '.'))

model_dir = join(application_path, '')
full_path = glob(join(model_dir, "HeatSource_Control.*"))

# checks to make sure HeatSource_Control.xlsx or HeatSource_Control.csv exists
if len(full_path) == 0:
    raise Exception("HeatSource_Control file not found. \
        Move the executable or place the control file in \
        this directory: {0}.".format(model_dir))

if len(full_path) > 1:
    raise Exception("There is more than one file named 'HeatSource_Control.' \
        in this directory: {0}. \
        Only one file can exist.".format(model_dir))

control_file = split(full_path[0])[1]
control_ext = splitext(control_file)[1]

if control_ext not in [".xlsx", ".csv"]:
    raise Exception("{0} must be an Excel '.xlsx' or '.csv' file.".format(control_file))

# Write blank input files,
BigRedButton.setup_mi(model_dir, control_file,
                      use_timestamp=True, overwrite=False)
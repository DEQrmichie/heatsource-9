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

"""

from heatsource9 import BigRedButton
import sys
from os.path import abspath
from os.path import dirname
from os.path import join
from os.path import realpath

if getattr(sys, 'frozen', False):
    # path to the directory where the exe is being executed from
    application_path = dirname(sys.executable)
else:
    # path to the directory where the script is being executed from
    application_path = abspath(join(dirname(realpath(__file__)), '.'))

model_dir = join(application_path, '')
control_file = 'HeatSource_Control.csv'

if not exists(join(model_dir, control_file)):
    raise Exception("HeatSource_Control.csv not found. \
    Move the executable or place the control file in \
    this directory: {0}.".format(model_dir))

# Write blank input files,
BigRedButton.setup_mi(model_dir, control_file,
                      use_timestamp=True, overwrite=False)
#!/usr/bin/python3

"""This script imports the heatsource module and executes the
hydraulic model routines. It is used by pyinstaller to build the
standalone executable. This script can also be used
directly. It must be located in the same directory as
HeatSource_control.csv. Run by opening a windows command
prompt and type:

cd path/to/this/script/
py -m hs9_run_hydraulics

Or create a .bat file that contains the heat source command:
hs run -hy

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
import sys
from os.path import abspath
from os.path import dirname
from os.path import exists
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

# Run Heat Source Hydraulics only, run_type = 2
BigRedButton.run_hydraulics(model_dir, control_file)
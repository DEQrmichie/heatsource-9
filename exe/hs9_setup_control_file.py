#!/usr/bin/python3

"""This script writes a control file with empty values.
It is used by pyinstaller to build the
standalone executable. This script can also be used
directly. It must be located in the same directory as
HeatSource_control.csv. Run by opening a windows command
prompt and type:

cd path/to/this/script/
py -m hs9_setup_model_inputs

Or create a .bat file that contains the heat source command:
hs setup -cf

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
  -csv, --csv-mode     Use -c to write a csv (Unicode UTF-8) formatted control instead of .xlsx. Default is .xlsx.
  
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

BigRedButton.setup_cf(model_dir, control_file,
                      use_timestamp=False, overwrite=False)
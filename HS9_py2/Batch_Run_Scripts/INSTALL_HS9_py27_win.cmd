#!/bin/sh

C:

echo "This script will install Heat Source"
echo ""

::#PIP INSTALL
::pip install heatsource-9

::#GIT REMOTE INSTALL
::#Install from git (need to have git installed)
::#This gets installed in user home directory
::git clone git://github.com/rmichie/heatsource-9.git

::#rename the directory
::mv heatsource-9 heatsource9
::cd heatsource9

::#Move the setup.py script out of the
::#heat source folder up one directory
::find . -name "*setup*" -print
::mv "setup.py" ~
::cd

::python setup.py install clean --all
::rm -rf heatsource9
::rm setup.py

::#GIT LOCAL INSTALL
cd ..
cd Workspace/Github/heatsource-9/HS9_py2
python setup.py install clean --all



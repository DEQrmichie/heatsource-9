#!/bin/sh

echo "This script will install Heat Source for python 3"
echo ""

pip3 install pandas
pip3 install xlrd

# Install from pypi USING pip (need to have pip installed)
#pip3 install heatsource9

# Install from git (need to have git installed)
# This gets installed in user home directory
git clone git://github.com/rmichie/heatsource-9.git

# rename the directory
#mv heatsource-9 heatsource9
cd heatsource-9/HS9_py3

python3 setup.py install clean --all
cd
rm -rf heatsource-9



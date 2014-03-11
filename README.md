heatsource-9
============

Current Version: heatsource 9.0.0b5 (beta 5 build 339)

================================================================================================
ABOUT 

Heat Source is a computer model used by the Oregon Department of Environmental Quality to simulate
stream thermodynamics and hydraulic routing. It was originally developed by Matt Boyd in 1996 as a 
Masters Thesis at Oregon State University in the Departments of Bioresource Engineering and Civil Engineering. 
Since then it has grown and changed significantly. Oregon DEQ currently maintains the Heat Source methodology and computer programming. Appropriate model use and application are the sole responsibility of the user.

http://www.deq.state.or.us/wq/TMDLs/tools.htm

Authors: Matt Boyd, Brian Kasper, John Metta, Ryan Michie, Dan Turner

Contact: Ryan Michie, michie.ryan@deq.state.or.us

================================================================================================
INSTALL

See OS specific instructions for step by step install.

Heat Source requires python 2.7 and the following python packages:

-numpy 1.8
-Bottleneck
-numexpr
-python-dateutil 1.5
-pytz
-xlrd
-pandas

================================================================================================
LICENSE

GNU General Public License v3 (GPLv3)

Heat Source, Copyright (C) 2000-2014, Oregon Department of Environmental Quality

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

================================================================================================
UPDATES

Major changes already implemented from last version (8.0.8)

- Heatsource package updated for use in python 2.7.6 with some design considerations made for eventual move to python 3.
- Heatsource package is compatible with Mac OS 10.6 and above.
- Removed excel interface for inputs and model control.
- ExcelDocument module is gone.
- All model inputs are formatted as csv files including the model control parameters. The inputs are similar to the format of  
heatsource 8 with some minor changes.
- csv files are read using the pandas package. Pandas was chosen over the standard python read methods because read/write speed is  
much faster, pandas is platform independent, and it contains desired methods which may be implemented in the future such as reading/writing excel (with xlrd) and other data manipulation and analysis routines.
- The pandas csv reading methods are implemented in the new CSVInterface module which includes most of the methods from the old  
ExcelInterface.
- Removed psyco optimization. Psyco is no longer maintained and not supported for python 2.7.
- C module is gone although it may be reintroduced at a future date with updates.
- Model setup, hydraulics, solar, and temperature routines are initiated via executable apps in Mac OSX (written in applescript). In  
windows they are executed via python scripts. Note the windows scripts are basically hacks and it would be nice to make them  
true executables.
- Reworked __version__ to conform to PEP 386 and added documentation.
- Heat source is now licensed under the GNU General Public License v3.
- Outputs now get stamped with model version, simulation name, and user text in addition to the output parameter, units and date/time.
- Implements variable number of radial samples (use 999 in control file for version 8 method).

================================================================================================
ROADMAP

Roadmap for this version
- Make the windows run scripts into executables.
- Develop routines to convert heatsource v8 inputs to v9 csv format using pandas/xlrd package.
- Switch outputs from txt to csv and write using pandas instead.
- Implement a stop button possibly using methods from pygame or pykeylogger packages.
- Fix/look into the Krieter bug.
- Reformulate Muskingum calcs using updated methods such as Moramarco et al 2006.
- Update the user manual.
- Check elevation is a number.
- General improvement to the QA/QC of certain inputs (see commented lines in CSVinterface).
- More user control of light extinction coefficients/density parameter, plus integrate more recent formulations for beer's law. Proposal to use LAI/PAI (simple change).
- Register heatsource with pypi so it can be automatically downloaded and installed via python package managers (pip/easy_install).
- Make a heatsource package mpkg installer for mac.
- Decide if cloudiness is a required input for solar runs.


Future Roadmap
- Allow variable timeseries input and utilize pandas interpolation methods.
- Fix the bug where the model ignores tributaries and accretion flows at first node.
- Look into issues with including evaporation losses (much higher rates in the first reach than rest of the model).
- Review cloudiness routines.
- Output hyporheic energy flux.
- Consider methodology change for hyporheic to consider recent publications.
- Output longitudinal landcover data (i.e. like vegematic in heatsource 7).
- User control over bottom reflection.
- Input option for solar radiation measurements.
- Implement a conservative tracer.
- Provide flexibility on longwave routines.
- Fix/implement the logger.
- Make a pause button.
- Post processing tools (graphing, analysis, etc). Some already written in R. Consider pandas and matplotlib for this as part of heatsource.
- Better output messages on errors.
- Allow user specified output dT or sub-hourly time step.
- Simplify the CSV interface class by separating read methods from node building methods.
- Implement seamless transition to different watershed/stream network scales.
- Scope implementing solar routines as a separate independent package or something less complex so it can be accessed via GIS routines.
- Consider using long format for landcover data as opposed to the current wide format (for GIS based solar modeling).
- Consider moving stream elevation into landcover data (so all GIS sampled data goes in the same location).

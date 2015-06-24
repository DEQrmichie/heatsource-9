Heat Source 9
============

Current Version: heatsource 9.0.0b13 (beta 13)

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

Heat Source 9 works on Windows, Mac, and Linux.

Requires python 2.7
https://www.python.org/downloads/

Download the zip file
```git clone git://github.com/rmichie/heatsource-9.git```

navigate to setup.py and install
```python setup.py install clean --all```

================================================================================================
RUNNING THE MODEL

1. Put the control file (HeatSource_Control.csv) and the model run scripts 
   in the same directory. You can generate a template control file by double 
   clicking on HS9_Setup_Control_File.py.
2. Open the control file and parameterize it with your model information.
   The control file must be named HeatSource_Control.csv
3. Use HS9_Setup_Model_Inputs.py to build template input files (they will 
   be saved to your input file directory that is specified in the control file).
5. Edit the template csv files with your input data. You can use excel. 
   Save it as a csv.
6. Run the model by double clicking one of the following model run scripts:
   HS9_Run_Hydraulics_Only.py,
   HS9_Run_Solar_Only.py,
   HS9_Run_Temperature.py,
   A console should open and you should see the model running.
7. Outputs are saved in the output directory (specified in the control file).

================================================================================================
INPUT FILES - GENERAL INFORMATION

1. The heat source control file must be named HeatSource_Control.csv
2. The other input csv files can be named whatever you want (file names are specified in the control file).
3. Do not change the parameter names in the control file. Only change the VALUE column (column 3).
4. The column header names can be changed but the data needs to be in 
   the correct column number (see below).
5. Use the specified unit and unit format identified in the control files 
   and input files. Example mm/dd/yyyy hh:mm is 07/01/2001 16:00
6. A parameter value that is not applicable may be left blank.
7. Tributary and climate data can have multiple files (or just one file).
   Separate file names and site stream km with commas (see below).
   
================================================================================================
CONTROL FILE
HeatSource_Control.csv

Below are all the input parameters that must be included in the control file.
Technically the parameters can be in any row as long as the names in column 2 don't change and the parameter values are in column 3.
If "USE HEAT SOURCE 8 LANDCOVER METHODS" is TRUE, the model will revert back to heat source 8 methods for landcover sampling methods.
If FALSE, the direction for each the radial samples (n) measured in degrees from 0 (0 = north) will correspond to (total radial samples / 360) x radial sample n. 
In addition the first landcover sample in each direction will be at the stream node.

 
|ROW COLUMN| [1] | [2]                                               | [3]  |

|[1]       |LINE |PARAMETER                                          |VALUE |
|:---------|----:|:--------------------------------------------------|----- |
|[2]       |1    |USER TEXT |      |
|[3]       |2    |SIMULATION NAME |      |
|[4]       |3    |STREAM LENGTH (KILOMETERS) |      |
|[5]       |4    |OUTPUT PATH |      |
|[6]       |5    |INPUT PATH |      |
|[7]       |6    |DATA START DATE (mm/dd/yyyy) |      |
|[8]       |7    |MODELING START DATE (mm/dd/yyyy) |      |
|[9]       |8    |MODELING END DATE (mm/dd/yyyy) |      |
|[10]      |9   |DATA END DATE (mm/dd/yyyy) |      |
|[11]      |10  |FLUSH INITIAL CONDITION (DAYS) |      |
|[12]      |11  |TIME OFFSET FROM UTC (HOURS) |      |
|[13]      |12  |MODEL TIME STEP - DT (MIN) |      |
|[14]      |13  |MODEL DISTANCE STEP - DX (METERS) |      |
|[15]      |14  |LONGITUDINAL STREAM SAMPLE DISTANCE (METERS) |      |
|[16]      |15  |BOUNDARY CONDITION FILE NAME |      |
|[17]      |16  |TRIBUTARY SITES |      |
|[18]      |17  |TRIBUTARY INPUT FILE NAMES |      |
|[19]      |18  |TRIBUTARY MODEL KM |      |
|[20]      |19  |ACCRETION INPUT FILE NAME |      |
|[21]      |20  |CLIMATE DATA SITES |      |
|[22]      |21  |CLIMATE INPUT FILE NAMES |      |
|[23]      |22  |CLIMATE MODEL KM |      |
|[24]      |23  |INCLUDE EVAPORATION LOSSES FROM FLOW (TRUE/FALSE) |      |
|[25]      |24  |EVAPORATION METHOD (Mass Transfer/Penman) |      |
|[26]      |25  |WIND FUNCTION COEFFICIENT A |      |
|[27]      |26  |WIND FUNCTION COEFFICIENT B |      |
|[28]      |27  |INCLUDE DEEP ALLUVIUM TEMPERATURE (TRUE/FALSE) |      |
|[29]      |28  |DEEP ALLUVIUM TEMPERATURE (*C) |      |
|[30]      |29  |MORPHOLOGY DATA FILE NAME |      |
|[31]      |30  |LANDCOVER DATA FILE NAME |      |
|[32]      |31  |LANDCOVER CODES FILE NAME |      |
|[33]      |32  |NUMBER OF TRANSECTS PER NODE |      |
|[34]      |33  |NUMBER OF SAMPLES PER TRANSECT |      |
|[35]      |34  |DISTANCE BETWEEN TRANSESCT SAMPLES (METERS) |      |
|[36]      |35  |ACCOUNT FOR EMERGENT VEG SHADING (TRUE/FALSE) |      |
|[37]      |36  |LANDCOVER DATA INPUT TYPE (Codes/Values) |      |
|[38]      |37  |CANOPY DATA TYPE (LAI/CanopyCover) |      |
|[39]      |38  |VEGETATION ANGLE CALCULATION METHOD (point/zone) |      |
|[40]      |39  |USE HEAT SOURCE 8 LANDCOVER METHODS (TRUE/FALSE) |      |

================================================================================================
ACCRETION INPUT FILE
UserDefinedFileName.csv

The temperature and flow rates of accretion are defined in this file. 

Accretion flows are inflows that enter the stream over more than one 
stream data node, and typically are subsurface seeps that occur over 
longer distances than discrete subsurface inflows (i.e. a spring).  
 
When accretion flows are close enough so that more than one occurs 
in a model distance step, the accretion flow rates will be summed and a 
flow based average accretion temperature will be derived and used 
in the mixing calculations

|ROW/COLUMN| [1]      | [2]     | [3]      | [4]    | [5]         | [6]     |

|[1]       |STREAM_ID | NODE_ID |STREAM_KM | INFLOW | TEMPERATURE | OUTFLOW |
|:---------|----------|---------|----------|--------|-------------|---------|
|[2]       | Value    | Value   | Value    | Value  | Value       | Value   |

COLUMN
[1] STREAM_ID (Optional)
[2] NODE_ID (Optional)
[3] STREAM_KM (headwaters at the top)
[4] INFLOW (Accretion flow cms)
[5] TEMPERATURE (Accretion flow temperature Celsius)
[6] OUTFLOW (withdrawal flows cms)

================================================================================================
BOUNDARY CONDITION FILE
UserDefinedFileName.csv

The stream flow and temperature condtions at the upstream model boundary are defined in this file.
The boundary condtions are defined at an hourly timestep.

|ROW/COLUMN| [1]      | [2]     | [3]       |

|[1]       | DATETIME | FLOW    |TEMPERATURE|
|:---------|----------|---------|-----------|
|[2]       | Value    | Value   | Value     |

COLUMN
[1] DATETIME (MM/DD/YYYY HH:MM)
[2] FLOW (cms)
[3] TEMPERATURE (Celsius)

================================================================================================
CLIMATE INPUT FILE/S 
(formally called Continuous data in heat source 8)
UserDefinedFileName.csv


|ROW/COLUMN| [1]      | [2]        | [3]       | [4]               | [5]             |

|[1]       | DATETIME | CLOUDINESS |WIND_SPEED | RELATIVE_HUMIDITY | AIR_TEMPERATURE |
|:---------|----------|------------|-----------|-------------------|-----------------|
|[2]       | Value    | Value      | Value     | Value             | Value           |

COLUMN
[1] DATETIME (MM/DD/YYYY HH:MM)
[2] CLOUDINESS (0-1)
[3] WIND_SPEED (meters/second)
[4] RELATIVE_HUMIDITY (0-1)
[5] AIR_TEMPERATURE (Celsius)

Note - multiple csv files may be used for each set of climate inputs with the 
format above or all data can be saved in the same file like the example below. 
This is controlled in the control file by designating the number of inputs and 
the input stream km.

COLUMN
[1] DATETIME (MM/DD/YYYY HH:MM)
[2] CLOUDINESS 1 (0-1)
[3] WIND_SPEED 1 (meters/second)
[4] RELATIVE_HUMIDITY 1 (0-1)
[5] AIR_TEMPERATURE 1 (Celsius)
[6] CLOUDINESS 2 (0-1)
[7] WIND_SPEED 2 (meters/second)
[8] RELATIVE_HUMIDITY 2 (0-1)
[9] AIR_TEMPERATURE 2 (Celsius)
[n] CLOUDINESS n (0-1)
[n] WIND_SPEED 2 (meters/second)
[n] RELATIVE_HUMIDITY n (0-1)
[N] AIR_TEMPERATURE n (Celsius)

================================================================================================
TRIBUTARY INPUT FILE/S
Can also be outflows. Use negative flows.
UserDefinedFileName.csv

The tributary input files define the inflow/outflow rates and temperatures
at different points along the model stream. Inflows refers to localized 
(non-accretion) type flows such as tributaries, springs, returns, point sources, etc.
Outflows can be varous types of water withdrawals. They are input with a negative flow rate.
Temperatures for outflows are not used by the model.

The the number and stream km of the inflow/outflows is defined in the control file.
The flow and temperature are defined at an hourly timestep.  

|ROW/COLUMN| [1]      | [2]     | [3]       |

|[1]       | DATETIME | FLOW    |TEMPERATURE|
|----------|----------|---------|-----------|
|[2]       | Value    | Value   | Value     |

COLUMN
[1] DATETIME (MM/DD/YYYY HH:MM)
[2] FLOW (inflow/outflow cms)
[3] TEMPERATURE (Celsius)

Note - multiple csv files may be created for each tributary input with the
format above or all data can be saved in the same file like in the example below.
This is controlled in the control file by designating the number of inputs and
the input stream km.

COLUMN
[1] DATETIME (MM/DD/YYYY HH:MM)
[2] FLOW 1 (inflow/outflow cms)
[3] TEMPERATURE 1 (Celsius)
[4] FLOW 2 (inflow/outflow cms)
[5] TEMPERATURE 2 (Celsius)
[n] FLOW n (inflow/outflow cms)
[N] TEMPERATURE n (Celsius)

================================================================================================
Landcover Codes
UserDefinedFileName.csv

The landcover codes file contains the physical attribute information 
associated with each land cover code. Land cover codes can be alpha 
numeric values. zero should be avoided as a land cover code. The physical 
attribute such as height, canopy cover, LAI or overhang,  must be numeric.
There cannot be skipped rows 
(i.e. rows without information in between rows with information) 
because the model routines see a blank row as the end of the data sequence.

Canopy Type

land cover canopy information can be input as either canopy cover or 
effective leaf area index. This option is specficed in the control file.

Input file formatting when using canopy cover.

|ROW/COLUMN| [1]      | [2]     | [3]      | [4]          | [5]         |

|[1]       | NAME     | CODE    |HEIGHT    | CANOPY_COVER | OVERHANG    |
|----------|----------|---------|----------|--------------|-------------|
|[2]       | Value    | Value   | Value    | Value        | Value       |

COLUMN
[1] NAME (Landcover Name)
[2] CODE (alpha or numeric code)
[3] HEIGHT (meters)
[4] CANOPY_COVER (0-1)
[5] OVERHANG (meters)

Input file formatting when using LAI.

|ROW/COLUMN| [1]      | [2]     | [3]      | [4]          | [5]         | [6]      |

|[1]       | NAME     | CODE    |HEIGHT    | LAI          | k           | OVERHANG |
|----------|----------|---------|----------|--------------|-------------|----------|
|[2]       | Value    | Value   | Value    | Value        | Value       | Value    |


COLUMN
[1] NAME (Landcover Name)
[2] CODE (alpha or numeric code)
[3] HEIGHT (meters)
[4] LAI (Effective Leaf Area Index)
[5] k Extinction coefficient  (dimensionless)
[6] OVERHANG (meters)

================================================================================================
LANDCOVER DATA
(formally called TTools in heatsource 8)
UserDefinedFileName.csv

This file defines land cover information. This data can be derived 
from geospatial data using TTools.

|ROW/COLUMN| [1]      | [2]     | [3]      | [4]       | [5]      | [6]    | [7]    | [8]    | [9] -[N] |

|[1]       |STREAM_ID | NODE_ID |STREAM_KM | LONGITUDE | LATITUDE | TOPO_W | TOPO_S | TOPO_E | ....     |
|----------|----------|---------|----------|-----------|----------|--------|--------|--------|----------|
|[2]       | Value    | Value   | Value    | Value     | Value    | Value  | Value  | Value  | ....     |

COLUMN
[1] STREAM_ID (Optional)
[2] NODE_ID (Optional)
[3] STREAM_KM  (headwaters at the top)
[4] LONGITUDE (decimal degrees)
[5] LATITUDE (decimal degrees)
[6] TOPO_W (degrees)
[7] TOPO_S (degrees)
[8] TOPO_E (degrees)
[9-n] Landcover Samples (code or height meters)
[n-n] Elevation Samples (meters)
[n-n] Canopy/LAI Samples (blank if not used)
[n-N] k Extinction Coefficient (dimensionless) blank if not used.
Note - the column numbers for the landcover, elevation, and caonpy cover/LAI samples are dependent on user specified information in the control file.

================================================================================================
MORPHOLOGY DATA FILE
UserDefinedFileName.csv

This file defines channel morphology and substrate information.
Refer to the user manual for more information about each parameter.

COLUMN
[1] Stream ID (Optional)
[2] Node ID (Optional)
[3] STREAM_KM (headwaters at the top)
[4] ELEVATION (meters)
[5] GRADIENT (meters/meters)
[6] BOTTOM_WIDTH (meters)
[7] CHANNEL_ANGLE_Z
[8] MANNINGS_n (unitless)
[8] SED_THERMAL_CONDUCTIVITY (Watts/meter/Celsius)
[10] SED_THERMAL_DIFFUSIVITY (square cm/second)
[11] SED_HYPORHEIC_THICKNESSS - Sediment/hyporheic zone thickness (meters)
[12] %HYPORHEIC_EXCHANGE - Percent Hyporheic exchange (0-1)
[13] POROSITY (0-1)

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

- Updated for use in python 2.7
- Now compatible with Linux and Mac OS 10.6 and above.
- Removed excel interface along with ExcelDocument module.
- All model inputs and outputs are formatted as csv files including the model control parameters. The inputs are similar to the format of  
heatsource 8 with some minor changes.
- Model setup, hydraulics, solar, and temperature routines are initiated via python exectuable scripts.
- Removed psyco optimization. Psyco is no longer maintained or supported for python 2.7.
- C module was removed although it may be reintroduced at a future date with updates.
- Reworked __version__ to conform to PEP 386 and added documentation.
- Heat source is now licensed under the GNU General Public License v3.
- Outputs files get stamped with model version, simulation name, user text in addition to the output parameter, units and date/time.
- Implements variable number of radial samples (use heatsource 8 switch in control file for version 8 method).
- Added effective LAI option
- Added solar azimuth and solar elevation outputs

================================================================================================
ROADMAP

Roadmap for this version
- Develop routines to convert heatsource v8 inputs to v9 csv format using pandas/xlrd package.
- Fix/look into the Krieter bug.
- Reformulate Muskingum calcs using updated methods such as Moramarco et al 2006.
- Update the user manual.
- Check elevation is a number.
- General improvement to the QA/QC of certain inputs (see commented lines in CSVinterface).
- Register heat source with pypi
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
- Make a pause button.
- Post processing tools (graphing, analysis, etc). Some already written in R. Consider pandas and matplotlib for this as part of heatsource.
- Better output messages on errors.
- Allow user specified output dT or sub-hourly time step.
- Simplify the ModelSetup class by separating read methods from node building methods.
- Implement seamless transition to different watershed/stream network scales.
- Scope implementing solar routines as a separate independent package or something less complex so it can be accessed via GIS routines.
- Consider using long format for landcover data as opposed to the current wide format (for GIS based solar modeling).
- Consider moving stream elevation into landcover data (so all GIS sampled data goes in the same location).

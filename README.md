Heat Source 9
============

Current Version: heatsource 9.0.0b20 (beta 20)

================================================================================================
## ABOUT 

Heat Source is a computer model used by the Oregon Department of 
Environmental Quality to simulate stream thermodynamics and hydraulic 
routing. It was originally developed by Matt Boyd in 1996 as a [Masters 
Thesis][1] at Oregon State University in the Departments of Bioresource 
Engineering and Civil Engineering. Since then it has grown and changed 
significantly. Oregon DEQ currently maintains the Heat Source methodology
and computer programming. Appropriate model use and application are 
the sole responsibility of the user.

Heat Source 8 and user manual:
http://www.oregon.gov/deq/wq/tmdls/Pages/TMDLs-Tools.aspx

Authors: Matt Boyd, Brian Kasper, John Metta, Ryan Michie, Dan Turner

Contact: Ryan Michie, michie.ryan@deq.state.or.us

[1]: http://ir.library.oregonstate.edu/xmlui/handle/1957/27036

================================================================================================
## INSTALL

Heat Source 9 works on Windows, Mac, and Linux.

Requires python 2.7
https://www.python.org/downloads/

Download  
```git clone git://github.com/rmichie/heatsource-9.git```

navigate to setup.py and install  
```python setup.py install```

================================================================================================
## RUNNING THE MODEL

1. Place the control file (HeatSource_Control.csv) and the model run
   scripts in the same directory. You can generate a template control 
   file by executing HS9_Setup_Control_File.py.
2. Open the control file and parameterize it with your model information. 
   The control file must be named HeatSource_Control.csv
3. Use HS9_Setup_Model_Inputs.py to build template input files (they will
   be saved to your input file directory that is specified in the control file).
5. Edit the template csv files with your input data. You can use excel. 
 Save it as a csv.
6. Run the model by executing one of the following model run scripts:
   HS9_Run_Hydraulics_Only.py,
   HS9_Run_Solar_Only.py,
   HS9_Run_Temperature.py,
   A console should open and you should see the model running.
7. Outputs are saved in the output directory (specified in the control file).

================================================================================================
## INPUT FILES - GENERAL INFORMATION

1. Control and input files are ASCII comma delimited files.
2. The heat source control file must be named HeatSource_Control.csv
3. The other input files can be named whatever you want 
 (file names are specified in the control file).
4. Do not change the parameter names in the control file. Only change 
   the VALUE column (column 4).
5. The column header names can be changed but the data needs to be in 
 the correct column number (see below).
6. Use the specified unit and data format identified in the control files 
 and input files. Example mm/dd/yyyy hh:mm is 07/01/2001 16:00
7. A parameter value that is not applicable may be left blank although all values with float
	data type will be assigned as zero.
8. Tributary and climate data can have multiple files (or just one file).
 In the control file, separate file names and site stream km with 
 commas wrapped in quotes.

================================================================================================
### CONTROL FILE  
HeatSource_Control.csv

The control file is where most of the model operation and initial parameterization is set.


```
control_file = 'HeatSource_Control.csv'
model_dir = r'C://path/to/model_directory/'

# create an input object
inputs = Inputs(model_dir, control_file)

# Write a blank control file
inputs.parameterize_cf(overwrite=False)
```

Below are all the input parameters that must be included in the control file.

Separate tributary and climate file names and/or site stream km with 
 commas wrapped in quotes.
```"inflow1.csv, inflow2.csv"```


|LINE |PARAMETER                                          |KEY    |VALUE |
|----:|:--------------------------------------------------|-------|----- |
|    2|Model Description/User Notes						  |usertxt|	     |
|	 3|Simulation Name            						  | name  |	 	 |
|	 4|Input Directory Path|inputdir|	|
|	 5|Output Directory Path|outputdir|	|
|	 6|Stream Length (kilometers)|length|	|
|	 7|Output Stream Kilometers|outputkm|	|
|	 8|Data Start Date (mm/dd/yyyy)|datastart|	|
|	 9|Modeling Start Date (mm/dd/yyyy)|modelstart|	|
|	10|Modeling End Date (mm/dd/yyyy)|modelend|	|
|	11|Data End Date (mm/dd/yyyy)|dataend|	|
|	12|Flush Initial Condition (days)|flushdays|	|
|	13|Time Offset From UTC (hours)|offset|	|
|	14|Model Time Step (minutes)|dt|	|
|	15|Model Distance Step (meters)|dx|	|
|	16|Longitudinal Stream Sample Distance (meters)|longsample|	|
|	17|Boundary Condition Input File Name|bcfile|	|
|	18|Tributary Inflow Sites|inflowsites|	|
|	19|Tributary Inflow Input File Name|inflowinfiles|	|
|	20|Tributary Inflow Model kilometers|inflowkm|	|
|	21|Accretion Input File Name|accretionfile|	|
|	22|Meteorological Data Sites|metsites|	|
|	23|Meteorological Data Input File Name|metfiles|	|
|	24|Meteorological Data Model kilometers|metkm|	|
|	25|Include Evaporation Losses From Flow (True/False)|calcevap|	|
|	26|Evaporation Method (Mass Transfer/Penman)|evapmethod|	|
|	27|Wind Function Coefficient a|wind_a|	|
|	28|Wind Function Coefficient b|wind_b|	|
|	29|Include Deep Alluvium Temperature (True/False)|calcalluvium|	|
|	30|Deep Alluvium Temperature (Celsius)|alluviumtemp|	|
|	31|Morphology Input Data File Name|morphfile|	|
|	32|Land Cover Input Data File Name|lcdatafile|	|
|	33|Land Cover Codes Input File Name|lccodefile|	|
|	34|Number Of Transects Per Node|trans_count|	|
|	35|Number Of Samples Per Transect|transsample_count|	|
|	36|Distance Between Transect Samples (meters)|transsample_distance|	|
|	37|Account For Emergent Veg Shading (True/False)|emergent|	|
|	38|Land Cover Data Input Type (Codes/Values)|lcdatainput|	|
|	39|Canopy Data Type (LAI/CanopyCover)|canopy_data|	|
|	40|Land Cover Sample Method (point/zone)|lcsampmethod|	|
|	41|Use Heat Source 8 Land Cover Methods (True/False)|heatsource8|	|

================================================================================================
### ACCRETION INPUT FILE  
UserDefinedFileName.csv

The temperature and flow rates of accretion are defined in this file. 

Accretion flows are inflows that enter the stream over more than one 
stream data node, and typically are subsurface seeps that occur over 
longer distances than discrete subsurface inflows (i.e. a spring).  

When accretion flows are close enough so that more than one occurs 
in a model distance step, the accretion flow rates will be summed and a 
flow based average accretion temperature will be derived and used 
in the mixing calculations

|STREAM_ID | NODE_ID |STREAM_KM   | INFLOW | TEMPERATURE | OUTFLOW |
|----------|---------|------------|--------|-------------|---------|
| Value    | Value   | upstream   | Value  | Value       | Value   |
| ...      | ...     | ...        | ...    | ...         | ...     |
| Value    | Value   | downstream | Value  | Value       | Value   |

COLUMN  
 (1) STREAM_ID (Optional)
 (2) NODE_ID (Optional)
 (3) STREAM_KM (headwaters at the top)  
 (4) INFLOW (Accretion flow cms)  
 (5) TEMPERATURE (Accretion flow temperature Celsius)  
 (6) OUTFLOW (withdrawal flows cms)  

================================================================================================
### BOUNDARY CONDITION FILE  
UserDefinedFileName.csv

The stream flow and temperature conditions at the upstream model boundary 
are defined in this file. The boundary conditions are defined at an 
hourly timestep.

| DATETIME | FLOW    |TEMPERATURE|
|----------|---------|-----------|
| start    | Value   | Value     |
| ...      |...      | ...       |
| end      | Value   | Value     |

COLUMN  
 (1) DATETIME (MM/DD/YYYY HH:MM)  
 (2) FLOW (cms)  
 (3) TEMPERATURE (Celsius)  

================================================================================================
### METEOROLOGICAL INPUT FILE/S   
(formally called Continuous data in heat source 8)
UserDefinedFileName.csv

| DATETIME | CLOUDINESS1 |WIND_SPEED1 | RELATIVE_HUMIDITY1 | AIR_TEMPERATURE1 |
|----------|-------------|------------|--------------------|------------------|
| start    | Value       | Value      | Value              | Value            |
| ...      | ...         | ...        | ...                | ...              |
| end      | Value       | Value      | Value              | Value            |

COLUMN
 (1) DATETIME (MM/DD/YYYY HH:MM)
 (2) CLOUDINESS (0-1)
 (3) WIND_SPEED (meters/second)
 (4) RELATIVE_HUMIDITY (0-1)
 (5) AIR_TEMPERATURE (Celsius)

Note - multiple csv files may be used for each set of meteorological inputs with the 
format above or all data can be saved in the same file like the example below. 
This is controlled in the control file by designating the number of inputs and 
the input stream km.

COLUMN  
 (1) DATETIME (MM/DD/YYYY HH:MM)  
 (2) CLOUDINESS 1 (0-1)  
 (3) WIND_SPEED 1 (meters/second)  
 (4) RELATIVE_HUMIDITY 1 (0-1)  
 (5) AIR_TEMPERATURE 1 (Celsius)  
 (6) CLOUDINESS 2 (0-1)  
 (7) WIND_SPEED 2 (meters/second)  
 (8) RELATIVE_HUMIDITY 2 (0-1)  
 (9) AIR_TEMPERATURE 2 (Celsius)  
 (n) CLOUDINESS n (0-1)  
 (n) WIND_SPEED 2 (meters/second)  
 (n) RELATIVE_HUMIDITY n (0-1)  
 (N) AIR_TEMPERATURE n (Celsius)  

================================================================================================
### TRIBUTARY INPUT FILE/S  
Can also be outflows. Use negative flows.  
UserDefinedFileName.csv  

The tributary input files define the inflow/outflow rates and temperatures
at different points along the model stream. Inflows refers to localized 
(non-accretion) type flows such as tributaries, springs, returns, point 
sources, etc. Outflows can be various types of water withdrawals. They 
are input with a negative flow rate. Temperatures for outflows are not 
used by the model.

The number and stream km of the inflow/outflows is defined in the control file.
The flow and temperature are defined at an hourly timestep.  

| DATETIME | FLOW1   | TEMPERATURE1 |
|----------|---------|--------------|
| start    | Value   | Value        |
| ...      |...      | ...          |
| end      | Value   | Value        |

COLUMN  
 (1) DATETIME (MM/DD/YYYY HH:MM)  
 (2) FLOW (inflow/outflow cms)  
 (3) TEMPERATURE (Celsius)  

Note - multiple csv files may be created for each tributary input with the
format above or all data can be saved in the same file like in the example below.
This is controlled in the control file by designating the number of inputs and
the input stream km.

COLUMN  
 (1) DATETIME (MM/DD/YYYY HH:MM)  
 (2) FLOW 1 (inflow/outflow cms)  
 (3) TEMPERATURE 1 (Celsius)  
 (4) FLOW 2 (inflow/outflow cms)  
 (5) TEMPERATURE 2 (Celsius)  
 (n) FLOW n (inflow/outflow cms)  
 (N) TEMPERATURE n (Celsius)  

================================================================================================
### LANDCOVER CODES FILE  
UserDefinedFileName.csv  

The landcover codes file contains the physical attribute information 
associated with each land cover code. Land cover codes can be alpha 
numeric values. zero should be avoided as a land cover code. The physical 
attribute such as height, canopy cover, LAI or overhang,  must be numeric.
There cannot be skipped rows 
(i.e. rows without information in between rows with information) 
because the model routines see a blank row as the end of the data sequence.

#### Canopy Type

land cover canopy information can be input as either canopy cover or 
effective leaf area index. This option is specified in the control file.

##### Canopy Cover  
Input file formatting when using canopy cover.

| NAME     | CODE    |HEIGHT    | CANOPY_COVER | OVERHANG    |
|----------|---------|----------|--------------|-------------|
| Value    | Value   | Value    | Value        | Value       |

COLUMN  
 (1) NAME (Landcover Name)  
 (2) CODE (alpha or numeric code)  
 (3) HEIGHT (meters)  
 (4) CANOPY_COVER (0-1)  
 (5) OVERHANG (meters)  

##### LAI  
Input file formatting when using LAI.

| NAME     | CODE    |HEIGHT    | LAI          | K           | OVERHANG |
|----------|---------|----------|--------------|-------------|----------|
| Value    | Value   | Value    | Value        | Value       | Value    |


COLUMN  
 (1) NAME (Landcover Name)  
 (2) CODE (alpha or numeric code)  
 (3) HEIGHT (meters)  
 (4) LAI (Effective Leaf Area Index)  
 (5) K Extinction coefficient  (dimensionless)  
 (6) OVERHANG (meters)  

================================================================================================
### LANDCOVER DATA  
(formally called TTools in heatsource 8)  
UserDefinedFileName.csv  

This file defines land cover information. This data can be derived 
from geospatial data using TTools.

|STREAM_ID | NODE_ID |STREAM_KM   | LONGITUDE | LATITUDE | TOPO_W | TOPO_S | TOPO_E | ...   |
|----------|---------|------------|-----------|----------|--------|--------|--------|-------|
| Value    | Value   | upstream   | Value     | Value    | Value  | Value  | Value  | ...   |
| ...      | ...     | ...        | ...       | ...      | ...    | ...    | ...    | ...   |
| Value    | Value   | downstream | Value     | Value    | Value  | Value  | Value  | ...   |

COLUMN  
 (1) STREAM_ID (Optional)  
 (2) NODE_ID (Optional)  
 (3) STREAM_KM  (headwaters at the top)  
 (4) LONGITUDE (decimal degrees)  
 (5) LATITUDE (decimal degrees)  
 (6) TOPO_W (degrees)  
 (7) TOPO_S (degrees)  
 (8) TOPO_E (degrees)  
 (9-n) Landcover Samples (code or height meters)  
 (n-n) Elevation Samples (meters)  
 (n-n) Canopy/LAI Samples (blank if not used)  
 (n-N) K Extinction Coefficient (dimensionless) blank if not used.  

Note - the column numbers for the landcover, elevation, and canopy cover/LAI samples are dependent on user specified information in the control file.

================================================================================================
### MORPHOLOGY DATA FILE  
UserDefinedFileName.csv  

This file defines channel morphology and substrate information.
Refer to the user manual for more information about each parameter.

COLUMN  
 (1) STREAM_ID (Optional)  
 (2) NODE_ID (Optional)  
 (3) STREAM_KM (headwaters at the top)  
 (4) ELEVATION (meters)  
 (5) GRADIENT (meters/meters)  
 (6) BOTTOM_WIDTH (meters)  
 (7) CHANNEL_ANGLE_Z  
 (8) MANNINGS_n (unitless)  
 (8) SED_THERMAL_CONDUCTIVITY (Watts/meter/Celsius)  
 (10) SED_THERMAL_DIFFUSIVITY (square cm/second)  
 (11) SED_HYPORHEIC_THICKNESSS - Sediment/hyporheic zone thickness (meters)  
 (12) %HYPORHEIC_EXCHANGE - Percent Hyporheic exchange (0-1)  
 (13) POROSITY (0-1)  

================================================================================================
## LICENSE

GNU General Public License v3 (GPLv3)

Heat Source, Copyright (C) 2000-2017, Oregon Department of Environmental Quality

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
## ROADMAP

Roadmap for this version
- Develop routines to convert heatsource v8 inputs to v9 csv format.
- Fix/look into the Krieter bug.
- Reformulate Muskingum calcs using updated methods such as Moramarco et al 2006.
- Update the user manual.
- Register heat source with pypi
- Make a heatsource package mpkg installer for mac.

Future Roadmap
- Allow variable timeseries input and utilize pandas interpolation methods.
- Fix the bug where the model ignores tributaries and accretion flows at first node.
- Look into issues with including evaporation losses (much higher rates in the first reach than rest of the model).
- Review cloudiness routines.
- Output hyporheic energy flux.
- Consider methodology change for hyporheic
- Output longitudinal landcover data (i.e. like vegematic in heatsource 7).
- User control over bottom reflection.
- Input option for solar radiation measurements.
- Implement a conservative tracer.
- Provide flexibility on longwave routines.
- Make a pause button.
- Post processing tools (graphing, analysis, etc). Some already written in R. Consider pandas and matplotlib for this as part of heatsource.
- Better output messages on errors.
- Allow user specified output dT or sub-hourly time step.

- Implement seamless transition to different watershed/stream network scales.
- Scope implementing solar routines as a separate independent package or something less complex so it can be accessed via GIS routines.
- Consider using long format for landcover data as opposed to the current wide format (for GIS based solar modeling).
- Consider moving stream elevation into landcover data (so all GIS sampled data goes in the same location).

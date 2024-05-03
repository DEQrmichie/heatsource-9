![Wheel Build](https://github.com/DEQrmichie/heatsource/actions/workflows/build_wheels.yml/badge.svg?branch=main)

Heat Source 9
-------------
Current Version: heatsource 9.0.0b27 (beta 27)

## ABOUT 
Heat Source is a computer model used by the Oregon Department of 
Environmental Quality to simulate stream thermodynamics and hydraulic 
routing. It was originally developed by Matt Boyd in 1996 as a [Masters 
Thesis][1] at Oregon State University in the Departments of Bioresource 
Engineering and Civil Engineering. Since then, it has grown and changed 
significantly. Oregon DEQ currently maintains the Heat Source methodology
and computer programming. Appropriate model use and application are 
the sole responsibility of the user. 

Heat Source 7-8 and user manual: 
http://www.oregon.gov/deq/wq/tmdls/Pages/TMDLs-Tools.aspx

Authors: Matt Boyd, Brian Kasper, John Metta, Ryan Michie, Dan Turner

Contact: Ryan Michie, ryan.michie @ deq.oregon.gov

[1]: http://ir.library.oregonstate.edu/xmlui/handle/1957/27036

## INSTALL

Download the [heat source python wheel][3] appropriate to your OS platform and python version. 
Python wheels have been built to support Windows, Mac, and Linux.
DEQ uses windows so other platforms have limited testing. Users of these platforms 
have reported success. 
Requires install of Python 3.8, 3.9, 3.10, 3.11, or 3.12.

https://www.python.org/downloads/

Install the wheel from command line using pip:
```shell
# These commands are for windows
cd path\to\directory_where_you_saved_the_heatsource9_wheel\
py -m pip install <name of wheel file>

# Installs the Python 3.12 heatsource wheel for windows in the local directory
py -m pip install heatsource9-9.0.0b27-cp312-cp312-win32.whl --user

# Installs the Python 3.12 heatsource wheel for windows in the global directory
py -m pip install heatsource9-9.0.0b27-cp312-cp312-win32.whl


 ```
[3]: https://github.com/DEQrmichie/heatsource-9/releases/tag/v9.0.0b27

## QUICK STEPS TO GET GOING

1. Place the control file (HeatSource_Control.csv) and the model run
   scripts in the same directory. You can generate a template control 
   file by executing *hs_setup_control_file.py* or by using commend line.
   ```shell
   cd path\to\model_directory
   hs setup -cf
   ```
2. Open the control file and parameterize it with your model information. 
   The control file must be named HeatSource_Control.csv 
    
3. Use hs9_setup_model_inputs.py to build template input files or by using commend line. The input files will
   be saved to the input file directory that is specified in the control file).
   ```shell
   cd path\to\model_directory
   hs setup -mi
   ```   
5. Edit the template csv files with your input data. You can use excel although 
   make sure the datetimes are formatted correctly. Save the files as a csv. 

6. Run the model by executing one of the following model python scripts/executables:
   * hs9_run_hydraulics
   * hs9_run_solar 
   * hs9_run_temperature
   
   Or use command line to run the model:
   ```shell
   cd path\to\model_directory
   hs run -t
   ```
   A console should open and you should see the model running.
   
7. Outputs are saved in the output directory (specified in the control file).

## INPUT FILES - GENERAL INFORMATION

1. Control and input files are ASCII comma delimited files.
2. The heat source control file must be named `HeatSource_Control.csv`.
3. The other input files can be named whatever you want 
 (file names are specified in the control file).
4. Do not change the key names in the control file. Only change 
   the VALUE column (column 4).
5. The column header names can be changed but the data needs to be in 
 the correct column number (see below).
6. Use the specified unit and data format identified in the control file 
 and input files. Example mm/dd/yyyy hh:mm is 07/01/2001 16:00
7. An input parameter value that is not applicable may be left blank although all values with float
	data type will be assigned as zero.

Key to model input information:<br>
Required:  Input value required<br>
Required (Not Used):  Input value required but not used other than for model setup (may be changed in future versions)<br>
Optional 1:  Input value optional (can be left blank).<br>
Optional 2:  Input value may be optional based on control file parameterization.<br>

To write blank input files from a python script:
```python
# requires a parameterized control file 
from heatsource9 import BigRedButton

control_file = 'HeatSource_Control.csv'
model_dir = r'C://path/to/model_directory/'

BigRedButton.setup_mi(model_dir, control_file,
                      use_timestamp=True, overwrite=False)
```
To write a blank inputs files from command line:
```shell
hs setup -mi
```
## Using Command Line 

Heat Source can be setup and run directly from command line.

usage: hs <command> [options]
commands:<br>
<div style="padding-left: 10px;">
run Command to run a model with arguments -t | -s | -hy <br>
setup Command to setup a model with arguments -cf | -mi <br>
</div>

run options:<br>
&nbsp; -h, --help &nbsp;&nbsp; show this help message<br>
&nbsp; -t, --temperature &nbsp;&nbsp; Runs a temperature model.<br>
&nbsp; -s, --solar &nbsp;&nbsp; Runs solar routines only.<br>
&nbsp; -hy, --hydraulics &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Runs hydraulics only.<br>
  
setup options:<br>
&nbsp; -h, --help           show this help message<br>
&nbsp; -cf, --control-file  Writes a blank control file.<br>
&nbsp; -mi, --model-inputs  Write blank input files. Control file must already be parameterized.<br>
&nbsp; -t, --timestamp      Use -t to add a timestamp to the file name.<br>
&nbsp; -o, --overwrite      Use -o to overwrite any existing file.
  
other options<br>
  -h, --help            show this help message<br>
  -v                    The heat source version and install directory.<br>
  -md [MODEL_DIR], --model-dir [MODEL_DIR]<br>
                        Path to the model directory. If not used the default is current<br>
                        working directory.

## CONTROL FILE  
HeatSource_Control.csv

The control file is where most of the model operation and initial parameterization is set.

To write a blank template control file from a python script:
```python
from heatsource9 import BigRedButton

control_file = 'HeatSource_Control.csv'
model_dir = r'C://path/to/model_directory/'

BigRedButton.setup_cf(model_dir, control_file,
                      use_timestamp=False, overwrite=False)
```
To write a blank template control file from command line:
```Batchfile
cd path\to\model_directory
hs setup -cf
```

You can also parameterize the control file in python directly using `**kwargs`. Any control file key arguments passed will be parameterized into the control file.
```python
from heatsource9 import BigRedButton
from os.path import join

control_file = 'HeatSource_Control.csv'
model_dir = r'C://path/to/model_directory/'

# Parameterize the control file and write to csv
BigRedButton.setup_cf(model_dir, control_file, use_timestamp=True, overwrite=False,
                      usertxt="This model is an example model",
                      name="example model",
                      inputdir=join(model_dir, "inputs", ""),
                      outputdir=join(model_dir, "outputs", ""),
                      length=1.8,
                      outputkm="all",
                      datastart="05/06/2003",
                      modelstart="07/01/2003",
                      modelend="07/14/2003",
                      dataend="09/21/2003",
                      flushdays=1,
                      offset=-7,
                      dt=1,
                      dx=30,
                      longsample=50,
                      bcfile="bc.csv",
                      inflowsites=4,
                      inflowinfiles="inflow1.csv, inflow2.csv, inflow3.csv, inflow4.csv",
                      inflowkm="1.65, 1.5, 1.3, 0.85",
                      accretionfile="accretion.csv",
                      metsites=4,
                      metfiles="met1.csv, met2.csv, met3.csv, met4.csv",
                      metkm="1.75, 1.45, 1.10, 0.9",
                      calcevap="False",
                      evapmethod="Mass Transfer",
                      wind_a=1.51E-09,
                      wind_b=1.6E-09,
                      calcalluvium="True",
                      alluviumtemp=12.0,
                      morphfile="morphology.csv",
                      lcdatafile="lcdata.csv",
                      lccodefile="lccodes.csv",
                      trans_count=8,
                      transsample_count=4,
                      transsample_distance=8,
                      emergent="True",
                      lcdatainput="Codes",
                      canopy_data="CanopyCover",
                      lcsampmethod="point",
                      heatsource8="False")
```

Below are all the input parameters that must be included in the control file.

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


## ACCRETION INPUT FILE  
UserDefinedFileName.csv

The temperature and flow rates of accretion are defined in this file. 

Accretion flows are inflows that enter the stream over more than one 
stream data node, and typically are subsurface seeps that occur over 
longer distances than discrete subsurface inflows (i.e. a spring).  

When accretion flows are close enough so that more than one occurs 
in a model distance step, the accretion flow rates will be summed and a 
flow based average accretion temperature will be derived and used 
in the mixing calculations.

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`STREAM_ID`|Stream ID|N/A|string|Optional 1|Optional 1|Optional 1|
|2|`NODE_ID`|Node ID|N/A|integer|Required|Required|Required|
|3|`STREAM_KM`|Stream km|kilometers|float|Required|Required|Required|
|4|`INFLOW`|Accretion Inflow|cubic meters/second|float|Optional 1|Required|Required|
|5|`TEMPERATURE`|Accretion Temperature|degrees Celsius|float|Optional 1|Required|Required|
|6|`OUTFLOW`|Withdrawal flow|cubic meters/second|float|Optional 1|Required|Required|


## BOUNDARY CONDITION FILE  
UserDefinedFileName.csv

The stream flow and temperature conditions at the upstream model boundary 
are defined in this file. The boundary conditions are defined at an 
hourly timestep.

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`DATETIME`|The date/time|mm/dd/yyyy hh:mm|string|Optional 1|Required|Required|
|2|`FLOW`|Boundary condition flow|cubic meters/second|float|Optional 1|Required|Required|
|3|`TEMPERATURE`|Boundary condition temperature|degrees Celsius|float|Optional 1|Required|Required|


## METEOROLOGICAL INPUT FILE/S
(formally called Continuous data in heat source 8)
UserDefinedFileName.csv

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`DATETIME`|The date/time|mm/dd/yyyy hh:mm|string|Optional 1|Required|Required|
|2|`CLOUDINESS1`|Cloudiness|decimal fraction (0-1)|float|Optional 1|Optional 2|Optional 2|
|3|`WIND_SPEED1`|Wind Speed|meters/second|float|Optional 1|Optional 2|Optional 2|
|4|`RELATIVE_HUMIDITY1`|Relative Humidity|decimal fraction (0-1)|float|Optional 1|Optional 2|Optional 2|
|5|`AIR_TEMPERATURE1`|Air Temperature|degrees Celsius|float|Optional 1|Optional 2|Optional 2|

Note - multiple csv files may be used for each set of meteorological inputs with the 
format above or all data can be saved in the same file like the example below. 
This is controlled in the control file by designating the number of inputs and 
the input stream km.

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`DATETIME`|The date/time|mm/dd/yyyy hh:mm|string|Optional 1|Required|Required|
|2|`CLOUDINESS1`|Cloudiness at site 1|decimal fraction (0-1)|float|Optional 1|Optional 2|Optional 2|
|3|`WIND_SPEED1`|Wind Speed at site 1|meters/second|float|Optional 1|Optional 2|Optional 2|
|4|`RELATIVE_HUMIDITY1`|Relative Humidity at site 1|decimal fraction (0-1)|float|Optional 1|Optional 2|Optional 2|
|5|`AIR_TEMPERATURE1`|Air Temperature at site 1|degrees Celsius|float|Optional 1|Optional 2|Optional 2|
|6|`CLOUDINESS2`|Cloudiness at site 2|decimal fraction (0-1)|float|Optional 1|Optional 2|Optional 2|
|7|`WIND_SPEED2`|Wind Speed at site 2|meters/second|float|Optional 1|Optional 2|Optional 2|
|8|`RELATIVE_HUMIDITY2`|Relative Humidity at site 2|decimal fraction (0-1)|float|Optional 1|Optional 2|Optional 2|
|9|`AIR_TEMPERATURE2`|Air Temperature at site 2|degrees Celsius|float|Optional 1|Optional 2|Optional 2|

## TRIBUTARY INPUT FILE/S  
Can also be outflows. Use negative flows.  
UserDefinedFileName.csv  

The tributary input files define the inflow/outflow rates and temperatures
at different points along the model stream. Inflows refers to localized 
(non-accretion) type flows such as tributaries, springs, returns, point 
sources, etc. Outflows can be various types of water withdrawals. Outflows 
are input with a negative flow rate. Temperatures for outflows are not 
used by the model.

The number and stream km of the inflow/outflows is defined in the control file.
The flow and temperature are defined at an hourly timestep.  

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`DATETIME`|The date/time|mm/dd/yyyy hh:mm|string|Optional 1|Required|Required|
|2|`FLOW1`|Tributary flow|cubic meters/second|float|Optional 1|Required|Required|
|3|`TEMPERATURE1`|Tributary Temperature|degrees Celsius|float|Optional 1|Required|Required|

Note - multiple csv files may be created for each tributary input with the
format above or all data can be saved in the same file like in the example below.
This is controlled in the control file by designating the number of inputs and
the input stream km.

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`DATETIME`|The date/time|mm/dd/yyyy hh:mm|string|Optional 1|Required|Required|
|2|`FLOW1`|Tributary 1 flow|cubic meters/second|float|Optional 1|Required|Required|
|3|`TEMPERATURE1`|Tributary 1 Temperature|degrees Celsius|float|Optional 1|Required|Required|
|4|`FLOW2`|Tributary 2 flow|cubic meters/second|float|Optional 1|Required|Required|
|5|`TEMPERATURE2`|Tributary 2 Temperature|degrees Celsius|float|Optional 1|Required|Required|

## LANDCOVER CODES FILE  
UserDefinedFileName.csv  

The landcover codes file contains the physical attribute information 
associated with each land cover code. Land cover codes can be alphanumeric 
values. Zero should be avoided as a land cover code. The physical 
attribute such as height, canopy cover, LAI, or overhang must be numeric.
There cannot be skipped rows 
(i.e. rows without information in between rows with information) 
because the model routines see a blank row as the end of the data sequence.

### Canopy Type

land cover canopy information can be input as either canopy cover or 
effective leaf area index. This option is specified in the control file using the key ```canopy_data```.

#### Canopy Cover
Input file formatting when ```canopy_data = "CanopyCover"``` in the control file.

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`NAME`|Land cover Name|N/A|string|Optional 1|Optional 1|Optional 1|
|2|`CODE`|Land cover code|N/A|string|Required|Required (Not Used)|Required|
|3|`HEIGHT`|Land cover height|meters|float|Required|Required (Not Used)|Required|
|4|`CANOPY`|Canopy cover|decimal fraction (0-1)|float|Optional 2|Optional 2|Optional 2|
|5|`OVERHANG`|Overhang|meters|float|Required|Required (Not Used)|Required|


#### LAI
Input file formatting when ```canopy_data = "LAI"``` in the control file.

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`NAME`|Land cover Name|N/A|string|Optional 1|Optional 1|Optional 1|
|2|`CODE`|Land cover code|N/A|string|Required|Required (Not Used)|Required|
|3|`HEIGHT`|Land cover height|meters|float|Required|Required (Not Used)|Required|
|4|`LAI`|Effective Leaf Area Index|dimensionless|float|Optional 2|Optional 2|Optional 2|
|5|`k`|k extinction coefficient|dimensionless|float|Optional 2|Optional 2|Optional 2|

The landcover codes file can be parameterized from script.
```python
from heatsource9.ModelSetup.Inputs import Inputs
from heatsource9.Dieties.IniParamsDiety import IniParams
from heatsource9.Dieties.IniParamsDiety import dtype

control_file = 'HeatSource_Control.csv'
model_dir = r'C://path/to/model_directory/'

# create an input object
inputs = Inputs(model_dir, control_file)

# imports the control file into input object
inputs.import_control_file()

# Parameterize the lccodes input. Uses canopy closure data.
lccodes = [('Active River Channel',100,0,0,0), 
           ('Barren - Clearcut',127,0,0,0), 
           ('Brush',128,1,0.4,0), 
           ('Dominate Coniferous',133,32,0.7,1.5), 
           ('Dominate Broadleaf (Riparian)',149,32,0.5,2), 
           ('Dominate Broadleaf (Upland)',150,32,0.5,2), 
           ('Road Unpaved',255,0,0,0)]

inputs.parameterize_lccodes(lccodes, overwrite=True)
```

## LANDCOVER DATA  
(formally called TTools in heatsource 8)  
UserDefinedFileName.csv  

This file defines land cover information. This data can be derived 
from geospatial data using TTools.

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`STREAM_ID`|Stream ID|N/A|string|Optional 1|Optional 1|Optional 1|
|2|`NODE_ID`|Node ID|N/A|integer|Required|Required|Required|
|3|`STREAM_KM`|Stream km|kilometer|float|Required|Required (Not Used)|Required|
|4|`LONGITUDE`|Node Longitude|decimal degrees|float|Required|Required (Not Used)|Required|
|5|`LATITUDE`|Node Latitude|decimal degrees|float|Required|Required (Not Used)|Required|
|6|`TOPO_W`|Topographic shade angle to the west|degrees|float|Required|Required (Not Used)|Required|
|7|`TOPO_S`|Topographic shade angle to the south|degrees|float|Required|Required (Not Used)|Required|
|8|`TOPO_E`|Topographic shade angle to the east|degrees|float|Required|Required (Not Used)|Required|

### Land Cover Data Input Type
Land cover information can be input into the model in two different ways:
Using codes or values. If using codes unique land cover attribute information is represented as a unique code.
The land cover attribute  associated to each code is parameterized in the Land Cover Codes file in long format.
If using values, the land cover attribute information for each transect sample is explicitly in the land cover data file.

The land cover data input type is identified in the control file.

#### Codes
When ```lcdatainput = "Codes"```, the following columns will be used after column 8:

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|multiple|`LC_T_S`|Landcover code on transect T for sample S|N/A|string|Required|Required (Not Used)|Required|
|multiple|`ELE_T_S`|Elevation on transect T for sample S|meters|float|Required|Required (Not Used)|Required|

#### Values
When ```lcdatainput = "Values"```, and ```canopy_data = "CanopyClosure"``` the following columns will be used after column 8:

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|multiple|`HT_T_S`|Land cover height on transect T for sample S|N/A|string|Required|Required (Not Used)|Required|
|multiple|`ELE_T_S`|Elevation on transect T for sample S|meters|float|Required|Required (Not Used)|Required|
|multiple|`CAN_T_S`|Canopy cover on transect T for sample S|decimal fraction (0-1)|float|Required|Required (Not Used)|Required|
|multiple|`OH_T_S`|Overhang on transect T for sample S|meters|float|Required|Required (Not Used)|Required|

When ```lcdatainput = "Values"```, and ```canopy_data = "LAI"``` the following columns will be used after column 8:

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|multiple|`HT_T_S`|Land cover height on transect T for sample S|N/A|string|Required|Required (Not Used)|Required|
|multiple|`ELE_T_S`|Elevation on transect T for sample S|meters|float|Required|Required (Not Used)|Required|
|multiple|`LAI_T_S`|Effective Leaf Area Index on transect T for sample S|dimensionless|float|Required|Required (Not Used)|Required|
|multiple|`k_T_S`|k extinction coefficient on transect T for sample S|dimensionless|float|Required|Required (Not Used)|Required|
|multiple|`OH_T_S`|Overhang on transect T for sample S|meters|float|Required|Required (Not Used)|Required|

Note - the number of columns are dependent on the number of transects and samples specified in the control file.

## MORPHOLOGY DATA FILE  
UserDefinedFileName.csv  

This file defines channel morphology and substrate information.
Refer to the user manual for more information about each parameter.

|COLUMN NUMBER |COLUMN NAME |DESCRIPTION |UNITS |DATA TYPE |SOLAR RUNS |HYDRAULIC RUNS|TEMPERATURE RUNS |
|-------------:|:-----------|:-----------|:-----|:---------|:---------:|:------------:|:---------------:|
|1|`STREAM_ID`|Stream ID|N/A|string|Optional 1|Optional 1|Optional 1|
|2|`NODE_ID`|Node ID|N/A|integer|Required|Required|Required|
|3|`STREAM_KM`|Stream km|kilometers|float|Required|Required|Required|
|4|`ELEVATION`|Stream Elevation|meters|float|Required|Required (Not Used)|Required|
|5|`GRADIENT`|Stream Gradient|meters/meters|float|Optional 1|Required|Required|
|6|`BOTTOM_WIDTH`|Bottom Width|meters|float|Optional 1|Required|Required|
|7|`CHANNEL_ANGLE_Z`|Channel Angle z |meters/meters|float|Optional 1|Required|Required|
|8|`MANNINGS_n`|Manning's n|seconds/meter|float|Optional 1|Required|Required|
|9|`SED_THERMAL_CONDUCTIVITY`|Sediment Thermal Conductivity|watts/meters/degrees Celsius|float|Optional 1|Required|Required|
|10|`SED_THERMAL_DIFFUSIVITY`|Sediment Thermal Diffusivity|square centimeters/second|float|Optional 1|Required|Required|
|11|`SED_HYPORHEIC_THICKNESSS`|Hyporheic Zone Thickness|meters|float|Optional 1|Required|Required|
|12|`HYPORHEIC_PERCENT`|Percent Hyporheic Exchange|decimal fraction (0-1)|float|Optional 1|Required|Required|
|13|`POROSITY`|Porosity|decimal fraction (0-1)|float|Optional 1|Required|Required|

## LICENSE
GNU General Public License v3 (GPLv3)

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

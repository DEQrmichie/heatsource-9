![Wheel Build](https://github.com/DEQrmichie/heatsource-9/actions/workflows/build_wheels.yml/badge.svg)

Heat Source 9
-------------
Current Version: heatsource 9.0.0

## 1.0 ABOUT 
Heat Source is a computer model used by the Oregon Department of  Environmental Quality to simulate stream 
thermodynamics and hydraulic routing. It was originally developed by Matt Boyd in 1996 as a [Master's Thesis][1] 
at Oregon State University in the Departments of Bioresource Engineering and Civil Engineering. Since then, 
it has grown and changed significantly. Oregon DEQ currently maintains the Heat Source methodology and computer 
programming. Appropriate model use and application are the sole responsibility of the user. 

Heat Source 7-8 and user manual: 
http://www.oregon.gov/deq/wq/tmdls/Pages/TMDLs-Tools.aspx

Authors: Matt Boyd, Brian Kasper, Terra Metta, Ryan Michie, Dan Turner

Contact: Ryan Michie, ryan.michie @ deq.oregon.gov

[1]: http://ir.library.oregonstate.edu/xmlui/handle/1957/27036

## 2.0 INSTALL

There are two options for installing and running Heat Source 9:

1. Download the [windows executables][2]. Place the executables in the directory with your model files. Double-click 
the executable to run the model. That's it. Python installation is not required. Optionally add hs.exe to the
user PATH and run the model from the command line (see section 5.1). These executables were developed on Windows 10. 
They have not been tested on other versions of Windows. 

2. Install the model as a Python package. Requires install of Python 3.10 or higher.
https://www.python.org/downloads/. This may be better option for Mac or Linux users.

After Python has been installed, install the Heat Source package from the command line using pip.
```shell
# This command installs heat source version 9.0.0 directly from the GitHub repository.
pip install "git+https://github.com/DEQrmichie/heatsource-9@v9.0.0"
```
Alternatively, the package can be installed by downloading the [heat source python wheel][3] appropriate to 
your OS platform and python version. Python wheels have been built to support Windows, Mac, and Linux.
DEQ uses Windows so other platforms have limited testing. We've heard folks having success on both Mac and Linux.
After downloading the wheel, install from the command line using pip.
```shell
# These commands are for windows
cd path\to\directory_where_the_heatsource9_wheel_was_saved\
py -m pip install <name of wheel file>

# Installs the Python 3.12 heat source wheel for windows in the local directory
py -m pip3 install heatsource9-9.0.0-cp312-cp312-win_amd64.whl --user

# Installs the Python 3.12 heat source wheel for windows in the global directory
py -m pip3 install heatsource9-9.0.0-cp312-cp312-win_amd64.whl
 ```
[2]: https://github.com/DEQrmichie/heatsource-9/releases/download/v9.0.0/HS9_Windows_Executables_v9.0.0.zip
[3]: https://github.com/DEQrmichie/heatsource-9/releases/tag/v9.0.0

## 3.0 QUICK STEPS TO GET GOING

1. Generate a template control file by executing *hs_setup_control_file* or by using command line.
   ```shell
   cd path\to\model_directory
   hs setup -cf
   ```
   Make sure the control file (HeatSource_Control.[xlsx|csv]) and the model run
   executables are in the same directory.

2. Open the control file and parameterize it with your model information. 
   The control file must be named HeatSource_Control.[xlsx|csv]. Make sure to set the number of meteorological sites with `metsites`
   and the number of tributary inflow sites with `tribsites`.
    
3. Use *hs9_setup_model_inputs* to build template input files or by using command line. The input files will
   be saved to the input file directory that is specified in the control file. This step also creates the meteorological sites file
   (HeatSource_Met_Sites) and tributary sites file
   (HeatSource_Tributary_Sites) in the model directory when they do not already exist.
   ```shell
   cd path\to\model_directory
   hs setup -mi
   ```   
4. Open HeatSource_Met_Sites.[xlsx|csv] and HeatSource_Tributary_Sites.[xlsx|csv] and enter the site settings
   for each meteorological site and tributary site. 

5. Add model input data to the required template input files.

6. Run the model by executing one of the following model executables or
   python scripts:
   * `hs9_run_hydraulics`
   * `hs9_run_solar`
   * `hs9_run_temperature`
   
   Or use the command line:
   ```shell
   cd path\to\model_directory
   hs run -hy
   hs run -s
   hs run -t
   ```
   
7. Outputs are saved in the output directory (specified in the control file).

## 4.0 CSV MODE
By default, the template control file and input files are written as Excel files (.xlsx). Some users may want to
use CSV formatted files instead. To switch from Excel to CSV, resave the control file as a CSV (UTF-8 Unicode).
A template CSV control file can also be generated from command line using the following:
```shell
hs setup -cf -csv
```
The model will read and write input files using the same format as the control file. If the control file is formatted
as a CSV, the input files must also be CSV. If the control file is an Excel file (.xlsx), the input files must also be
saved as Excel files. Model output files are always written as CSV (UTF-8 Unicode) files.

## 5.0 MODEL RUN OPTIONS

There are multiple workflows available for setting up and running a model.

### 5.1 Windows Executables

[Windows executables][2] have been compiled from the source code and can be used to run the model on Windows machines. 
Place the executables in the model directory with your model files. Double-click the executable to setup or run 
the model.

Setup the model executing using one of the following:
   * `hs9_setup_control_file.exe` writes a blank template control file.
   * `hs9_setup_model_inputs.exe` writes blank template input files from a parameterized control file.

Run the model by executing one of the following:
   * `hs9_run_hydraulics.exe` runs a hydraulics only model.
   * `hs9_run_solar.exe` runs a solar and shade only model.
   * `hs9_run_temperature.exe` runs a temperature model, which include both hydraulics and solar.

### 5.2 Command Line

Heat Source can be set up and run directly from the command line. Using command line requires the Python package 
be installed or having hs.exe on PATH.

Usage: 
``` shell
hs <command> [options]
```
Command forms:
``` shell
hs [-md MODEL_DIR] run -t | -s | -hy
hs [-md MODEL_DIR] setup (-cf | -mi) [setup options]
```

Global options:
``` shell
-h / --help              : show help
-md / --model-dir        : path to the model directory
-v / --version           : print Heat Source version
```

Run types:
``` shell
-t / --temperature       : run the full temperature model
-s / --solar             : run solar routines only
-hy / --hydraulics       : run hydraulics only
```

Setup types:
``` shell
-cf / --control-file     : write a blank control file template
-mi / --model-inputs     : write blank input files from a parameterized control file
```

Setup options:
``` shell
-csv / --csv-mode        : with -cf, write a CSV control file instead of XLSX
-set / --set KEY=VALUE   : with -cf, set one control file key, repeat as needed
-t / --timestamp         : add a timestamp to the file name
-o / --overwrite         : overwrite existing files (default is to keep existing files)
```

Usage examples:
``` shell
hs setup -cf -md /path/to/model_dir
hs setup -cf -set name="Example Model" -set length=7.5
hs setup -mi -o
hs run -t
```

### 5.3 Using Python

Heat Source can be set up and run directly using python scripts. The heatsource9 Python package must
be installed.

To write a blank XLSX control file from a python script:
```python
from heatsource9 import setup

model_dir = r"C://path/to/model_directory/"

setup.setup_cf(model_dir)
```

To parameterize a control file from a python script:
```python
from heatsource9 import setup

model_dir = r"C://path/to/model_directory/"

setup.setup_cf(
    model_dir,
    name="Example Model",
    length=7.5,
    dt=1,
)
```

If the control file already exists and `overwrite=False`, `setup_cf(...)` leaves it unchanged. Use `overwrite=True` to rewrite the control file from a blank template and then apply keyword values.

To write blank XLSX input files from a python script:
```python
# requires a parameterized control file 
from heatsource9 import setup

control_file = 'HeatSource_Control.xlsx'
model_dir = r'C://path/to/model_directory/'

setup.setup_mi(model_dir=model_dir, control_file=control_file,
               use_timestamp=True, overwrite=False)
```
The file format for the input templates are CSV or XLSX from the control file extension.


To run a temperature model from a python script:
```python
from heatsource9 import run

control_file = "HeatSource_Control.xlsx"
model_dir = r"C://path/to/model_directory/"

run.temperature(model_dir, control_file)
```

## 6.0 MODEL INPUT FILES
The following table summarizes what input files are needed to run each type of model. More specific details about the 
format and content of each input file is included in the sections below.

Field details:

| INPUT FILE          | FILE NAME                  | xlsx SHEET NAME      |
|:--------------------|:---------------------------|:---------------------|
| CONTROL FILE        | HeatSource_Control         | Control Settings     |
| MET SITE FILE       | HeatSource_Met_Sites       | Meteorological Sites |
| TRIBUTARY SITE FILE | HeatSource_Tributary_Sites | Tributary Sites      |
| ACCRETION           | User Defined               | Accretion Flow       |
| BOUNDARY CONDITION  | User Defined               | Boundary Conditions  |
| METEOROLOGICAL DATA | User Defined               | Meteorological Data  |
| TRIBUTARY DATA      | User Defined               | Tributary Data       |
| LAND COVER CODES    | User Defined               | Land Cover Codes     |
| LAND COVER DATA     | User Defined               | Land Cover Data      |
| MORPHOLOGY DATA     | User Defined               | Morphology Data      |

Model run requirements:

| INPUT FILE          | SOLAR RUNS           | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:--------------------|:---------------------|:---------------|:-----------------|
| CONTROL FILE        | Required             | Required       | Required         |
| MET SITE FILE       | Required             | Optional       | Required         |
| TRIBUTARY SITE FILE | Optional             | Required       | Required         |
| ACCRETION           | Optional             | Optional       | Optional         |
| BOUNDARY CONDITION  | Optional             | Required       | Required         |
| METEOROLOGICAL DATA | Required<sup>1</sup> | Optional       | Required         |
| TRIBUTARY DATA      | Optional             | Required       | Required         |
| LAND COVER CODES    | Required             | Optional       | Required         |
| LAND COVER DATA     | Required             | Optional       | Required         |
| MORPHOLOGY DATA     | Required             | Required       | Required         |

Key to model input information:

 -   Required: Input value required.
 -   Optional: File or input value optional. Note control file value for `tribsites` must be 0 if there are no files.
 
 <sup>1 Cloudiness required, other met fields can be blank.</sup>

General Information
1. The control file and input files are set up as Excel (`.xlsx`) files by default. They can also be CSV (UTF-8) comma delimited files.
2. The Heat Source control file must be named HeatSource_Control.[xlsx|csv]. The meteorological sites file must be 
   named HeatSource_Met_Sites.[xlsx|csv]. The tributary sites file must be named HeatSource_Tributary_Sites.[xlsx|csv]. 
   The other input files can use any name.
3. The column header names can be changed but the data needs to be in the correct column number.
4. Use the specified units and data formats identified in the control file and input files.
   Example `yyyy-mm-dd hh:mm` is `2001-07-01 16:00`.
5. An input parameter value that is optional may be left blank. Optional parameters with float data type are assigned as zero and do not impact model results.
6. Control key `accretionfile` is optional for hydraulics and temperature runs, when the key is set, section 6.6 field requirements still apply and `NODE_ID` and `STREAM_KM` remain required.

### 6.1 CONTROL FILE  
File name: `HeatSource_Control.[xlsx|csv]`

xlsx sheet name: `Control Settings`

The control file is where most of the model operation and initial parameterization is set. 
Do not change the key names in the control file. Only change the VALUE column (column 4). The heat source control file 
must be named HeatSource_Control.[xlsx|csv].

The control file also defines the number of meteorological and tributary sites used in HeatSource_Met_Sites.[xlsx|csv] 
and HeatSource_Tributary_Sites.[xlsx|csv]. The `metsites` and `tribsites` values must be set before creating these 
files with the setup executable or by using the command `hs setup -mi`.

To write a blank template control file from a python script:
```python
from heatsource9 import setup

control_file = 'HeatSource_Control.xlsx'
model_dir = r'C://path/to/model_directory/'

setup.setup_cf(model_dir=model_dir, control_file=control_file,
               use_timestamp=False, overwrite=False)
```
To write a blank template control file from command line:
```Batchfile
cd path\to\model_directory
hs setup -cf
```

The control file can also be parameterized in python directly using `**kwargs`.
Any control file key arguments passed will be written into the output control file when the control file is first created, or when `overwrite=True`.
Unknown keys raise an error.
```python
from heatsource9 import setup
from os.path import join

control_file = 'HeatSource_Control.xlsx'
model_dir = r'C://path/to/model_directory/'

# Parameterize the control file and write to XLSX
setup.setup_cf(model_dir=model_dir, control_file=control_file,
               use_timestamp=True, overwrite=False,
               usertxt="This model is an example model",
               name="example model",
               inputdir=join(model_dir, "inputs", ""),
               outputdir=join(model_dir, "outputs", ""),
               length=1.8,
               outputkm="all",
               datastart="2003-05-06",
               modelstart="2003-07-01",
               modelend="2003-07-14",
               dataend="2003-09-21",
               flushdays=1,
               offset=-7,
               dt=1,
               dx=30,
               outputdt=60,
               longsample=50,
               bcfile="bc.csv",
               tribsites=4,
               accretionfile="accretion.csv",
               metsites=4,
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
               canopy_data="CanopyCover",
               lcsampmethod="point",
               heatsource8="False")
```

To parameterize a control file from command line:
```shell
hs setup -cf -set name="My Model" -set length=10.5 -set dt=15
```

Below are all the input parameters that must be included in the control file.

| LINE | PARAMETER                                         | KEY                  | VALUE |
|-----:|:--------------------------------------------------|:---------------------|-------|
|    2 | Model Description/User Notes                      | usertxt              |       |
|    3 | Simulation Name                                   | name                 |       |
|    4 | Input Directory Path                              | inputdir             |       |
|    5 | Output Directory Path                             | outputdir            |       |
|    6 | Stream Length (kilometers)                        | length               |       |
|    7 | Output Stream Kilometers                          | outputkm             |       |
|    8 | Data Start Date (yyyy-mm-dd)                      | datastart            |       |
|    9 | Modeling Start Date (yyyy-mm-dd)                  | modelstart           |       |
|   10 | Modeling End Date (yyyy-mm-dd)                    | modelend             |       |
|   11 | Data End Date (yyyy-mm-dd)                        | dataend              |       |
|   12 | Flush Initial Condition (days)                    | flushdays            |       |
|   13 | Time Offset From UTC (hours)                      | offset               |       |
|   14 | Model Time Step (minutes)                         | dt                   |       |
|   15 | Model Distance Step (meters)                      | dx                   |       |
|   16 | Output Time Step (minutes)                        | outputdt             |       |
|   17 | Longitudinal Stream Sample Distance (meters)      | longsample           |       |
|   18 | Boundary Condition Input File Name                | bcfile               |       |
|   19 | Tributary Inflow Sites                            | tribsites            |       |
|   20 | Accretion Input File Name                         | accretionfile        |       |
|   21 | Meteorological Data Sites                         | metsites             |       |
|   22 | Include Evaporation Losses From Flow (True/False) | calcevap             |       |
|   23 | Evaporation Method (Mass Transfer/Penman)         | evapmethod           |       |
|   24 | Wind Function Coefficient a                       | wind_a               |       |
|   25 | Wind Function Coefficient b                       | wind_b               |       |
|   26 | Include Deep Alluvium Temperature (True/False)    | calcalluvium         |       |
|   27 | Deep Alluvium Temperature (Celsius)               | alluviumtemp         |       |
|   28 | Morphology Input Data File Name                   | morphfile            |       |
|   29 | Land Cover Input Data File Name                   | lcdatafile           |       |
|   30 | Land Cover Codes Input File Name                  | lccodefile           |       |
|   31 | Number Of Transects Per Node                      | trans_count          |       |
|   32 | Number Of Samples Per Transect                    | transsample_count    |       |
|   33 | Distance Between Transect Samples (meters)        | transsample_distance |       |
|   34 | Account For Emergent Veg Shading (True/False)     | emergent             |       |
|   35 | Canopy Data Type (LAI/CanopyCover)                | canopy_data          |       |
|   36 | Land Cover Sample Method (point/zone)             | lcsampmethod         |       |
|   37 | Use Heat Source 8 Land Cover Methods (True/False) | heatsource8          |       |

Some control file keys are optional and may be left blank depending on the model run type.
Model run requirements:

| KEY                    | INPUT FILE SETUP | SOLAR RUNS           | HYDRAULIC RUNS       | TEMPERATURE RUNS        |
|:-----------------------|:-----------------|:---------------------|:---------------------|:------------------------|
| `usertxt`              | Optional         | Optional             | Optional             | Optional                |
| `name`                 | Optional         | Optional             | Optional             | Optional                |
| `inputdir`             | Required         | Required             | Required             | Required                |
| `outputdir`            | Optional         | Required             | Required             | Required                |
| `length`               | Required         | Required             | Required             | Required                |
| `outputkm`             | Optional         | Required             | Required             | Required                |
| `datastart`            | Required         | Required             | Required             | Required                |
| `modelstart`           | Optional         | Required             | Required             | Required                |
| `modelend`             | Optional         | Required             | Required             | Required                |
| `dataend`              | Required         | Required             | Required             | Required                |
| `flushdays`            | Optional         | Optional             | Required             | Required                |
| `offset`               | Optional         | Required             | Optional             | Required                |
| `dt`                   | Optional         | Required             | Required             | Required                |
| `dx`                   | Optional         | Required             | Required             | Required                |
| `outputdt`             | Optional         | Optional<sup>1</sup> | Optional<sup>1</sup> | Optional<sup>1</sup>    |
| `longsample`           | Required         | Required             | Required             | Required                |
| `bcfile`               | Optional         | Optional             | Required             | Required                |
| `tribsites`            | Optional         | Optional             | Required             | Required                |
| `accretionfile`        | Optional         | Optional             | Optional             | Optional                |
| `metsites`             | Optional         | Required             | Optional             | Required                |
| `calcevap`             | Optional         | Optional             | Optional             | Required                |
| `evapmethod`           | Optional         | Optional             | Optional             | Required                |
| `wind_a`               | Optional         | Optional             | Optional             | Required                |
| `wind_b`               | Optional         | Optional             | Optional             | Required                |
| `calcalluvium`         | Optional         | Optional             | Optional             | Required                |
| `alluviumtemp`         | Optional         | Optional             | Optional             | Conditional<sup>2</sup> |
| `morphfile`            | Optional         | Required             | Required             | Required                |
| `lcdatafile`           | Optional         | Required             | Optional             | Required                |
| `lccodefile`           | Optional         | Required             | Optional             | Required                |
| `trans_count`          | Optional         | Required             | Optional             | Required                |
| `transsample_count`    | Optional         | Required             | Optional             | Required                |
| `transsample_distance` | Optional         | Required             | Optional             | Required                |
| `emergent`             | Optional         | Required             | Optional             | Required                |
| `canopy_data`          | Optional         | Required             | Optional             | Required                |
| `lcsampmethod`         | Optional         | Required             | Optional             | Required                |
| `heatsource8`          | Optional         | Required             | Optional             | Required                |

<sup>1 `outputdt` defaults to 60 if left blank.</sup>
<sup>2 `alluviumtemp` needs to be set if ```calcalluvium = "True"```, otherwise if can be left blank.</sup>

### 6.2 METEOROLOGICAL SITE FILE
File name: `HeatSource_Met_Sites.[xlsx|csv]`

xlsx sheet name: `Meteorological Sites`

The met sites file defines the meteorological input locations used by the model. Each row represents one 
meteorological site and links that site to a meteorological input file, a stream kilometer, and a wind
measurement height.

Field details:

| COLUMN NUMBER | COLUMN NAME  | DESCRIPTION                               | UNITS      | DATA TYPE |
|:-------------:|:-------------|:------------------------------------------|:-----------|:----------|
|       1       | `COLID`      | Internal column ordering ID for each site | N/A        | integer   |
|       2       | `MET_NAME`   | Meteorological site name or details       | N/A        | string    |
|       3       | `STREAM_KM`  | Stream kilometer for the met site         | kilometers | float     |
|       4       | `FILE_NAME`  | Meteorological input data file name       | N/A        | string    |
|       5       | `MET_HEIGHT` | Wind measurement height above ground      | meters     | float     |

Model run requirements:

| COLUMN NUMBER | COLUMN NAME  | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:-------------|:----------:|:--------------:|:----------------:|
|       1       | `COLID`      |  Required  |    Optional    |     Required     |
|       2       | `MET_NAME`   |  Optional  |    Optional    |     Optional     |
|       3       | `STREAM_KM`  |  Required  |    Optional    |     Required     |
|       4       | `FILE_NAME`  |  Required  |    Optional    |     Required     |
|       5       | `MET_HEIGHT` |  Optional  |    Optional    |     Required     |

Example of met sites file setup with multiple excel files for each meteorological input.

| COLID | MET_NAME   | STREAM_KM | FILE_NAME      | MET_HEIGHT |
|:-----:|:-----------|----------:|:---------------|-----------:|
|   1   | Upper Site |      5.45 | met_site1.xlsx |        2.0 |
|   2   | Lower Site |      0.30 | met_site2.xlsx |        2.0 |

Example of met sites file setup with single excel file for all meteorological inputs.

| COLID | MET_NAME   | STREAM_KM | FILE_NAME | MET_HEIGHT |
|:-----:|:-----------|----------:|:----------|-----------:|
|   1   | Upper Site |      5.45 | met.xlsx  |        2.0 |
|   2   | Lower Site |      0.30 | met.xlsx  |        2.0 |

Notes:
- `COLID` is used to define the column ordering for data from each meteorological site when a single input file is used. 
For example, data for `COLID = 1` will be stored in the first group of meteorological data columns in the input file, 
while `COLID = 2` will be stored in the second group.
- Each row defines one meteorological site.
- `FILE_NAME` identifies the meteorological input file used for that site.

### 6.3 METEOROLOGICAL DATA INPUT FILE/S
File name: `UserDefinedFileName.[xlsx|csv]`

xlsx sheet name: `Meteorological Data`
(formally called Continuous data in Heat Source 8)

The meteorological data input file contains the hourly meteorological data used by the model. 
Only solar and temperature model runs require meteorological data. Solar runs only require cloudiness.

Model run requirements:

| COLUMN NUMBER | COLUMN NAME          | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:---------------------|:----------:|:--------------:|:----------------:|
|       1       | `DATETIME`           |  Required  |    Optional    |     Required     |
|       2       | `CLOUDINESS1`        |  Required  |    Optional    |     Required     |
|       3       | `WIND_SPEED1`        |  Optional  |    Optional    |     Required     |
|       4       | `RELATIVE_HUMIDITY1` |  Optional  |    Optional    |     Required     |
|       5       | `AIR_TEMPERATURE1`   |  Optional  |    Optional    |     Required     |

The number of meteorological input sites is defined in the control file using `metsites`. The input meteorological 
file names  and number of files is defined in `HeatSource_Met_Sites` with the `FILE_NAME` column. There can be 
one site per file, or all sites can be included in a single input file. The examples below show the different setup.

The example below shows the columns needed if one file is used for each meteorological input.

Field details:

| COLUMN NUMBER | COLUMN NAME          | DESCRIPTION       | UNITS                  | DATA TYPE |
|:-------------:|:---------------------|:------------------|:-----------------------|:----------|
|       1       | `DATETIME`           | The date/time     | yyyy-mm-dd hh:mm       | string    |
|       2       | `CLOUDINESS1`        | Cloudiness        | decimal fraction (0-1) | float     |
|       3       | `WIND_SPEED1`        | Wind speed        | meters/second          | float     |
|       4       | `RELATIVE_HUMIDITY1` | Relative humidity | decimal fraction (0-1) | float     |
|       5       | `AIR_TEMPERATURE1`   | Air temperature   | degrees Celsius        | float     |

The example below shows the columns needed if one input data file is used for all meteorological 
sites (two sites in this case).

Field details:

| COLUMN NUMBER | COLUMN NAME          | DESCRIPTION                 | UNITS                  | DATA TYPE |
|:-------------:|:---------------------|:----------------------------|:-----------------------|:----------|
|       1       | `DATETIME`           | The date/time               | yyyy-mm-dd hh:mm       |  string   |
|       2       | `CLOUDINESS1`        | Cloudiness at site 1        | decimal fraction (0-1) |   float   |
|       3       | `WIND_SPEED1`        | Wind Speed at site 1        | meters/second          |   float   |
|       4       | `RELATIVE_HUMIDITY1` | Relative Humidity at site 1 | decimal fraction (0-1) |   float   |
|       5       | `AIR_TEMPERATURE1`   | Air Temperature at site 1   | degrees Celsius        |   float   |
|       6       | `CLOUDINESS2`        | Cloudiness at site 2        | decimal fraction (0-1) |   float   |
|       7       | `WIND_SPEED2`        | Wind Speed at site 2        | meters/second          |   float   |
|       8       | `RELATIVE_HUMIDITY2` | Relative Humidity at site 2 | decimal fraction (0-1) |   float   |
|       9       | `AIR_TEMPERATURE2`   | Air Temperature at site 2   | degrees Celsius        |   float   |

### 6.4 TRIBUTARY SITE FILE
File name: `HeatSource_Tributary_Sites.[xlsx|csv]`

xlsx sheet name: `Tributary Sites`

The tributary sites file defines the tributary input locations used by the model. The number of 
tributary input sites is defined in the control file using `tribsites`. A tributary input can be used to
configure a point source discharge, or any other inflow. A tributary site can also have a negative flow and be used 
to represent a timeseries of water withdrawals. Each row in the site file represents one tributary site and links 
that site to a tributary input file and a stream kilometer.

Field details:

| COLUMN NUMBER | COLUMN NAME | DESCRIPTION                               | UNITS      | DATA TYPE |
|:-------------:|:------------|:------------------------------------------|:-----------|:----------|
|       1       | `COLID`     | Internal column ordering ID for each site | N/A        | integer   |
|       2       | `TRIB_NAME` | Tributary site name or details            | N/A        | string    |
|       3       | `STREAM_KM` | Stream kilometer for the tributary site   | kilometers | float     |
|       4       | `FILE_NAME` | Tributary input data file name            | N/A        | string    |

Model run requirements:

| COLUMN NUMBER | COLUMN NAME | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:------------|:----------:|:--------------:|:----------------:|
|       1       | `COLID`     |  Optional  |    Required    |     Required     |
|       2       | `TRIB_NAME` |  Optional  |    Optional    |     Optional     |
|       3       | `STREAM_KM` |  Optional  |    Required    |     Required     |
|       4       | `FILE_NAME` |  Optional  |    Required    |     Required     |

Example of tributary sites file setup with multiple excel files for each tributary input.

| COLID | TRIB_NAME      | STREAM_KM | FILE_NAME  |
|:-----:|:---------------|----------:|:-----------|
|   1   | Quarter Branch |      4.05 | trib1.xlsx |
|   2   | Dry Creek      |      3.35 | trib2.xlsx |

Example of tributary sites file setup with single excel file for all tributary inputs.

| COLID | TRIB_NAME      | STREAM_KM | FILE_NAME  |
|:-----:|:---------------|----------:|:-----------|
|   1   | Quarter Branch |      4.05 | tribs.xlsx |
|   2   | Dry Creek      |      3.35 | tribs.xlsx |

Notes:
- `COLID` is used to define the column ordering for data from each tributary site when a single input file is used. 
For example, data for `COLID = 1` will be stored in the first group of tributary data columns in the input file, 
while `COLID = 2` will be stored in the second group.
- `FILE_NAME` identifies the tributary input file used for that site.

### 6.5 TRIBUTARY DATA INPUT FILE/S
File name: `UserDefinedFileName.[xlsx|csv]`

xlsx sheet name: `Tributary Data`

The tributary input files define the inflow and outflow rates and temperatures at different points along the 
model stream. Inflows refer to localized, non accretion type flows such as tributaries, springs, returns, and 
point sources. Outflows can be various types of water withdrawals. Outflows are input with a negative flow rate. 
Temperatures for outflows are not used by the model. The flow and temperature are defined at an hourly timestep.

If `tribsites` > 0 in the control file, hydraulic and temperature model runs require the tributary data be setup.

Model run requirements if `tribsites` > 0:

| COLUMN NUMBER | COLUMN NAME    | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:---------------|:----------:|:--------------:|:----------------:|
|       1       | `DATETIME`     |  Optional  |    Required    |     Required     |
|       2       | `FLOW1`        |  Optional  |    Required    |     Required     |
|       3       | `TEMPERATURE1` |  Optional  |    Required    |     Required     |

The number of tributary input sites is defined in the control file using `tribsites`. The input tributary file names 
and number of files is defined in `HeatSource_Tributary_Sites` with the `FILE_NAME` column. There can be 
one site per file, or all sites can be included in a single input file. The examples below show the different setup.

The example below shows the columns needed if one file is used for each meteorological input.

Field details:

| COLUMN NUMBER | COLUMN NAME    | DESCRIPTION           | UNITS               | DATA TYPE |
|:-------------:|:---------------|:----------------------|:--------------------|:----------|
|       1       | `DATETIME`     | The date/time         | yyyy-mm-dd hh:mm    | string    |
|       2       | `FLOW1`        | Tributary flow        | cubic meters/second | float     |
|       3       | `TEMPERATURE1` | Tributary temperature | degrees Celsius     | float     |


The example below shows the columns needed if one input data file is used for all tributary 
sites (two sites in this case).

Field details:

| COLUMN NUMBER | COLUMN NAME    | DESCRIPTION             | UNITS               | DATA TYPE |
|:-------------:|:---------------|:------------------------|:--------------------|:----------|
|       1       | `DATETIME`     | The date/time           | yyyy-mm-dd hh:mm    |  string   |
|       2       | `FLOW1`        | Tributary 1 flow        | cubic meters/second |   float   |
|       3       | `TEMPERATURE1` | Tributary 1 temperature | degrees Celsius     |   float   |
|       4       | `FLOW2`        | Tributary 2 flow        | cubic meters/second |   float   |
|       5       | `TEMPERATURE2` | Tributary 2 temperature | degrees Celsius     |   float   |

### 6.6 ACCRETION INPUT FILE  
File name: `UserDefinedFileName.[xlsx|csv]`

xlsx sheet name: `Accretion Flow`

The temperature and flow rates of accretion are defined in this file.

Accretion flows are inflows that enter the stream over more than one 
stream data node, and typically are subsurface seeps that occur over 
longer distances than discrete subsurface inflows (i.e. a spring).  

When accretion flows are close enough so that more than one occurs 
in a model distance step, the accretion flow rates will be summed and a 
flow based average accretion temperature will be derived and used 
in the mixing calculations.

Field details:

| COLUMN NUMBER | COLUMN NAME   | DESCRIPTION           | UNITS               | DATA TYPE |
|:-------------:|:--------------|:----------------------|:--------------------|:----------|
|       1       | `STREAM_ID`   | Stream ID             | N/A                 | string    |
|       2       | `NODE_ID`     | Node ID               | N/A                 | integer   |
|       3       | `STREAM_KM`   | Stream km             | kilometers          | float     |
|       4       | `INFLOW`      | Accretion Inflow      | cubic meters/second | float     |
|       5       | `TEMPERATURE` | Accretion Temperature | degrees Celsius     | float     |
|       6       | `OUTFLOW`     | Withdrawal flow       | cubic meters/second | float     |

Model run requirements:

| COLUMN NUMBER | COLUMN NAME   | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:--------------|:----------:|:--------------:|:----------------:|
|       1       | `STREAM_ID`   |  Optional  |    Optional    |     Optional     |
|       2       | `NODE_ID`     |  Optional  |    Required    |     Required     |
|       3       | `STREAM_KM`   |  Optional  |    Required    |     Required     |
|       4       | `INFLOW`      |  Optional  |    Required    |     Required     |
|       5       | `TEMPERATURE` |  Optional  |    Required    |     Required     |
|       6       | `OUTFLOW`     |  Optional  |    Required    |     Required     |


### 6.7 BOUNDARY CONDITION FILE  
File name: `UserDefinedFileName.[xlsx|csv]`

xlsx sheet name: `Boundary Conditions`

The stream flow and temperature conditions at the upstream model boundary
are defined in this file. The boundary conditions are defined at an 
hourly timestep.

Field details:

| COLUMN NUMBER | COLUMN NAME   | DESCRIPTION                    | UNITS               | DATA TYPE |
|:-------------:|:--------------|:-------------------------------|:--------------------|:----------|
|       1       | `DATETIME`    | The date/time                  | yyyy-mm-dd hh:mm    | string    |
|       2       | `FLOW`        | Boundary condition flow        | cubic meters/second | float     |
|       3       | `TEMPERATURE` | Boundary condition temperature | degrees Celsius     | float     |

Model run requirements:

| COLUMN NUMBER | COLUMN NAME   | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:--------------|:----------:|:--------------:|:----------------:|
|       1       | `DATETIME`    |  Optional  |    Required    |     Required     |
|       2       | `FLOW`        |  Optional  |    Required    |     Required     |
|       3       | `TEMPERATURE` |  Optional  |    Required    |     Required     |


### 6.8 LAND COVER CODES FILE  
File name: `UserDefinedFileName.[xlsx|csv]`

xlsx sheet name: `Land Cover Codes`

The land cover codes file contains the physical attribute information 
associated with each land cover code. Land cover codes can be alphanumeric 
values. Zero should be avoided as a land cover code. The physical 
attribute such as height, canopy cover, LAI, or overhang must be numeric.
There cannot be skipped rows 
(i.e. rows without information in between rows with information) 
because the model routines see a blank row as the end of the data sequence.

#### 6.8.1 Canopy Type

Land cover canopy information can be input as either canopy cover or 
effective leaf area index. This option is specified in the control file using the key ```canopy_data```.

##### 6.8.1.1 Canopy Cover
Input file formatting when ```canopy_data = "CanopyCover"``` in the control file.

Field details:

| COLUMN NUMBER | COLUMN NAME    | DESCRIPTION       | UNITS                  | DATA TYPE |
|:-------------:|:---------------|:------------------|:-----------------------|:----------|
|       1       | `NAME`         | Land cover Name   | N/A                    |  string   |
|       2       | `CODE`         | Land cover code   | N/A                    |  string   |
|       3       | `HEIGHT`       | Land cover height | meters                 |   float   |
|       4       | `CANOPY`       | Canopy cover      | decimal fraction (0-1) |   float   |
|       5       | `OVERHANG`     | Overhang          | meters                 |   float   |
|       6       | `CANOPY_DEPTH` | Canopy depth      | meters                 |   float   |

Model run requirements:

| COLUMN NUMBER | COLUMN NAME    | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:---------------|:----------:|:--------------:|:----------------:|
|       1       | `NAME`         |  Optional  |    Optional    |     Optional     |
|       2       | `CODE`         |  Required  |    Optional    |     Required     |
|       3       | `HEIGHT`       |  Required  |    Optional    |     Required     |
|       4       | `CANOPY`       |  Required  |    Optional    |     Required     |
|       5       | `OVERHANG`     |  Required  |    Optional    |     Required     |
|       6       | `CANOPY_DEPTH` |  Required  |    Optional    |     Required     |


- Note `CANOPY_DEPTH` is a new column in heat source 9. It was not present in the landcover codes files 
for previous versions of heat source models. If you are updating the model files and lack canopy depth data, 
using the vegetation `HEIGHT` as the canopy depth may be a reasonable approximation. 
When the control file key heatsource8=True, the solar flux and shading methods revert to heat source 8 
methods and do not use the values in the `CANOPY_DEPTH` column. However, values in 
these columns are still required in the current version. See discussion on model updates 
in the documentation for further details. 


##### 6.8.1.2 LAI
Input file formatting when ```canopy_data = "LAI"``` in the control file.

Field details:

| COLUMN NUMBER | COLUMN NAME | DESCRIPTION               | UNITS         | DATA TYPE |
|:-------------:|:------------|:--------------------------|:--------------|:----------|
|       1       | `NAME`         | Land cover Name           | N/A           |  string   |
|       2       | `CODE`         | Land cover code           | N/A           |  string   |
|       3       | `HEIGHT`       | Land cover height         | meters        |   float   |
|       4       | `LAI`          | Effective Leaf Area Index | dimensionless |   float   |
|       5       | `k`            | k extinction coefficient  | dimensionless |   float   |
|       6       | `OVERHANG`     | Overhang                  | meters        |   float   |
|       7       | `CANOPY_DEPTH` | Canopy depth              | meters        |   float   |

Model run requirements:

| COLUMN NUMBER | COLUMN NAME    | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:---------------|:----------:|:--------------:|:----------------:|
|       1       | `NAME`         |  Optional  |    Optional    |     Optional     |
|       2       | `CODE`         |  Required  |    Optional    |     Required     |
|       3       | `HEIGHT`       |  Required  |    Optional    |     Required     |
|       4       | `LAI`          |  Required  |    Optional    |     Required     |
|       5       | `k`            |  Required  |    Optional    |     Required     |
|       6       | `OVERHANG`     |  Required  |    Optional    |     Required     |
|       7       | `CANOPY_DEPTH` |  Required  |    Optional    |     Required     |

- Note `CANOPY_DEPTH` is a new column in heat source 9. It was not present in the landcover codes files 
for previous versions of heat source models. If you are updating the model files and lack canopy depth data, 
using the vegetation `HEIGHT` as the canopy depth may be a reasonable approximation. See discussion on 
model updates in the documentation for further details.

### 6.9 LAND COVER DATA  
File name: `UserDefinedFileName.[xlsx|csv]`

xlsx sheet name: `Land Cover Data` (formally called TTools in Heat Source 8) 

This file defines land cover information. This data can be derived 
from geospatial data using TTools.

Field details:

| COLUMN NUMBER | COLUMN NAME | DESCRIPTION                                 | UNITS           | DATA TYPE |
|:-------------:|:------------|:--------------------------------------------|:----------------|:----------|
|       1       | `STREAM_ID` | Stream ID                                   | N/A             | string    |
|       2       | `NODE_ID`   | Node ID                                     | N/A             | integer   |
|       3       | `STREAM_KM` | Stream km                                   | kilometer       | float     |
|       4       | `LONGITUDE` | Node Longitude                              | decimal degrees | float     |
|       5       | `LATITUDE`  | Node Latitude                               | decimal degrees | float     |
|       6       | `TOPO_W`    | Topographic shade angle to the west         | degrees         | float     |
|       7       | `TOPO_S`    | Topographic shade angle to the south        | degrees         | float     |
|       8       | `TOPO_E`    | Topographic shade angle to the east         | degrees         | float     |
|   multiple    | `LC_T#_S#`  | Land cover code on transect T for sample S  | N/A             | string    |
|   multiple    | `ELE_T#_S#` | Ground elevation on transect T for sample S | meters          | float     |

Model run requirements:

| COLUMN NUMBER | COLUMN NAME  | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:-------------|:----------:|:--------------:|:----------------:|
|       1       | `STREAM_ID`  |  Optional  |    Optional    |     Optional     |
|       2       | `NODE_ID`    |  Required  |    Optional    |     Required     |
|       3       | `STREAM_KM`  |  Required  |    Optional    |     Required     |
|       4       | `LONGITUDE`  |  Required  |    Optional    |     Required     |
|       5       | `LATITUDE`   |  Required  |    Optional    |     Required     |
|       6       | `TOPO_W`     |  Required  |    Optional    |     Required     |
|       7       | `TOPO_S`     |  Required  |    Optional    |     Required     |
|       8       | `TOPO_E`     |  Required  |    Optional    |     Required     |
|   multiple    | `LC_T#_S#`   |  Required  |    Optional    |     Required     |
|   multiple    | `ELE_T#_S#`  |  Required  |    Optional    |     Required     |

The number of `LC_T#_S#` and `ELE_T#_S#` columns is dependent on the number of land cover transects and samples specified 
in the control file. The '#' in the column name will be a number and refers to the specific transect (T) number or 
sample (S) number. For example LC_T2_S4 refers to the land cover sample on transect 2, sample number 4.

The values that go into the `LC_T#_S#` columns are the land cover codes representing land cover attributes at each node 
and sample location. The values that go into the `ELE_T#_S#` columns are the ground elevation at each sample. Each land 
cover code represents a unique combination of vegetation height, canopy cover or LAI, overhang, and canopy depth. The 
attributes for each code are defined in the Land Cover Codes file in long format.

In previous versions of Heat Source, a `values` option was supported for setting up the land cover data file where
the height, density, and overhang at each sample were input. In heat source 9 beta this was set using 
`lcdatainput = values`. This option is no longer supported and all land cover attributes must be setup using the land cover
codes file. Existing value style land cover data files may be converted using a built-in helper function as shown below.

```python
from heatsource9.setup import lc_values_to_codes

lcdatafile_vals = "path/to/lcdatafile_using_values.csv"
lcdatafile_codes_new = "path/to/new_lcdatafile_using_codes.csv"
lccodesfile_new = "path/to/new_lccodesfile_using_codes.csv"

# Convert an old values style land cover data file to land cover codes.
lc_values_to_codes(lcdatafile = lcdatafile_vals,
                   output_lccodefile = lccodesfile_new,
                   output_lcdatafile = lcdatafile_codes_new,
                   canopy_data = "CanopyCover",
                   trans_count = 8,
                   transsample_count = 4,
                   heatsource8 = False,
                   overwrite = False)
```
### 6.10 MORPHOLOGY DATA FILE  
File name: `UserDefinedFileName.[xlsx|csv]`

xlsx sheet name: `Morphology Data`

This file defines channel morphology and substrate information.
Refer to the user manual for more information about each parameter.

Field details:

| COLUMN NUMBER | COLUMN NAME                | DESCRIPTION                   | UNITS                        | DATA TYPE |
|--------------:|:---------------------------|:------------------------------|:-----------------------------|:----------|
|             1 | `STREAM_ID`                | Stream ID                     | N/A                          | string    |
|             2 | `NODE_ID`                  | Node ID                       | N/A                          | integer   |
|             3 | `STREAM_KM`                | Stream km                     | kilometers                   | float     |
|             4 | `ELEVATION`                | Stream Elevation              | meters                       | float     |
|             5 | `GRADIENT`                 | Stream Gradient               | meters/meters                | float     |
|             6 | `BOTTOM_WIDTH`             | Bottom Width                  | meters                       | float     |
|             7 | `CHANNEL_ANGLE_Z`          | Channel Angle z               | meters/meters                | float     |
|             8 | `MANNINGS_n`               | Manning's n                   | seconds/meter                | float     |
|             9 | `SED_THERMAL_CONDUCTIVITY` | Sediment Thermal Conductivity | watts/meters/degrees Celsius | float     |
|            10 | `SED_THERMAL_DIFFUSIVITY`  | Sediment Thermal Diffusivity  | square centimeters/second    | float     |
|            11 | `SED_HYPORHEIC_THICKNESS`  | Hyporheic Zone Thickness      | meters                       | float     |
|            12 | `HYPORHEIC_PERCENT`        | Percent Hyporheic Exchange    | decimal fraction (0-1)       | float     |
|            13 | `POROSITY`                 | Porosity                      | decimal fraction (0-1)       | float     |

Model run requirements:

| COLUMN NUMBER | COLUMN NAME                | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-------------:|:---------------------------|:----------:|:--------------:|:----------------:|
|       1       | `STREAM_ID`                |  Optional  |    Optional    |     Optional     |
|       2       | `NODE_ID`                  |  Required  |    Required    |     Required     |
|       3       | `STREAM_KM`                |  Required  |    Required    |     Required     |
|       4       | `ELEVATION`                |  Required  |    Required    |     Required     |
|       5       | `GRADIENT`                 |  Optional  |    Required    |     Required     |
|       6       | `BOTTOM_WIDTH`             |  Optional  |    Required    |     Required     |
|       7       | `CHANNEL_ANGLE_Z`          |  Optional  |    Required    |     Required     |
|       8       | `MANNINGS_n`               |  Optional  |    Required    |     Required     |
|       9       | `SED_THERMAL_CONDUCTIVITY` |  Optional  |    Optional    |     Required     |
|      10       | `SED_THERMAL_DIFFUSIVITY`  |  Optional  |    Optional    |     Required     |
|      11       | `SED_HYPORHEIC_THICKNESS`  |  Optional  |    Optional    |     Required     |
|      12       | `HYPORHEIC_PERCENT`        |  Optional  |    Optional    |     Optional     |
|      13       | `POROSITY`                 |  Optional  |    Optional    |     Required     |

## 7.0 MODEL OUTPUT FILES
Model outputs are written as CSV files. The first six rows of every output file 
provide the following details:

- `File Created`: Date/time the output file was created.
- `Heat Source Version`: The heat source version number
- `Simulation Name`: User entered text from the `name` field in the control file.
- `User Text`:  User entered text from the `usertxt` field in the control file.
- `Output file description`
- `blank space

The headers for each output file are in the seventh row. For most files, the first header column is
`Datetime`. The output datetimes are formatted as numeric Excel serial date and time values (or the OLE Automation date). 
Specifically, the datetime values are days since December 30, 1899 at 00:00:00. Starting in the second column, all 
the headers are the output stream km. 

Here's a snippet of a stream temperature output file: `Temp_H2O.csv`.
```CSV
File Created:,Tue Mar  3 14:18:18 2026
Heat Source Version:,9.0.0
Simulation Name:,Example River - HS9_example_model_xslx
User Text:,This is an example model using xlsx files.
Output:,Stream Temperature (Celsius)
""
Datetime,10.100,10.050,10.000,9.950,...
37078.0000000,20.1000,19.9912,19.9103,19.8013,...
37078.0416667,19.6000,19.5326,19.4866,19.4099,...
37078.0833333,19.1000,19.0583,19.0309,18.9761,...

```
The table below provides a summary of all the model outputs.

Output Files:

| OUTPUT NAME     | DESCRIPTION                                       | UNITS                  |
|:----------------|:--------------------------------------------------|:-----------------------|
| `Heat_SR1.csv`  | Solar Radiation Flux above Topographic Features   | watts/square meter     |
| `Heat_SR2.csv`  | Solar Radiation Flux below Topographic Features   | watts/square meter     |
| `Heat_SR3.csv`  | Solar Radiation Flux below Land Cover             | watts/square meter     |
| `Heat_SR3b.csv` | Solar Radiation Flux blocked by Land Cover Sample | watts/square meter     |
| `Heat_SR4.csv`  | Solar Radiation Flux below Bank Shade & Emergent  | watts/square meter     |
| `Heat_SR5.csv`  | Solar Radiation Flux Entering Stream              | watts/square meter     |
| `Heat_SR6.csv`  | Solar Radiation Flux Received by Stream           | watts/square meter     |
| `Heat_SR7.csv`  | Solar Radiation Flux Received by Substrate        | watts/square meter     |
| `Heat_Cond.csv` | Streambed Conduction Flux                         | watts/square meter     |
| `Heat_Long.csv` | Longwave Flux                                     | watts/square meter     |
| `Heat_Conv.csv` | Convection Flux                                   | watts/square meter     |
| `Heat_Evap.csv` | Evaporation Flux                                  | watts/square meter     |
| `Rate_Evap.csv` | Evaporation Rate                                  | mm/hour                |
| `Hyd_Disp.csv`  | Hydraulic Dispersion                              | mm/hour                |
| `Hyd_DA.csv`    | Average Depth                                     | meters                 |
| `Hyd_DM.csv`    | Max Depth                                         | meters                 |
| `Hyd_Flow.csv`  | Flow Rate                                         | meters                 |
| `Hyd_Hyp.csv`   | Hyporheic Exchange                                | cubic meters/second    |
| `Hyd_Vel.csv`   | Flow Velocity                                     | meters/second          |
| `Hyd_WT.csv`    | Top Width                                         | square meters/second   |
| `Temp_H2O.csv`  | Stream Temperature                                | Celsius                |
| `Temp_Sed.csv`  | Sediment Temperature                              | Celsius                |
| `Temp_Hyp.csv`  | Hyporheic Return Water Temperature                | Celsius                |
| `Shade.csv`     | Effective Shade                                   | decimal fraction (0-1) |
| `VTS.csv`       | View to Sky                                       | decimal fraction (0-1) |

The table below summarizes the outputs written by model run type. 
YES means the output file is written for that model run type.

| OUTPUT NAME      | SOLAR RUNS | HYDRAULIC RUNS | TEMPERATURE RUNS |
|:-----------------|:----------:|:--------------:|:----------------:|
| `Heat_SR1.csv`   |    YES     |       NO       |       YES        |
| `Heat_SR2.csv`   |    YES     |       NO       |       YES        |
| `Heat_SR3.csv`   |    YES     |       NO       |       YES        |
| `Heat_SR3b.csv`  |    YES     |       NO       |       YES        |
| `Heat_SR4.csv`   |    YES     |       NO       |       YES        |
| `Heat_SR5.csv`   |    YES     |       NO       |       YES        |
| `Heat_SR6.csv`   |     NO     |       NO       |       YES        |
| `Heat_SR7.csv`   |     NO     |       NO       |       YES        |
| `Heat_Cond.csv`  |     NO     |       NO       |       YES        |
| `Heat_Long.csv`  |     NO     |       NO       |       YES        |
| `Heat_Conv.csv`  |     NO     |       NO       |       YES        |
| `Heat_Evap.csv`  |     NO     |       NO       |       YES        |
| `Rate_Evap.csv`  |     NO     |       NO       |       YES        |
| `Hyd_Disp.csv`   |     NO     |       NO       |       YES        |
| `Hyd_DA.csv`     |     NO     |       YES      |       YES        |
| `Hyd_DM.csv`     |     NO     |       YES      |       YES        |
| `Hyd_Flow.csv`   |     NO     |       YES      |       YES        |
| `Hyd_Hyp.csv`    |     NO     |       YES      |       YES        |
| `Hyd_Vel.csv`    |     NO     |       YES      |       YES        |
| `Hyd_WT.csv`     |     NO     |       YES      |       YES        |
| `Temp_H2O.csv`   |     NO     |       NO       |       YES        |
| `Temp_Sed.csv`   |     NO     |       NO       |       YES        |
| `Temp_Hyp.csv`   |     NO     |       NO       |       YES        |
| `Shade.csv`      |    YES     |       NO       |       YES        |
| `VTS.csv`        |    YES     |       NO       |       YES        |

## 8.0 LICENSE
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

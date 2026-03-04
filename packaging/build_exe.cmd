@echo OFF
:: This script builds standalone heat source executables for windows
:: Requires pyinstaller

setlocal
set "SCRIPT_DIR=%~dp0"
pushd "%SCRIPT_DIR%\.." >nul

:: Set the full path to the folder that has all the pyinstaller spec files.
set "PYI_DIR=%CD%\packaging\pyinstaller"

:: Set the folder where the .exe will be saved.
set "DIST_DIR=%CD%\dist\_dev"

:: Set the build directory.
set "BUILD_DIR=%CD%\build\pyinstaller"

:: Set the folder where the pyinstaller

py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --workpath "%BUILD_DIR%" --distpath "%DIST_DIR%" --specpath "%PYI_DIR%" "%PYI_DIR%"\hs9_run_hydraulics.py --hidden-import=heatsource9.Stream.PyHeatsource
py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --workpath "%BUILD_DIR%" --distpath "%DIST_DIR%" --specpath "%PYI_DIR%" "%PYI_DIR%"\hs9_run_solar.py --hidden-import=heatsource9.Stream.PyHeatsource
py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --workpath "%BUILD_DIR%" --distpath "%DIST_DIR%" --specpath "%PYI_DIR%" "%PYI_DIR%"\hs9_run_temperature.py --hidden-import=heatsource9.Stream.PyHeatsource
py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --workpath "%BUILD_DIR%" --distpath "%DIST_DIR%" --specpath "%PYI_DIR%" "%PYI_DIR%"\hs9_setup_control_file.py --hidden-import=heatsource9.Stream.PyHeatsource
py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --workpath "%BUILD_DIR%" --distpath "%DIST_DIR%" --specpath "%PYI_DIR%" "%PYI_DIR%"\hs9_setup_model_inputs.py --hidden-import=heatsource9.Stream.PyHeatsource

echo .
echo .
echo 'dist' folder contains executables
pause

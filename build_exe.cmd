@echo OFF

:: This script builds standalone heat source executables for windows
:: Requires pyinstaller

py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --specpath exe\ exe\hs9_run_hydraulics.py --hidden-import=heatsource9.Stream.PyHeatsource
py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --specpath exe\ exe\hs9_run_solar.py --hidden-import=heatsource9.Stream.PyHeatsource
py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --specpath exe\ exe\hs9_run_temperature.py --hidden-import=heatsource9.Stream.PyHeatsource
py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --specpath exe\ exe\hs9_setup_control_file.py --hidden-import=heatsource9.Stream.PyHeatsource
py -3 -m PyInstaller -F --version-file exe_metadata.txt -i hslogo256.ico --specpath exe\ exe\hs9_setup_model_inputs.py --hidden-import=heatsource9.Stream.PyHeatsource

echo .
echo .
echo 'dist' folder contains executables
pause

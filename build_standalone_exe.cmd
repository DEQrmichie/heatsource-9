@echo OFF

:: This script builds heat source executables for windows 
:: Requires pyinstaller

py -2.7 setup.py build_ext --inplace

py -2.7 -m cython src\heatsource9\Stream\PyHeatsource.pyx --embed
py -2.7 -m cython src\heatsource9\Stream\StreamNode.pyx --embed

py -2.7 -m PyInstaller -F exe\hs9_run_hydraulics.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
py -2.7 -m PyInstaller -F exe\hs9_run_solar.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
py -2.7 -m PyInstaller -F exe\hs9_run_temperature.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
py -2.7 -m PyInstaller -F exe\hs9_setup_control_file.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
py -2.7 -m PyInstaller -F exe\hs9_setup_model_inputs.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
echo .
echo .
echo 'dist' folder contains executables
pause

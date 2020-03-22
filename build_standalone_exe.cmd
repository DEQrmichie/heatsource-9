@echo OFF

:: This script builds heat source exectubles for windows 
:: Requires pyinstaller

C:\Python27\python.exe setup.py build_ext --inplace

cython src\heatsource9\Stream\PyHeatsource.pyx --embed
cython src\heatsource9\Stream\StreamNode.pyx --embed

pyinstaller -F exe\hs9_run_hydraulics.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
pyinstaller -F exe\hs9_run_solar.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
pyinstaller -F exe\hs9_run_temperature.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
pyinstaller -F exe\hs9_setup_control_file.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
pyinstaller -F exe\hs9_setup_model_inputs.py --hidden-import=heatsource9.Stream.PyHeatsource --specpath exe\
echo .
echo .
echo 'dist' folder contains exectubles
pause

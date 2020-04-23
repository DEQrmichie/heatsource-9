@echo OFF

:: This script builds heat source exectubles for windows 
:: Requires pyinstaller
C:\Python27\python.exe setup.py build_ext --inplace

cython heatsource9\Stream\PyHeatsource.pyx --embed
cython heatsource9\Stream\StreamNode.pyx --embed

pyinstaller -F HS9_Run_Hydraulics_Only.py --hidden-import=heatsource9.Stream.PyHeatsource --icon="HSLogo16.ico"
pyinstaller -F HS9_Run_Solar_Only.py --hidden-import=heatsource9.Stream.PyHeatsource --icon="HSLogo16.ico"
pyinstaller -F HS9_Run_Temperature.py --hidden-import=heatsource9.Stream.PyHeatsource --icon="HSLogo16.ico"
pyinstaller -F HS9_Setup_Control_File.py --hidden-import=heatsource9.Stream.PyHeatsource --icon="HSLogo16.ico"
pyinstaller -F HS9_Setup_Model_Inputs.py --hidden-import=heatsource9.Stream.PyHeatsource --icon="HSLogo16.ico"
echo .
echo .
echo 'dist' folder contains exectubles
pause

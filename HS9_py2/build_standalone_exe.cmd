@echo OFF

:: This script builds heat source exectubles for windows 
:: Requires pyinstaller
pyinstaller -F HS9_Run_Hydraulics_Only.py
pyinstaller -F HS9_Run_Solar_Only.py
pyinstaller -F HS9_Run_Temperature.py
pyinstaller -F HS9_Setup_Control_File.py
pyinstaller -F HS9_Setup_Model_Inputs.py
echo .
echo .
echo 'dist' folder contains exectubles
pause

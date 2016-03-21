@echo OFF

:: This script builds an install exectuble for windows and 
:: then deletes the build folder
python setup.py bdist_wininst --target-version="2.7"
rd /S /q build
echo .
echo .
echo New 'dist' folder contains Heat Source 9 install exectuble
echo Use this executbale to complete installation
pause

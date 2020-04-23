@echo OFF

:: This script builds a sdist and wheel for windows and 
:: then deletes the build folder

py -3 setup.py sdist bdist_wheel

rd /S /q build
echo .
echo .
echo 'dist' folder contains Heat Source 9 sdist and wheel
pause

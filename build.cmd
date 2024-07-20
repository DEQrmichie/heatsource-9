@echo OFF

:: This script builds a sdist and wheel for windows
py -3 -m build --sdist
py -3 -m build --wheel

rd /S /q build
echo .
echo .
echo 'dist' folder contains Heat Source 9 sdist and wheel
pause

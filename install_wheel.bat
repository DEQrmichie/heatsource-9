@echo OFF

:: This script installs the heatsource wheel for windows in the local directory
py -m pip install --upgrade --force-reinstall dist/heatsource9-9.0.0b27-cp312-cp312-win32.whl --user

:: This script installs the heatsource wheel for windows in the global directory
py -m pip install --upgrade --force-reinstall dist/heatsource9-9.0.0b27-cp312-cp312-win32.whl
pause

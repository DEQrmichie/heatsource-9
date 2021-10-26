@echo OFF

:: This script installs the heatsource wheel for windows in the local directory
pip3 install --upgrade --force-reinstall dist/heatsource9-9.0.0b26-cp38-cp38-win32.whl --user

:: This script installs the heatsource wheel for windows in the global directory
pip3 install --upgrade --force-reinstall dist/heatsource9-9.0.0b26-cp38-cp38-win32.whl
pause

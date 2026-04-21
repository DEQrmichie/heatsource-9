@echo OFF
setlocal
set "SCRIPT_DIR=%~dp0"
pushd "%SCRIPT_DIR%\.." >nul

set "WHL_DIR=%CD%\dist\_dev\"
set "WHL=heatsource9-9.0.0b32-cp312-cp312-win_amd64.whl"
set "WHL_PATH=%WHL_DIR%%WHL%"

:: This script installs the heatsource wheel for windows in the local directory
py -m pip install --force-reinstall "%WHL_PATH%"
popd >nul
endlocal
pause

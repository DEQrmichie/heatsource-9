@echo OFF
:: This script builds a wheel for windows

setlocal
set "SCRIPT_DIR=%~dp0"
pushd "%SCRIPT_DIR%\.." >nul

:: Set the folder where the wheel will be saved.
set "DIST_DIR=%CD%\dist\_dev"

py -3 -m build --outdir "%DIST_DIR%" --wheel

rd /S /q build
echo .
echo .
echo "%DIST_DIR%" folder contains Heat Source 9 wheel
popd >nul
endlocal
pause

@echo OFF
:: This script builds a wheel for windows and recursively scans the project repo and deletes build artifact files
:: anything in dist\ or .venv\ is kept.
:: This script only works by double clicking the .cmd files.
:: If run from cmd prompt the double percent signs in the for loops %% need to be single %

setlocal
set "SCRIPT_DIR=%~dp0"
pushd "%SCRIPT_DIR%\.." >nul

:: Set the folder where the wheel will be saved.
set "DIST_DIR=%CD%\dist\_dev"

py -3 -m build --outdir "%DIST_DIR%" --wheel

:: Delete root level artifacts, keep dist and .venv.
echo Removing build artifacts... one minute

for %%D in (
  build
  src\heatsource9.egg-info
  .pytest_cache
  pip-wheel-metadata
) do (
  if exist "%%~D" rd /S /q "%%~D"
)

:: Delete compiled files outside of dist and .venv
for %%E in (pyc pyo pyd) do (
  for /r %%F in (*.%%E) do (
    echo(%%~fF| findstr /i /c:"\dist\" /c:"\.venv\" >nul
    if errorlevel 1 if exist "%%~fF" del /F /Q "%%~fF"
  )
)

:: Delete __pycache__ directories outside of dist and .venv
for /d /r %%D in (__pycache__) do (
  echo(%%~fD| findstr /i /c:"\dist\" /c:"\.venv\" >nul
  if errorlevel 1 if exist "%%~fD" rd /S /Q "%%~fD"
)

echo .
echo .
echo Done. "%DIST_DIR%" folder contains Heat Source 9 wheel
popd >nul
endlocal
pause

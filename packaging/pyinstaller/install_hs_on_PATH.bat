@echo OFF
:: This installs hs.exe from this folder into the versioned folder at
:: C:\Users\username\AppData\Local\Programs\HeatSource
:: then copies the executable into C:\Users\username\AppData\Local\Programs\HeatSource\current
:: then adds the 'current' folder to the user's PATH by calling path.ps1.

setlocal
set "SCRIPT_DIR=%~dp0"
set "PS_SCRIPT=%SCRIPT_DIR%path.ps1"
set "HS_EXE=%SCRIPT_DIR%hs.exe"

if not exist "%PS_SCRIPT%" (
    echo Could not find %PS_SCRIPT%
    pause
    exit /b 1
)

powershell -NoProfile -ExecutionPolicy Bypass -File "%PS_SCRIPT%" -SourceExePath "%HS_EXE%"
if errorlevel 1 (
    echo.
    echo Installation failed.
    pause
    exit /b 1
)

echo.
echo Installation complete.
pause
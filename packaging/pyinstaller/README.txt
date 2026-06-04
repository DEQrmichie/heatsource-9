Heat Source 9 Windows install

Install steps

1. Extract the content of the zip file to a folder. Do not run files from inside the zip.
   If you are reading this you probably already completed this step.
2. In the extracted folder, double click install_hs_on_PATH.bat.
3. Wait for the installer window to finish.
4. Open a new Command Prompt, PowerShell, or Windows Terminal window.
5. Run:

   hs -v

If Windows recognized the hs command and prints the model version number, it worked.

Typical use

1. Change to your model directory. e.g. cd C:\path\to\model_directory
2. Run a temperature model using the command:

   hs run -t

Notes

- install_hs_on_PATH.bat copies hs.exe into the versioned folder at C:\Users\username\AppData\Local\Programs\HeatSource
  then copies the executable into C:\Users\username\AppData\Local\Programs\HeatSource\current
  then adds the 'current' folder to the user's PATH by calling path.ps1.
- You may need to close and reopen your shell before the `hs` command is recognized.
- The other hs9_*.exe files can still be used directly if needed. They must be in the same folder as the control file.
- use `hs setup -h` for setup command help and `hs run -h` for model run command help.
- See the documentation on GitHub or the user manual for more info. https://github.com/DEQrmichie/heatsource-9
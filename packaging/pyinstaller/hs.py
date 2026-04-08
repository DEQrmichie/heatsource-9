#!/usr/bin/python3

"""Heat Source entry point executable.
This script is used by PyInstaller to build a Windows executable which can be used to for entry point commands.

e.g hs run -t.

The resulting hs.exe must be included on PATH for command line to work.
"""

from heatsource9.__main__ import main


if __name__ == "__main__":
    raise SystemExit(main())

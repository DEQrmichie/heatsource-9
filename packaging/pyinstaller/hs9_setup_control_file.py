#!/usr/bin/python3

"""Run Heat Source control file setup from the local model directory.

This script is used by PyInstaller to build a Windows executable and can
also be run directly with Python. It determines the model directory from:
1) the executable location when running as an executable, or
2) this script's location when running as a normal Python script.

It then writes a blank `HeatSource_Control.xlsx` file in that directory.
"""
import sys
from pathlib import Path

from heatsource9.setup import setup_cf


def main():
    # Use executable folder for PyInstaller builds, otherwise script folder.
    if getattr(sys, "frozen", False):
        model_dir = Path(sys.executable).resolve().parent
    else:
        model_dir = Path(__file__).resolve().parent

    # Write a blank control file in XLSX format.
    setup_cf(
        model_dir=model_dir,
        control_file="HeatSource_Control.xlsx",
        use_timestamp=False,
        overwrite=False,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

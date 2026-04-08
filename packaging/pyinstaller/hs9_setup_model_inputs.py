#!/usr/bin/python3

"""Run Heat Source model input setup from the local model directory.

This script is used by PyInstaller to build a Windows executable and can
also be run directly with Python. It determines the model directory from:
1) the executable location when running as an executable, or
2) this script's location when running as a normal Python script.

It then finds that `HeatSource_Control.[xlsx|csv]` is in that directory and
writes blank model input files based on that control file.
"""
import sys
from pathlib import Path

from heatsource9.io.control_file import cf_path
from heatsource9.setup import write_mi


def main():
    # Use executable folder for PyInstaller builds, otherwise script folder.
    if getattr(sys, "frozen", False):
        model_dir = Path(sys.executable).resolve().parent
    else:
        model_dir = Path(__file__).resolve().parent

    # Find and validate HeatSource_Control.xlsx or HeatSource_Control.csv.
    control_path = cf_path(model_dir)

    # Write blank model input files from the selected control file.
    write_mi(
        model_dir=model_dir,
        control_file=str(control_path),
        use_timestamp=False,
        overwrite=False,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

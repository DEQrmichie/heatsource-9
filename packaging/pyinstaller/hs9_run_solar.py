#!/usr/bin/python3

"""Run the Heat Source solar model from the local model directory.

This script is used by PyInstaller to build a Windows executable and can
also be run directly with Python. It determines the model directory from:
1) the executable location when running as an executable, or
2) this script's location when running as a normal Python script.

It then finds that `HeatSource_Control.[xlsx|csv]` is in that directory and
runs the solar model.
"""
import sys
from pathlib import Path

import heatsource9.run as hs_run


def main():
    # Use executable folder for PyInstaller builds, otherwise script folder.
    if getattr(sys, "frozen", False):
        model_dir = Path(sys.executable).resolve().parent
    else:
        model_dir = Path(__file__).resolve().parent

    # Run Heat Source Solar.
    hs_run.solar(model_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

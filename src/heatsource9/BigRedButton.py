"""
Backward compatible functions for older Heat Source 8 and beta 
version Heat Source 9 scripts. Preserves the BRB entrypoints 
(e.g., `setup_cf`, `setup_mi`, `run`, `run_solar`, `run_hydraulics`). 
It passes execution to heatsource9.setup and heatsource9.run.
"""

from pathlib import Path

import heatsource9.run as hs_run
from heatsource9.setup import setup_cf as _setup_cf_brb, write_mi


def setup_cf(
    model_dir,
    control_file = "HeatSource_Control.xlsx",
    use_timestamp = False,
    overwrite = False,
    **kwargs,
):
    """
    Backward compatible entry for control file setup and parameterization.
    """

    return _setup_cf_brb(
        model_dir=model_dir,
        control_file=control_file,
        use_timestamp=use_timestamp,
        overwrite=overwrite,
        strict=True,
        **kwargs)

def setup_mi(model_dir, control_file = "HeatSource_Control.xlsx", 
             use_timestamp = False, overwrite = False,
):
    """Backward compatible entry for blank model input file generation."""

    model_path = Path(model_dir)

    mi = write_mi(model_dir=model_path, control_file=control_file, 
                  use_timestamp=use_timestamp, overwrite=overwrite)
    
    return mi

def run(model_dir, control_file = None):
    """Backward compatible entry for temperature runs."""
    hs_run.temperature(model_dir, control_file)

def run_solar(model_dir, control_file = None):
    """Backward compatible entry for solar runs."""
    hs_run.solar(model_dir, control_file)

def run_hydraulics(model_dir, control_file = None):
    """Backward compatible entry for hydraulic only runs."""
    hs_run.hydraulics(model_dir, control_file)

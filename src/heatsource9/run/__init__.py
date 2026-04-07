from heatsource9.io.control_file import cf_path
from heatsource9.run.model_runner import ModelRunner


def temperature(model_dir, control_file = None):
    ModelRunner(
        control_file_path=cf_path(model_dir, control_file),
    ).temperature()


def solar(model_dir, control_file = None):
    ModelRunner(
        control_file_path=cf_path(model_dir, control_file),
    ).solar()


def hydraulics(model_dir, control_file = None):
    ModelRunner(
        control_file_path=cf_path(model_dir, control_file),
    ).hydraulics()


__all__ = [
    "temperature",
    "solar",
    "hydraulics",
    "ModelRunner",
]

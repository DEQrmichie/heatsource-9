import logging
from pathlib import Path

from heatsource9.io.control_file import cf_path
from heatsource9.io.logging_config import configure_logging
from heatsource9.run.model_runner import ModelRunner

logger = logging.getLogger(__name__)


def temperature(model_dir, control_file = None):
    model_path = Path(model_dir).expanduser().resolve()
    try:
        configure_logging(model_dir=model_path, overwrite=True)
        control_path = cf_path(model_path, control_file)
    except Exception as exc:
        msg = str(exc)
        logger.exception(msg)
        raise
    ModelRunner(control_file_path=control_path).temperature()


def solar(model_dir, control_file = None):
    model_path = Path(model_dir).expanduser().resolve()
    try:
        configure_logging(model_dir=model_path, overwrite=True)
        control_path = cf_path(model_path, control_file)
    except Exception as exc:
        msg = str(exc)
        logger.exception(msg)
        raise
    ModelRunner(control_file_path=control_path).solar()


def hydraulics(model_dir, control_file = None):
    model_path = Path(model_dir).expanduser().resolve()
    try:
        configure_logging(model_dir=model_path, overwrite=True)
        control_path = cf_path(model_path, control_file)
    except Exception as exc:
        msg = str(exc)
        logger.exception(msg)
        raise
    ModelRunner(control_file_path=control_path).hydraulics()


__all__ = [
    "temperature",
    "solar",
    "hydraulics",
    "ModelRunner",
]

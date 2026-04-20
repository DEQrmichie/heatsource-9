
import logging
from pathlib import Path

PROGRESS_LOGGER_NAME = "heatsource9.progress"


class _ExactLevelFilter(logging.Filter):
    def __init__(self, levelno):
        super().__init__()
        self.levelno = int(levelno)

    def filter(self, record):
        return record.levelno == self.levelno


class _MinimumLevelFilter(logging.Filter):
    def __init__(self, levelno):
        super().__init__()
        self.levelno = int(levelno)

    def filter(self, record):
        return record.levelno >= self.levelno


def configure_logging(model_dir, *, overwrite = True):
    """
    Configure logging for console and log files.

    Console output shows status, progress, warning, and error messages.

    File logging writes:
    - hs.log for all messages, including status, progress messages,
      warnings, errors, and exceptions
    - hs_warning.log for WARNING only
    - hs_error.log for ERROR and above
    """
    model_path = Path(model_dir).expanduser().resolve()
    log_path = model_path / "hs.log"
    warning_log_path = model_path / "hs_warning.log"
    error_log_path = model_path / "hs_error.log"

    # Reset progress logger before reconfiguring root handlers.
    progress_logger = logging.getLogger(PROGRESS_LOGGER_NAME)
    for handler in list(progress_logger.handlers):
        progress_logger.removeHandler(handler)
        try:
            handler.close()
        except Exception:
            pass

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(logging.Formatter("%(message)s"))

    mode = "w" if overwrite else "a"
    file_handler = logging.FileHandler(log_path, mode=mode, encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(message)s"))

    warning_handler = logging.FileHandler(warning_log_path, mode=mode, encoding="utf-8", delay=True)
    warning_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(message)s"))
    warning_handler.setLevel(logging.WARNING)
    warning_handler.addFilter(_ExactLevelFilter(logging.WARNING))

    error_handler = logging.FileHandler(error_log_path, mode=mode, encoding="utf-8", delay=True)
    error_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(message)s"))
    error_handler.setLevel(logging.ERROR)
    error_handler.addFilter(_MinimumLevelFilter(logging.ERROR))

    logging.basicConfig(
        level=logging.INFO,
        handlers=[stream_handler, file_handler, warning_handler, error_handler],
        force=True,
    )

    # Progress logger writes only to file so console progress bars do not duplicate.
    progress_logger.setLevel(logging.INFO)
    progress_logger.propagate = False
    progress_logger.addHandler(file_handler)

    return {
        "log": log_path,
        "warning_log": warning_log_path,
        "error_log": error_log_path,
    }

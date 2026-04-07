
import logging
from pathlib import Path

PROGRESS_LOGGER_NAME = "heatsource9.progress"


def configure_logging(model_dir, *, overwrite = True):
    """
    Configure logging for console and heatsource.log.

    Console prints progress message text only, file logging 
    track every itteration setp with a timestamp.
    The log file is overwritten by default for each run.
    """
    model_path = Path(model_dir).expanduser().resolve()
    log_path = model_path / "heatsource.log"

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

    logging.basicConfig(
        level=logging.INFO,
        handlers=[stream_handler, file_handler],
        force=True,
    )

    # Progress logger writes only to file so console progress bars do not duplicate.
    progress_logger.setLevel(logging.INFO)
    progress_logger.propagate = False
    progress_logger.addHandler(file_handler)

    return log_path

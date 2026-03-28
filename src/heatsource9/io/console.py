"""Console output and logging tools."""


import logging
import sys

from heatsource9.io.logging_config import PROGRESS_LOGGER_NAME

_LOGGER_NAME = "heatsource9"


class Console:
    """Class to print messages and progress updates to the console and logger."""

    def __init__(self, logger_name = _LOGGER_NAME):
        """
        This sets the logger and initializes model progress output.
        """
        self.log = logging.getLogger(logger_name)
        self.progress_log = logging.getLogger(PROGRESS_LOGGER_NAME)
        self._is_tty = bool(getattr(sys.stdout, "isatty", lambda: False)())

    def info(self, msg):
        self.log.info(str(msg))

    def warning(self, msg):
        self.log.warning(str(msg))

    def error(self, msg):
        self.log.error(str(msg))

    def debug(self, msg):
        self.log.debug(str(msg))

    def progress(self, msg, current, total):
        """
        Output a progress update to the console.

        Progress updates are emitted as carriage return lines so each step
        overwrites the prior console line until completion.
        Progress updates are always written to heatsource.log at every step.
        """
        total_i = int(total) if total not in (None, 0, "0") else 1
        cur_i = int(current) if current is not None else 0
        cur_i = max(0, min(cur_i, total_i))
        line = self.progress_line(str(msg), cur_i, total_i)
        log_line = self.progress_log_line(str(msg), cur_i, total_i)

        if self._is_tty:
            # \x1b[2K clears the full current terminal line so no old characters remain.
            # \r moves the cursor back to column 1 so the next progress text overwrites that line.
            ansi_prefix = "\x1b[2K\r"
            # Keep writing on the same line until the final step, then add \n 
            # to lock the line in place.
            end = "\n" if cur_i >= total_i else ""
            # Print one complete progress line.
            sys.stdout.write(ansi_prefix + line + end)
            # Flush immediately so each progress step appears in real time.
            sys.stdout.flush()
        else:
            sys.stdout.write(line + "\n")
            sys.stdout.flush()

        # Write progress to file without duplicating console formatting.
        self.progress_log.info(log_line)

    @staticmethod
    def progress_line(msg, current, total):
        percent = int(current / total * 100) if total else 0
        prefix = f"{msg} {current} | {total}"
        bar_start = " ["
        bar_end = f"] {percent}%"
        bar_size = max(1, 79 - len(prefix + bar_start + bar_end))
        amount = min(bar_size, int(current / (total / float(bar_size)))) if total else 0
        remain = bar_size - amount
        bar = "." * amount + " " * remain
        return prefix + bar_start + bar + bar_end

    @staticmethod
    def progress_log_line(msg, current, total):
        percent = int(current / total * 100) if total else 0
        return f"{msg} ({current} {total} {percent}%)"

console = Console()

def print_console(msg, progress_bar = False, current = 100, total = 100):
    """Print a message or progress update using the shared Console instance."""
    if progress_bar:
        console.progress(str(msg), current=current, total=total)
    else:
        console.info(str(msg))

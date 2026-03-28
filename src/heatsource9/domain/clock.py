"""
Clock class and time calculation utilities.
"""

from datetime import datetime, timezone
from time import gmtime


class Clock:
    """
    Manages model simulation time progression.

    Internal model time and time progression is measured in
    UNIX epoch seconds. Epoch seconds is the number of seconds 
    since the epoch (1970-01-01 00:00:00 UTC)

    start_time: Start time of the simulation (includes flushtime)
    end_time: End time of the simulation.
    timestep_seconds: Timestep in seconds.

    start_time and end_time can be unix epoch seconds (int/float)
    or a datetime (timezone aware; but if no timezone attached treated as UTC)

    Clock is a one way trip through time. There's no restart ability.
    A restart() method might need to be added if we ever want to 
    rerun the same Clock or the same Simulation object again without 
    manually restarting the model.
    
    Returns UNIX epoch seconds.
    """

    def __init__(self, start_time, end_time, timestep_seconds):
        if timestep_seconds <= 0:
            raise ValueError("timestep_seconds must be > 0")

        self._start = _to_epoch_seconds(start_time)
        self._end = _to_epoch_seconds(end_time)
        self._dt = float(int(timestep_seconds))

        self._current = float(self._start)

    def __iter__(self):
        """
        Yield epoch seconds for each simulation step, inclusive of end_time.
        """
        while self._current <= self._end:
            yield float(self._current)
            self._current += self._dt


def _to_epoch_seconds(t):
    """
    Convert either epoch seconds or datetime to UNIX epoch seconds.
    Epoch seconds is the number of seconds since the epoch (1970-01-01 00:00:00 UTC)
    """
    if isinstance(t, (int, float)):
        return float(t)
    if isinstance(t, datetime):
        dtu = t.astimezone(timezone.utc) if t.tzinfo else t.replace(tzinfo=timezone.utc)
        return float(dtu.timestamp())
    raise TypeError(f"Unsupported time type: {type(t).__name__}")


def time_parts(epoch_seconds):
    """
    Return time values calculated from UNIX epoch seconds.

    (hour, minute, second, JD, JC)

    Where:
      - JD: julian day (tm_yday from gmtime)
      - JC: julian century
    """
    model_time = gmtime(float(epoch_seconds))
    hour = int(model_time.tm_hour)
    minute = int(model_time.tm_min)
    second = int(model_time.tm_sec)
    jd = int(model_time.tm_yday)
    jc = julian_century(float(epoch_seconds))
    return hour, minute, second, jd, jc


def julian_century(epoch_seconds):
    """
    Returns Julian Century (JC) from epoch seconds.
    This same function returned JDC in older models (Julian Day Century).
    JC = JDC. Changed to JC because JDC is not a common term.

    Julian Century is centuries elapsed since Julian Day 2451545.0 (J2000):

        JC = (Julian Day - 2451545.0) / 36525.0

    The formulation follows the equations in Meeus (1988, 1991).

    Here JC is computed on demand for the current time.
    """
    model_time = gmtime(float(epoch_seconds))
    y = int(model_time.tm_year)
    m = int(model_time.tm_mon)
    d = int(model_time.tm_mday)

    if m < 3:
        m += 12
        y -= 1

    julian_day = int(365.25 * (y + 4716.0)) + int(30.6001 * (m + 1)) + d - 1524.5

    if julian_day > 2299160.0:
        a = int(y / 100)
        b = (2 - a + int(a / 4))
        julian_day += b

    return round((julian_day - 2451545.0) / 36525.0, 10)


def pretty_time(epoch_seconds):
    """Return model epoch seconds as a human readable UTC date/time string."""
    return datetime.fromtimestamp(float(epoch_seconds), tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

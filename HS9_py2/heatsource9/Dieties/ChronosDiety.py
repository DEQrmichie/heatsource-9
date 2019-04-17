from __future__ import division
# Heat Source, Copyright (C) 2000-2016, 
# Oregon Department of Environmental Quality

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from builtins import object
from past.utils import old_div
from time import ctime, gmtime

class ChronosDiety(object):
    """This class provides a clock to be used in the model timestepping.

    The ChronosDiety controls the current time and steps through time
    in an organized fashion based on start/stop/step functionality.
    The time is stored as seconds from the epoch. Changes to this
    class to support more advanced timezone and datetime functionality
    (as an older version had) should be pushed into a subclass."""
    def __init__(self):
        self.second = 1
        self.minute = self.second * 60
        self.hour = self.minute * 60
        self.day = self.hour * 24
        self.week = self.day * 7
        self.__current = None # Current time
        self.__start = None
        self.__dt = None
        self.__stop = None

    def __call__(self,tick=False):
        """Return current time as seconds since epoch

        Returns the current time float if called with no arguments. If tick=True,
        advances the clock one timestep and returns new current time. If we are
        in a spin-up period, it tests whether the advance would move us to the
        next day and if so it reverts to the beginning of the spin-up day, thus
        running a single day over and over until the spin-up time is complete."""
        # Quickly return if we're not moving forward in time
        if not tick: return self.__current
        # First we increment the current time by dt
        self.__current += self.__dt
        if gmtime(self.__current)[2] != gmtime(self.__thisday)[2]:
            # Day numbers not equal, need to recalculate julian day
            self.__thisday = self.__current
            self.CalcJulianCentury()
        # Then test whether we're still spinning up

        if self.__current < self.__start:
            # If we're still in the spin-up period
#            print "True",
#            #Make sure we don't advance to next day (i.e. just run the first day over and over)
#            if gmtime(self.__spin_current+self.__dt)[2] != gmtime(self.__spin_start)[2]:
#                self.__spin_current = self.__spin_start # We would have advanced, so we start again on the first day
#            else:
#                print "False"
#                self.__spin_current += self.__dt # We're spinning up and haven't advanced, so use the current spin-up time
#            # Set TheTime according to either spin_current or current and then return it
            self.__spin_current += self.__dt
            return self.__spin_current
        else:
            return self.__current

    def __iter__(self):
        """Iterator support to allow 'for tm in Chronos' loops"""
        if not self.__start or not self.__dt:
            raise Exception("Must call %s with the Start() method before using." % self.__class__.__name__)
        while self.__current <= self.__stop:
            yield self.__current
            self(True)
    def __len__(self):
        """Length will report the number of timesteps

        Currently, this method will screw things up because it actually
        iterates through the sequence, cycling through time. This is a problem."""
        raise NotImplementedError("This needs some work- do we actually need it?")
        #return len([i for i in self])

    def PrettyTime(self): return ctime(self.__current)
    def Year(self): return gmtime(self.__current)[0]
    def Month(self): return gmtime(self.__current)[1]
    def Day(self): return gmtime(self.__current)[2]
    def TimeTuple(self):
        year, month, day, hour, minute, second, weekday, jday, offset = gmtime(self.__current)
        return year, month, day, hour, minute, second, jday, offset, self.__jdc
    #def ExcelTime(self): return float(pyTime(self.__current))

    def Start(self, start, dt=None, stop=None, spin=0, offset=0):
        """Initialize the clock to some default values and get ready to run.

        Initial values are starting and stopping times, timestep in seconds, number
        of days to spin up and the timezone offset in hours."""
        #Make sure the start, stop and dt values are datetime instances
        if (not isinstance(start, float) and not isinstance(start, int)) or \
            (stop and (not isinstance(stop, float) and not isinstance(stop, int))) or \
            (dt and (not isinstance(dt, float) and not isinstance(dt, int))):
            raise Exception("Start, stop times and timestep much be floating point numbers or integers.")
        # Some values used in internal calculations
        self.__offset = offset # The timezone offset, default to GMT
        self.__start = start
        self.__dt = dt or self.minute # There's a default dt of one minute
        self.__stop = stop or self.__start + self.day # There's a default runtime of one day
        self.__spin_start = self.__start- (spin*86400) if spin else self.__start # Start of the spin-up period
        self.__spin_current = self.__spin_start # Current time within the spinup period
        self.__current = self.__spin_current
        # Placeholder for deciding whether we have to recalculate the julian day
        self.__thisday = self.__current-self.__dt
        # Placeholder for making sure we don't leave the spin-up period until the right time
        self.__spinday = gmtime(self.__spin_current)[2]
        self.CalcJulianCentury()

    def CalcJulianCentury(self):
        # Then break out the time into a tuple
        y,m,d,H,M,S,day,wk,tz = gmtime(self.__current)
        dec_day = d + old_div((H + old_div((M + old_div(S,60)),60)),24)

        if m < 3:
            m += 12;
            y -= 1;

        julian_day = int(365.25*(y+4716.0)) + int(30.6001*(m+1)) + d - 1524.5;

        # This value should only be added if we fall after a certain date
        if julian_day > 2299160.0:
            a = int(old_div(y,100))
            b = (2 - a + int(old_div(a,4)))
            julian_day += b
        #This is the julian century
        self.__jdc = round((julian_day-2451545.0)/36525.0,10) # Eqn. 2-5 in HS Manual

    #####################################################
    # Properties to allow reading but no changes
    start = property(lambda self: self.__start)
    stop = property(lambda self: self.__stop)
    dt = property(lambda self: self.__dt)
    offset = property(lambda self: self.__offset)
    TheTime = property(lambda self: self.__current)
    JD = property(lambda self: (self.__jd, self.__jdc))

Chronos = ChronosDiety()

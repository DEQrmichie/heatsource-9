from __future__ import print_function, division
import unittest
from time import ctime, gmtime
import jdcal

from heatsource9.Dieties.ChronosDiety import Chronos

# This test evaluates if using true division will impact HS calculation of julian day.
# I pulled the relevant method out and comparing HS julian day to what is calculated by the jdcal package.
# Test resut = All good. HS method produces same result as jdcal package using true division and int division

def CalcJulianCentury(start):
    # Then break out the time into a tuple
    y, m, d, H, M, S, day, wk, tz = gmtime(start)
    
    print(y, m, d)

    z = jdcal.gcal2jd(year=y, month=m, day=d)

    if m < 3:
        m += 12
        y -= 1

    julian_day = int(365.25 * (y + 4716.0)) + int(30.6001 * (m + 1)) + d - 1524.5

    # This value should only be added if we fall after a certain date
    if julian_day > 2299160.0:
        a = int(y / 100)
        b = (2 - a + int(a / 4))
        julian_day += b


    julian_day2 = z[0] + z[1]
    if julian_day != julian_day2:
        print(julian_day, julian_day2)
    

    # This is the julian century
    jdc = round((julian_day - 2451545.0) / 36525.0, 10)  # Eqn. 2-5 in HS Manual

    return jdc

# Loop through a range of years
# -3786825600 1/1/1850
# 5711731200 12/31/2150
Chronos.Start(start=-3786825600, dt=86400.0,
              stop=5711731200, offset=-7)

# 993945600 7/1/2001
# 994550400 +24 hours
# Chronos.Start(start=993945600, dt=60.0,
#              stop=993945600 + (86400 * 1), offset=-7)

#year, month, day, hour, minute, second, JD, offset, JDC = Chronos.TimeTuple()

#jdcl = []
while Chronos.TheTime <= Chronos.stop:
    jdc = CalcJulianCentury(start=Chronos.TheTime)
#    jdcl.append(jdc)
    Chronos(tick=True)

print("done")
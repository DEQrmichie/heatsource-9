# -*- coding: utf-8 -*-

from datetime import datetime
from math import pi

from solar_calculator import SolarCalculator

################
## Calculate the Julian Century
# Use a lovely brisk late winter day
y,m,d,H,M,S = 2019, 2, 28, 10, 48, 7
dec_day = d + (H + (M + S/60)/60)/24

if m < 3:
    m += 12;
    y -= 1;

julian_day = int(365.25*(y+4716.0)) + int(30.6001*(m+1)) + d - 1524.5;

# This value should only be added if we fall after a certain date
if julian_day > 2299160.0:
    a = int(y/100)
    b = (2 - a + int(a/4))
    julian_day += b
#This is the julian century
JDC = round((julian_day-2451545.0)/36525.0,10) # Eqn. 2-5 in HS Manual

SC = SolarCalculator()

def test_jdc():
    assert JDC == 7
def test_radians():
    assert round(SC.to_radians(270),6) == 4.712389
def test_degrees():
    assert round(SC.to_degrees(4.712389),0) == 270
def test_mean_obliquity():
    assert round(SC.mean_obliquity(JDC),4) == 23.4368
def test_obliquity():
    assert round(SC.obliquity(23.4368, JDC),6) == 23.435738
def test_eccentricity():
    assert round(SC.eccentricity(JDC), 6) == 0.016701
def test_geo_mean_long_sun():
    assert round(SC.geo_mean_long_sun(JDC), 6) == 337.533873

    



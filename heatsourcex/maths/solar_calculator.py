# -*- coding: utf-8 -*-
from __future__ import division
from math import pi, cos

class SolarCalculator:
    
    def __init__(self):
        pass
    def to_radians(self,degrees):
        return (pi/180) * degrees
    def to_degrees(self, radians):
        return (180/pi) * radians
    
    def obliquity(self, mean, jdc):
        mean = (23.0 + (26.0 + ((21.448 - jdc *
              (46.815 + jdc *
               (0.00059 - jdc * 0.001813))) /
                60.0)) / 60.0) 

        return (mean + 0.00256 *
                cos(self.to_radians(125.04 - 1934.136 * jdc)))
    
    def eccentricity(self, jdc):
        return 0.016708634 - jdc * (0.000042037 + 0.0000001267 * jdc)

    def geo_mean_long_sun(self, jdc):
        geo_mean_long_sun = 280.46646 + jdc * (36000.76983 + 0.0003032 * jdc)

        while geo_mean_long_sun < 0:
            geo_mean_long_sun += 360
        while geo_mean_long_sun > 360:
            geo_mean_long_sun -= 360
        return geo_mean_long_sun
   def geo_mean_anomaly_sun(self, jdc):
       return 357.52911 + JDC * (35999.05029 - 0.0001537 * JDC)
   
  
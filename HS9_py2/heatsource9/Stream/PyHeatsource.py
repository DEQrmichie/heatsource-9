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

"""PyHeatsource is the core of the model and where the 
flux and hydraulic calculations are made."""

from __future__ import division, print_function
from math import pow, sqrt, log, log10, exp, pi
from math import atan, sin, cos, tan, acos, radians
from random import randint
from bisect import bisect

import logging
logger = logging.getLogger(__name__)

from ..Utils.Printer import Printer as print_console

def CalcSolarPosition(lat, lon, hour, min, sec, offset,
                      JDC, heatsource8, radial_count):
    toRadians = pi/180.0
    toDegrees = 180.0/pi
    
    MeanObliquity = (23.0 + (26.0 + ((21.448 - JDC *
                                      (46.815 + JDC *
                                       (0.00059 - JDC * 0.001813))) /
                                     60.0)) / 60.0)
    Obliquity = (MeanObliquity + 0.00256 *
                 cos(toRadians*(125.04 - 1934.136 * JDC)))
    
    Eccentricity = 0.016708634 - JDC * (0.000042037 + 0.0000001267 * JDC)
    GeoMeanLongSun = 280.46646 + JDC * (36000.76983 + 0.0003032 * JDC)

    while GeoMeanLongSun < 0:
        GeoMeanLongSun += 360
    while GeoMeanLongSun > 360:
        GeoMeanLongSun -= 360
    GeoMeanAnomalySun = 357.52911 + JDC * (35999.05029 - 0.0001537 * JDC)

    Dummy1 = toRadians*GeoMeanAnomalySun
    Dummy2 = sin(Dummy1)
    Dummy3 = sin(Dummy2 * 2)
    Dummy4 = sin(Dummy3 * 3)
    
    SunEqofCenter = (Dummy2 * (1.914602 - JDC *
                               (0.004817 + 0.000014 * JDC)) + Dummy3 *
                     (0.019993 - 0.000101 * JDC) + Dummy4 * 0.000289)
    
    SunApparentLong = ((GeoMeanLongSun + SunEqofCenter) -
                       0.00569 - 0.00478 *
                       sin(toRadians*((125.04 - 1934.136 * JDC))))

    Dummy1 = sin(toRadians*Obliquity) * sin(toRadians*SunApparentLong)
    Declination = toDegrees*(atan(Dummy1 / sqrt(-Dummy1 * Dummy1 + 1)))

    SunRadVector = ((1.000001018 * (1 - pow(Eccentricity,2))) /
                    (1 + Eccentricity * cos(toRadians*
                                            (GeoMeanAnomalySun +
                                             SunEqofCenter))))

    #======================================================
    # Equation of time (minutes)
    Dummy = pow((tan(Obliquity * pi / 360)),2)
    Dummy1 = sin(toRadians*(2 * GeoMeanLongSun))
    Dummy2 = sin(toRadians*(GeoMeanAnomalySun))
    Dummy3 = cos(toRadians*(2 * GeoMeanLongSun))
    Dummy4 = sin(toRadians*(4 * GeoMeanLongSun))
    Dummy5 = sin(toRadians*(2 * GeoMeanAnomalySun))
    Et = (toDegrees*(4 * (Dummy * Dummy1 - 2 *
                         Eccentricity * Dummy2 + 4 * Eccentricity *
                         Dummy * Dummy2 * Dummy3 - 0.5 * pow(Dummy,2) *
                         Dummy4 - 1.25 * pow(Eccentricity,2) * Dummy5)))

    SolarTime = ((hour*60.0) + min + (sec/60.0) +
                 (Et - 4.0 * -lon - (offset*60.0)))

    while SolarTime > 1440.0:
        SolarTime -= 1440.0
    HourAngle = SolarTime / 4.0 - 180.0
    if HourAngle < -180.0:
        HourAngle += 360.0

    Dummy = (sin(toRadians*lat) * sin(toRadians*Declination) +
             cos(toRadians*lat) * cos(toRadians*Declination) *
             cos(toRadians*HourAngle))
    
    if Dummy > 1.0:
        Dummy = 1.0
    elif Dummy < -1.0:
        Dummy = -1.0

    Zenith = toDegrees*(acos(Dummy))
    Dummy = cos(toRadians*lat) * sin(toRadians*Zenith)
    if abs(Dummy) >= 0.000999:
        Azimuth = ((sin(toRadians*lat) * cos(toRadians*Zenith) -
                    sin(toRadians*Declination)) / Dummy)
        
        if abs(Azimuth) > 1.0:
            if Azimuth < 0:
                Azimuth = -1.0
            else:
                Azimuth = 1.0

        Azimuth = 180 - toDegrees*(acos(Azimuth))
        if HourAngle > 0:
            Azimuth *= -1.0
    else:
        if lat > 0:
            Azimuth = 180.0
        else:
            Azimuth = 0.0
    if Azimuth < 0:
        Azimuth += 360.0

    AtmElevation = 90 - Zenith
    if AtmElevation > 85:
        RefractionCorrection = 0
    else:
        Dummy = tan(toRadians*(AtmElevation))
        if AtmElevation > 5:
            
            RefractionCorrection = (58.1 / Dummy - 0.07 /
                                    pow(Dummy,3) + 0.000086 /
                                    pow(Dummy,5))
        elif AtmElevation > -0.575:
            
            RefractionCorrection = (1735 + AtmElevation *
                                    (-518.2 + AtmElevation *
                                     (103.4 + AtmElevation *
                                      (-12.79 + AtmElevation * 0.711))))
        else:
            RefractionCorrection = -20.774 / Dummy
        RefractionCorrection = RefractionCorrection / 3600

    Zenith = Zenith - RefractionCorrection
    Altitude = 90 - Zenith
    Daytime = 0
    if Altitude > 0.0:
        Daytime = 1
    
    # Determine which landcover transect direction corresponds to the sun azimuth
    if heatsource8:  # same as 8 directions but no north
        tran = bisect((0.0,67.5,112.5,157.5,202.5,247.5,292.5),Azimuth)-1
        Azimuth_mod = Azimuth
    else:        
        Angle_Incr = 360.0 / radial_count
        DirNumbers = range(1, radial_count+1)
        AngleStart = [x*Angle_Incr-Angle_Incr/2 for x in DirNumbers]
        if Azimuth < AngleStart[0]:
            Azimuth_mod = Azimuth + 360
        else:
            Azimuth_mod = Azimuth
        tran = bisect(AngleStart,Azimuth_mod)-1
    return Altitude, Zenith, Daytime, tran, Azimuth_mod

def GetStreamGeometry(Q_est, W_b, z, n, S, D_est, dx, dt):
    Converge = 10
    dy = 0.01
    count = 0
    power = 2/3
    
    # ASSUMPTION: Make bottom width 1 cm to prevent undefined numbers 
    # in the math.
    W_b = 0.01 if W_b == 0 else W_b
    
    if D_est == 0:
        # This is a secant iterative solution method. It uses the
        # equation for discharge at equality then adds a slight change
        # to the depth and solves it again. It should iterate for a 
        # solution to depth within about 5-6 solutions.
        while Converge > 1e-7:
            Fy = ((D_est * (W_b + z * D_est)) *
                  pow(((D_est * (W_b + z * D_est)) /
                       (W_b + 2 * D_est * sqrt(1+ pow(z,2)))),power) -
                  ((n * Q_est) / sqrt(S)))
            
            thed = D_est + dy
            Fyy = ((thed * (W_b + z * thed)) *
                   pow((thed * (W_b + z * thed)) /
                       (W_b + 2 * thed * sqrt(1+ pow(z,2))),power) -
                   (n * Q_est) / sqrt(S))
            
            dFy = (Fyy - Fy) / dy
            if dFy <= 0: dFy = 0.99
            D_est -= Fy / dFy

            if (D_est < 0) or (D_est > 5000) or (count > 10000):
                # Damn, missed it. There may be a local minimum confusing
                # us, so we choose another depth at random and try again.
                D_est = randint(1,100)
                Converge = 0
                count = 0
            Converge = abs(Fy/dFy)
            count += 1
    # Use the calculated wetted depth to calculate new 
    # channel characteristics
    A = (D_est * (W_b + z * D_est))
    Pw = (W_b + 2 * D_est * sqrt(1+ pow(z,2)))
    Rh = A/Pw
    Ww = W_b + 2 * z * D_est
    U = Q_est / A

    # This is a sheer velocity estimate, followed by an estimate 
    # of numerical dispersion
    if S == 0.0:
        Shear_Velocity = U
    else:
        Shear_Velocity = sqrt(9.8 * D_est * S)
    Dispersion = ((0.011 * pow(U,2.0) * pow(Ww,2.0)) /
                  (D_est * Shear_Velocity))
    
    if (Dispersion * dt / pow(dx,2.0)) > 0.5:
        Dispersion = (0.45 * pow(dx,2)) / dt
    #Dispersion = 50
    return D_est, A, Pw, Rh, Ww, U, Dispersion

def CalcMuskingum(Q_est, U, W_w, S, dx, dt):
    """Return the values for the Muskigum routing coefficients
    using current timestep and optional discharge"""
    # Calculate an initial geometry based on an estimated 
    # discharge (typically (t,x-1))
    # Taken from the VB source in Heat Source version 7.
    c_k = (5/3) * U  # Wave celerity
    X = 0.5 * (1 - Q_est / (W_w * S * dx * c_k))
    if X > 0.5: X = 0.5
    elif X < 0.0: X = 0.0
    K = dx / c_k
    dt = dt

    # Check the celerity to ensure stability. These tests are 
    # from the VB code.
    if dt >= (2 * K * (1 - X)):
        # Unstable: Decrease dt or increase dx
        dt_stable = (2 * K * (1 - X)) / 60
        msg = "Unstable timestep. Decrease dt or increase dx. \
        dT must be < {0}, K={1}, X={2}".format(dt_stable, K, X)
        logger.error(msg)
        raise Exception(msg)

    # These calculations are from Chow's "Applied Hydrology"
    D = K * (1 - X) + 0.5 * dt
    C1 = (0.5*dt - K * X) / D
    C2 = (0.5*dt + K * X) / D
    C3 = (K * (1 - X) - 0.5*dt) / D
    # TODO: reformulate this using an updated model, 
    # such as Moramarco, et.al., 2006
    return C1, C2, C3

def CalcFlows(U, W_w, W_b, S, dx, dt, z, n, D_est, Q, Q_up, Q_up_prev,
              inputs, Q_bc):
    if Q_bc >= 0:
        Q_new = Q_bc
    else:
        Q1 = Q_up + inputs
        Q2 = Q_up_prev + inputs
        C = CalcMuskingum(Q2, U, W_w, S, dx, dt)
        Q_new = C[0]*Q1 + C[1]*Q2 + C[2]*Q

    #if Q_new > 0.000:
    Geom = GetStreamGeometry(Q_new, W_b, z, n, S, D_est, dx, dt)
    return Q_new, Geom

def GetSolarFlux(hour, JD, Altitude, Zenith, cloud, d_w, W_b, elevation,
                 TopoFactor, ViewToSky, transsample_distance, transsample_count,
                 BeersData, phi, emergent, lc_canopy, lc_height, lc_k,
                 ShaderList, tran, heatsource8):
    """ """
    FullSunAngle, TopoShadeAngle, BankShadeAngle, VegetationAngle = ShaderList
    F_Direct = [0]*8
    F_Diffuse = [0]*8
    F_Solar = [0]*8

    # Make all math functions local to save time by preventing failed 
    # searches of local, class and global namespaces
    #======================================================
    # 0 - Edge of atmosphere
    # TODO: Original VB code's JulianDay calculation:
    # JulianDay = -DateDiff("d", theTime, DateSerial(year(theTime), 1, 1))
    # This calculation for Rad_Vec should be checked, with 
    # respect to the DST hour/24 part.
    Rad_Vec = 1 + 0.017 * cos((2 * pi / 365) * (186 - JD + hour / 24))
    Solar_Constant = 1367 # W/m2
    
    # Global Direct Solar Radiation
    F_Direct[0] = ((Solar_Constant / (Rad_Vec ** 2)) *
                   sin(radians(Altitude)))
    
    F_Diffuse[0] = 0
    #======================================================
    # 1 - Above Topography
    Air_Mass = (35 / sqrt(1224 * sin(radians(Altitude)) + 1)) * \
        exp(-0.0001184 * elevation)
    Trans_Air = 0.0685 * cos((2 * pi / 365) * (JD + 10)) + 0.8
    #Calculate Diffuse Fraction
    F_Direct[1] = F_Direct[0] * (Trans_Air ** Air_Mass) * (1 - 0.65 *
                                                           cloud ** 2)
    if F_Direct[0] == 0:
        Clearness_Index = 1
    else:
        Clearness_Index = F_Direct[1] / F_Direct[0]

    Dummy = F_Direct[1]
    Diffuse_Fraction = (0.938 + 1.071 * Clearness_Index) - \
        (5.14 * (Clearness_Index ** 2)) + \
        (2.98 * (Clearness_Index ** 3)) - \
        (sin(2 * pi * (JD - 40) / 365)) * \
        (0.009 - 0.078 * Clearness_Index)
    F_Direct[1] = Dummy * (1 - Diffuse_Fraction)
    F_Diffuse[1] = Dummy * (Diffuse_Fraction) * (1 - 0.65 * cloud ** 2)

    #======================================================
    # 2 - Below Topography
    # 3 - Below Landcover (Above Bank Shade & Emergent)
    
    # need to add 1 because zero is emergent in stack data,
    # TODO fix. this should be consistent across the codebase
    tran = tran + 1
    Solar_blocked_byVeg = [0]*transsample_count
    
    if Altitude <= TopoShadeAngle:
        # Topographic shade is occurring
        F_Direct[2] = 0
        F_Diffuse[2] = F_Diffuse[1] * (1 - TopoFactor)
        F_Direct[3] = 0
    elif Altitude >= FullSunAngle:
        # Full sun
        F_Direct[2] = F_Direct[1]
        F_Diffuse[2] = F_Diffuse[1] * (1 - TopoFactor)
        F_Direct[3] = F_Direct[2]
    else:
        #======================================================
        # Topographic Shade is not occurring and
        # Partial shade from veg
        F_Direct[2] = F_Direct[1]
        F_Diffuse[2] = F_Diffuse[1] * (1 - TopoFactor)
        
        # 3 - Below Landcover (Above Bank Shade & Emergent)
        Dummy1 = F_Direct[2]
        
        zone = transsample_count - 1
        
        # zonal sun path length
        PLz = transsample_distance / cos(radians(Altitude))
        
        if heatsource8:
            #---------  Boyd and Kasper 2007
            PL = 10
        else:
            #--------- (Norman and Welles 1983)
            # here the extinction coefficient 
            # represents the foliage density * k
            # Calculate the path length as the total sun vector distance
            # through all the canopy zones
            # start with zero
            PL = 0
        
            for z in range(transsample_count): 
                if Altitude < VegetationAngle[z]:
                    # sun vector is passing through the canopy
                    PL += PLz
            
            # Just testing hs8 methods
            #PL = 10
        
        # Now calculate the fraction of radiation passed through the canopy
        while zone >= 0:
            if Altitude >= VegetationAngle[zone]:
                # no shading
                fraction_passed = 1
            else:
                # shading is occurring from this zone
                if BeersData == "LAI":
                    # use LAI data
                    # k as an input
                    RipExtinction = lc_k[tran][zone]
                    
                    # Assumption here is that LAI is for the whole canopy.
                    # Since we are calculating the fraction passed for a 
                    # portion of the canopy (zone) I take the approach 
                    # described by Norman and Welles 1983.
                    # Use total path length to calculate foliage 
                    # density, then calculate fraction passed for the 
                    # zonal path length.
                    fraction_passed = exp(-1 *
                                          RipExtinction *
                                          (lc_canopy[tran][zone] / PL) *
                                          PLz)
                else:
                    # Use canopy closure to calculate 
                    # the riparian extinction value
                    try:
                        RipExtinction = -log(1- lc_canopy[tran][zone])/ PL
                        fraction_passed = exp(-1* RipExtinction * PLz)
                    except:
                        if lc_canopy[tran][zone] >= 1:
                            # can't take log of zero
                            fraction_passed = 0
                        else:
                            # some other error
                            msg="Unknown error when calculating riparian extinction value. transect={0} zone={1} canopy={2} PL={3} ".format(tran,zone,lc_canopy[tran][zone],PL)
                            logger.error(msg)
                            print_console(msg)
                            raise
                        
            Solar_blocked_byVeg[zone] = Dummy1 - (Dummy1 * fraction_passed)
            Dummy1 *= fraction_passed
            zone -= 1
        F_Direct[3] = Dummy1
        
    F_Diffuse[3] = F_Diffuse[2] * ViewToSky
    diffuse_blocked = F_Diffuse[2]-F_Diffuse[3]
    Solar_blocked_byVeg.append(diffuse_blocked)
            
    #=========================================================
    # 4 - At Stream Surface (Below Bank Shade & Emergent)
    # What a Solar Pathfinder measures
    
    if Altitude > TopoShadeAngle and Altitude <= BankShadeAngle:
        # Bank shade is occurring
        F_Direct[4] = 0
        F_Diffuse[4] = F_Diffuse[3]
    else:
        # bank shade is not occurring
        F_Direct[4] = F_Direct[3]
        F_Diffuse[4] = F_Diffuse[3]

    if emergent:
        # Account for emergent vegetation
        if (lc_height[0][0] <= 0) or (lc_canopy[0][0] == 0):
            # Set to one if no veg or no canopy
            fraction_passed = 1
        
        else:
            if heatsource8:
                #---------  Boyd and Kasper 2007
                PLe = 10
            else:
                #--------- (Norman and Welles 1983)
                # PLe is the path length of the sun vector through 
                # the vegetation in the emergent zone
                PLe =  lc_height[0][0] / sin(radians(Altitude))
                
                # PLe can't be longer than the wetted width.
                # this happens when altitude is lower than height of veg.
                # if this is a solar only run the wetted width is not 
                # calculated so W_b = 0. The zone width is used instead.
                if W_b == 0: W_b = transsample_distance
                if PLe >  W_b: PLe = W_b

            if BeersData == "LAI":
                # use LAI data
                # k as an input
                RipExtinction = lc_k[0][0]
                fraction_passed = exp(-1* lc_canopy[0][0] * RipExtinction * PLe)
                
            else:
                try:
                    # Use canopy closure to calculate 
                    # the riparian extinction value
                    RipExtinction = -log(1- lc_canopy[0][0]) / PLe
                    fraction_passed = exp(-1* RipExtinction * PLe)
                except:
                    if lc_canopy[0][0] >= 1:
                        # can't take log of zero
                        fraction_passed = 0
                    else:
                        # some other error
                        msg="Unknown error when calculating emergent riparian extinction value. canopy={0} PLe={1} ".format(lc_canopy[0][0],PLe)
                        logger.error(msg)
                        print_console(msg)
                        raise Exception
                
        F_Direct[4] = F_Direct[4] * fraction_passed
        F_Diffuse[4] = F_Diffuse[4] * fraction_passed

    #=========================================================
    # 5 - Entering Stream
    if Zenith > 80:
        Stream_Reflect = 0.0515 * (Zenith) - 3.636
    else:
        Stream_Reflect = 0.091 * (1 / cos(Zenith * pi / 180)) - 0.0386
    if abs(Stream_Reflect) > 1:
        Stream_Reflect = 0.0515 * (Zenith * pi / 180) - 3.636
    if abs(Stream_Reflect) > 1:
        Stream_Reflect = 0.091 * (1 / cos(Zenith * pi / 180)) - 0.0386
    F_Diffuse[5] = F_Diffuse[4] * 0.91
    F_Direct[5] = F_Direct[4] * (1 - Stream_Reflect)
    
    #=========================================================
    # 6 - Received by Water Column
    # 7 - Received by Bed
    # Jerlov (1976)
    Water_Path = (d_w / cos(atan((sin(radians(Zenith)) / 1.3333) /
                                 sqrt(-(sin(radians(Zenith)) / 1.3333) *
                                      (sin(radians(Zenith)) / 1.3333) + 1))))
    
    Trans_Stream = 0.415 - (0.194 * log10(Water_Path * 100))
    if Trans_Stream > 1:
        Trans_Stream = 1
    
    # Direct Solar Radiation attenuated on way down
    Dummy1 = F_Direct[5] * (1 - Trans_Stream)
    
    # Direct Solar Radiation Hitting Stream bed
    Dummy2 = F_Direct[5] - Dummy1
    
    # Reflection Coef. for Direct Solar
    Bed_Reflect = exp(0.0214 * (Zenith * pi / 180) - 1.941)
    BedRock = 1 - phi
    
    # Direct Solar Radiation Absorbed in Bed
    Dummy3 = Dummy2 * (1 - Bed_Reflect)                
    
    # Direct Solar Radiation Immediately Returned to Water Column as Heat
    Dummy4 = 0.53 * BedRock * Dummy3                   
    
    # Direct Solar Radiation Reflected off Bed
    Dummy5 = Dummy2 * Bed_Reflect                      
    
    # Direct Solar Radiation attenuated on way up
    Dummy6 = Dummy5 * (1 - Trans_Stream)               
    
    F_Direct[6] = Dummy1 + Dummy4 + Dummy6
    F_Direct[7] = Dummy3 - Dummy4
    Trans_Stream = 0.415 - (0.194 * log10(100 * d_w))
    if Trans_Stream > 1:
        Trans_Stream = 1
    
    # Diffuse Solar Radiation attenuated on way down
    Dummy1 = F_Diffuse[5] * (1 - Trans_Stream)
    
    # Diffuse Solar Radiation Hitting Stream bed
    Dummy2 = F_Diffuse[5] - Dummy1                  
    
    # Reflection Coef. for Diffuse Solar
    # TODO: The following ALWAYS becomes exp(-1.941)
    Bed_Reflect = exp(0.0214 * (0) - 1.941)
    
    # Diffuse Solar Radiation Absorbed in Bed
    Dummy3 = Dummy2 * (1 - Bed_Reflect)
    
    # Diffuse Solar Radiation Immediately Returned to Water Column as Heat
    Dummy4 = 0.53 * BedRock * Dummy3
    
    # Diffuse Solar Radiation Reflected off Bed
    Dummy5 = Dummy2 * Bed_Reflect
    
    # Diffuse Solar Radiation attenuated on way up
    Dummy6 = Dummy5 * (1 - Trans_Stream)
    F_Diffuse[6] = Dummy1 + Dummy4 + Dummy6
    F_Diffuse[7] = Dummy3 - Dummy4
    
    #=========================================================
    # Flux_Solar(x) and Flux_Diffuse = Solar flux at various positions
    # 0 - Edge of atmosphere
    # 1 - Above Topography
    # 2 - Below Topography
    # 3 - Below Landcover (Above Bank Shade & Emergent)
    # 4 - At Stream Surface (What a Solar Pathfinder Measures)
    # 5 - Entering Stream
    # 6 - Received by Water Column
    # 7 - Received by Bed
    
    F_Solar[0] = F_Diffuse[0] + F_Direct[0]
    F_Solar[1] = F_Diffuse[1] + F_Direct[1]
    F_Solar[2] = F_Diffuse[2] + F_Direct[2]
    F_Solar[3] = F_Diffuse[3] + F_Direct[3]
    F_Solar[4] = F_Diffuse[4] + F_Direct[4]
    F_Solar[5] = F_Diffuse[5] + F_Direct[5]
    F_Solar[6] = F_Diffuse[6] + F_Direct[6]
    F_Solar[7] = F_Diffuse[7] + F_Direct[7]
    
    return F_Solar, F_Diffuse, F_Direct, Solar_blocked_byVeg

def GetGroundFluxes(cloud, wind, humidity, T_air, elevation, phi,
                    lc_height, ViewToSky, SedDepth, dx, dt, SedThermCond,
                    SedThermDiff, calcalluv, T_alluv, P_w, W_w, emergent,
                    penman, wind_a, wind_b, calcevap, T_prev, T_sed,
                    Q_hyp, F_Solar5, F_Solar7):

    # SedThermCond units of W/(m *C)
    # SedThermDiff units of cm^2/sec

    SedRhoCp = SedThermCond / (SedThermDiff / 10000)
    # NOTE: SedRhoCp is the product of sediment density and heat capacity
    # since thermal conductivity is defined as 
    # density * heat capacity * diffusivity,
    # therefore (density * heat capacity) = (conductivity / diffusivity)
    # units of (J / m3 / *C)

    # Water Variable
    rhow = 1000  #density of water kg / m3
    H2O_HeatCapacity = 4187 #J/(kg *C)

    # Conduction flux (positive is heat into stream)
    # units of (W/m2)
    F_Cond = SedThermCond * (T_sed - T_prev) / (SedDepth / 2) 
    
    # Calculate the conduction flux between deeper alluvium 
    # & substrate conditionally
    Flux_Conduction_Alluvium = SedThermCond * (T_sed - T_alluv) / (SedDepth / 2) if calcalluv else 0.0

    # Hyporheic flux (negative is heat into sediment)
    F_hyp = Q_hyp * rhow * H2O_HeatCapacity * (T_sed - T_prev) / (W_w * dx)

    NetFlux_Sed = F_Solar7 - F_Cond - Flux_Conduction_Alluvium - F_hyp
    DT_Sed = NetFlux_Sed * dt / (SedDepth * SedRhoCp)
    T_sed_new = T_sed + DT_Sed
    if T_sed_new > 50 or T_sed_new < 0:
        msg = "Sediment temperature is {0}. must be bounded in 0<=temp<=50".format(T_sed_new)
        logger.error(msg)
        raise Exception() # TODO RM

    #=====================================================
    # Calculate Longwave FLUX
    #=====================================================
    # Atmospheric variables
    
    # mbar (Chapra p. 567)
    Sat_Vapor = 6.1275 * exp(17.27 * T_air / (237.3 + T_air)) 
    Air_Vapor = humidity * Sat_Vapor
    
    # Stefan-Boltzmann constant (W/m2 K4)    
    Sigma = 5.67e-8
    
    # Dingman p 282
    Emissivity = (1.72 * (((Air_Vapor * 0.1) /
                           (273.2 + T_air)) ** (1 / 7)) *
                  (1 + 0.22 * cloud ** 2))
    #======================================================
    # Calculate the atmospheric longwave flux
    F_LW_Atm = 0.96 * ViewToSky * Emissivity * Sigma * (T_air + 273.2) ** 4
    # Calculate the backradiation longwave flux
    F_LW_Stream = -0.96 * Sigma * (T_prev + 273.2) ** 4
    # Calculate the vegetation longwave flux
    F_LW_Veg = 0.96 * (1 - ViewToSky) * 0.96 * Sigma * (T_air + 273.2) ** 4
    # Calculate the net longwave flux
    F_Longwave = F_LW_Atm + F_LW_Stream + F_LW_Veg

    #===================================================
    # Calculate Evaporation FLUX
    #===================================================
    # Atmospheric Variables
    Pressure = 1013 - 0.1055 * elevation #mbar
    
    # mbar (Chapra p. 567)
    Sat_Vapor = 6.1275 * exp(17.27 * T_prev / (237.3 + T_prev))
    Air_Vapor = humidity * Sat_Vapor
    #===================================================
    # Calculate the frictional reduction in wind velocity
    if emergent and lc_height[0][0] > 0:
        
        Zm = 2
        
        # zm > Zd + Zo
        if lc_height[0][0] < 2.5:
            Zd = 0.7 * lc_height[0][0]
            Zo = 0.1 * lc_height[0][0]

        else:
            Zd = 1.75
            Zo = 0.25

        # Vertical wind Decay Rate using Prandtl-von-Karman 
        # universal velocity distribution (Dingman 2002 p. 594)
        Friction_Velocity = (1 / 0.4) * wind * log((Zm - Zd) / Zo)
    else:
        Zo = 0.00023 #Brustsaert (1982) p. 277 Dingman
        Zd = 0 #Brustsaert (1982) p. 277 Dingman
        Zm = 2
        Friction_Velocity = wind
    #===================================================
    # Wind Function f(w)
    #m/mbar/s
    Wind_Function = float(wind_a) + float(wind_b) * Friction_Velocity 
    # Wind_Function = 0.000000001505 + 0.0000000016 * Friction_Velocity #m/mbar/s

    #===================================================
    # Latent Heat of Vaporization
    LHV = 1000 * (2501.4 + (1.83 * T_prev)) #J/kg
    #===================================================
    # Use Jobson Wind Function
    if penman:
        #Calculate Evaporation FLUX
        P = 998.2 # kg/m3
        Gamma = 1003.5 * Pressure / (LHV * 0.62198) #mb/*C  Cuenca p 141
        Delta = (6.1275 * exp(17.27 * T_air / (237.3 + T_air)) -
                 6.1275 * exp(17.27 * (T_air - 1) / (237.3 + T_air - 1)))
        
        NetRadiation = F_Solar5 + F_Longwave  #J/m2/s
        
        if NetRadiation < 0:
            NetRadiation = 0 #J/m2/s
        Ea = Wind_Function * (Sat_Vapor - Air_Vapor)  #m/s
        Evap_Rate = (((NetRadiation * Delta / (P * LHV)) + Ea * Gamma) /
                     (Delta + Gamma))
        
        F_Evap = -Evap_Rate * LHV * P #W/m2
        # Calculate Convection FLUX
        if (Sat_Vapor - Air_Vapor) != 0:
            Bowen = Gamma * (T_prev - T_air) / (Sat_Vapor - Air_Vapor)
        else:
            Bowen = 1        
    else:
        #===================================================
        # Calculate Evaporation FLUX
        Evap_Rate = Wind_Function * (Sat_Vapor - Air_Vapor)  #m/s
        P = 998.2 # kg/m3
        F_Evap = -Evap_Rate * LHV * P #W/m2
        # Calculate Convection FLUX
        if (Sat_Vapor - Air_Vapor) != 0:
            Bowen = (0.61 * (Pressure / 1000) * (T_prev - T_air) /
                     (Sat_Vapor - Air_Vapor))
        else:
            Bowen = 1
            
    F_Conv = F_Evap * Bowen
    E = Evap_Rate*W_w if calcevap else 0
    return F_Cond, T_sed_new, F_Longwave, F_LW_Atm, F_LW_Stream, F_LW_Veg, F_Evap, F_Conv, E

def CalcMacCormick(dt, dx, U, T_sed, T_prev, Q_hyp, Q_tup, T_tup, Q_up,
                   Delta_T, Disp, S1, S1_value, T0, T1, T2, Q_accr,
                   T_accr, MixTDelta_dn):
    Q_in = 0
    T_in = 0
    T_up = T0
    numerator = 0
    for i in xrange(len(Q_tup)):
        Qitem = Q_tup[i]
        Titem = T_tup[i]
        # make sure there's a value for discharge. Temp can be 
        # blank if discharge is negative (withdrawal)
        if Qitem is None or (Qitem > 0 and Titem is None):
            raise HeatSourceError("Problem with null value in tributary discharge or temperature")
        if Qitem > 0:
            Q_in += Qitem
            numerator += Qitem*Titem
    if numerator and (Q_in > 0):
        T_in = numerator/Q_in
    # This is basically MixItUp from the VB code
    T_mix = ((Q_in * T_in) + (T_up * Q_up)) / (Q_up + Q_in)
    
    # Calculate temperature change from mass transfer from hyporheic zone
    T_mix = (((T_sed * Q_hyp) + (T_mix * (Q_up + Q_in))) /
             (Q_hyp + Q_up + Q_in))
    
    # Calculate temperature change from accretion inflows
    # Q_hyp is commented out because we are not currently sure if 
    # it should be added to the flow. This is because adding it will 
    # cause overestimation of the discharge if Q_hyp is not subtracted
    # from the total discharge (Q_in) somewhere else, which it is not. 
    # We should check this eventually.
    T_mix = (((Q_accr * T_accr) + (T_mix * (Q_up + Q_in + Q_hyp))) /
             (Q_accr + Q_up + Q_in + Q_hyp))
    
    T_mix -= T_up
    
    # We need to adjust the upstream temperature by the tributary mixing
    # so the longitudinal slope of change in T is not over predicted.
    T0 += T_mix

    # Similarly we need to adjust the downstream temperature (T2) 
    # to account for mixing in that reach.
    T2 -= MixTDelta_dn

    Dummy1 = -U * (T1 - T0) / dx
    Dummy2 = Disp * (T2 - 2 * T1 + T0) / (dx**2)
    S = Dummy1 + Dummy2 + Delta_T / dt
    if S1:
        Temp = T_prev + ((S1_value + S) / 2) * dt
    else:
        Temp = T1 + S * dt

    return Temp, S, T_mix

def CalcHeatFluxes(metData, C_args, d_w, area, P_w, W_w, U, Q_tribs,
                   T_tribs, T_prev, T_sed, Q_hyp, T_dn_prev, ShaderList,
                   tran, Disp, hour, JD, daytime, Altitude, Zenith,
                   Q_up_prev, T_up_prev, solar_only, MixTDelta_dn_prev,
                   heatsource8):
    
    cloud, wind, humidity, T_air = metData

    W_b, elevation, TopoFactor, ViewToSky, phi, lc_canopy, lc_height, \
        lc_k, SedDepth, dx, dt, SedThermCond, SedThermDiff, Q_accr, \
        T_accr, has_prev, transsample_distance, transsample_count, \
        BeersData, emergent, wind_a, wind_b, calcevap, penman, \
        calcalluv, T_alluv = C_args

    solar = [0]*8
    diffuse = [0]*8
    direct = [0]*8
    
    # plus one for diffuse blocked
    veg_block = [0]*len(ShaderList[3])+[0]
    
    if daytime:
        
        (solar, diffuse,
        direct, veg_block) = GetSolarFlux(hour, JD, Altitude, Zenith,
                                        cloud, d_w, W_b, elevation,
                                        TopoFactor, ViewToSky,
                                        transsample_distance,
                                        transsample_count,
                                        BeersData, phi, emergent,
                                        lc_canopy, lc_height,
                                        lc_k, ShaderList, tran,
                                        heatsource8)

    # We're only running shade, so return solar and some empty calories
    if solar_only:
        ground = [0]*9
        F_Total =  0.0
        Delta_T = 0.0
        Mac = [0]*3
        
        # Boundary node
        if not has_prev: return solar, diffuse, direct, veg_block, ground, F_Total, Delta_T
        
        # regular node
        else: return solar, diffuse, direct, veg_block, ground, F_Total, Delta_T, Mac

    ground = GetGroundFluxes(cloud, wind, humidity, T_air, elevation,
                    phi, lc_height, ViewToSky, SedDepth, dx,
                    dt, SedThermCond, SedThermDiff, calcalluv, T_alluv,
                    P_w, W_w, emergent, penman, wind_a, wind_b,
                    calcevap, T_prev, T_sed, Q_hyp, solar[5],
                    solar[7])

    F_Total =  solar[6] + ground[0] + ground[2] + ground[6] + ground[7]
    
    # Vars are Cp (J/kg *C) and P (kgS/m3)
    Delta_T = F_Total * dt / ((area / W_w) * 4182 * 998.2) 

    if not has_prev:
        # Boundary node
        return solar, diffuse, direct, veg_block, ground, F_Total, Delta_T

    Mac = CalcMacCormick(dt, dx, U, ground[1], T_prev, Q_hyp, Q_tribs,
                         T_tribs, Q_up_prev, Delta_T, Disp, 0, 0.0,
                         T_up_prev, T_prev, T_dn_prev, Q_accr, T_accr,
                         MixTDelta_dn_prev)

    # Mac includes Temp, S, T_mix
    return solar, diffuse, direct, veg_block, ground, F_Total, Delta_T, Mac
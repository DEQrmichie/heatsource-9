# Heat Source, Copyright (C) 2000-2019, 
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

from builtins import range
from builtins import object
from math import sqrt
import logging
logger = logging.getLogger(__name__)

from . import pyheatsource as py_HS
from heatsource9.model.interpolator import Interpolator
from heatsource9.domain.clock import pretty_time

class StreamNode(object):
    """Definition of an individual stream segment"""
    def __init__(self, **kwargs):
        __slots = ["latitude", "longitude", "Zs", # Geographic params
                "FLIR_Temp", "FLIR_Time", # FLIR data
                "T_sed", "T_hyp", "T_accr", "T_tribs", # Temperature attrs
                "lc_height_top", "lc_canopy_cover", "lc_oh", "lc_canopy_depth", "lc_lai", "lc_k", # Land cover params
                "lc_height_node_top", # landcover height relative to the stream node elevation
                "metData", # Meteorological data
                "zm", # Meteorological measurement height
                "Zone", "T_bc", # Initialization parameters, Zone and boundary conditions
                "Delta_T", # Current temperature calculated from only local fluxes
                "Mix_T_Delta", #Change in temperature due to tribs, gw, points sources, accretion
                "T", "T_prev", # Current and previous stream temperature
                "TopoFactor", # Topo_W + Topo_S + Topo_E / (90*3). From Above stream surface solar flux calculations
                "ViewToSky", # Total angle of full sun view
                "ShaderList", # List of angles and attributes to determine sun shading.
                "F_DailySum", "F_Total", "Solar_Blocked", # Specific sums of solar fluxes
                "F_Solar", # List of important solar fluxes
                "F_Direct", 
                "F_Diffuse",
                "Ksed", "Alpha_sed", "Dsed", # Sediment conduction values
                "Q_hyp_frac", # Fraction hyporheic exchange                
                "S",        # Slope
                "n",        # Manning's n
                "Z", # Channel side slope ratio, ratio of run to rise of the side of a trapezoidal channel
                "Dw", # Wetted depth. Calculated in GetWettedDepth()
                "d_w_prev", # Wetted depth for the previous timestep
                "d_cont", # Control depth
                "Wb", # Bottom width
                "Ww", # Wetted width, calculated as Wb + 2*Z*Dw
                "A", # Cross-sectional Area, calculated Dw * (Wb + Z * Dw)
                "Pw", # Wetted Perimeter, calculated as Wb + 2 * Dw * sqrt(1 + Z**2)
                "Rh", # Hydraulic Radius, calculated as A/Pw
                "dx", # Length of this stream reach.
                "U", # Velocity from Manning's relationship
                "Q", # Discharge, from Manning's relationship
                "Q_prev", # Discharge at previous timestep, previous space step is taken from another node
                "Q_cont", # Control discharge
                "Vw", # Total volume, based on current flow
                "Q_tribs", # Inputs from tribs.
                "Q_accr", # Inputs from "accretion" in cubic meters per second
                "Q_with", # Withdrawals from the stream, in cubic meters per second
                "Q_hyp", # Hyporheic flow
                "km", # River kilometer, from mouth
                "nodeID",  # Node ID
                "streamID",  # Stream ID
                "Q_bc", # Boundary conditions, in a TimeList class, for discharge.
                "Q_evap", # Evaporation flow loss in cubic meters per second
                "dt", # This is the timestep (for kinematic wave movement, etc.)
                "Eta", # Porosity of the bed
                "Log",  # Global logging class
                "Disp", "S1", # Dispersion due to shear stress and placeholder calculation variable
                "next_km", "prev_km", "head", # reference placeholders for next, previous nodes and headwater node
                "Q_mass", # Local mass balance variable (StreamChannel level)
                "F_Conduction","F_Convection","F_LW","F_Evaporation", # Ground flux variables
                "F_LW_stream", "F_LW_atm", "F_LW_veg", # Longwave fluxes
                "C_args", # tuple of variables that do not change during the model
                "CalcHeat", "CalcDischarge", # Reference to correct heat calculation method
                "SolarPos", "UTC_offset",
                "run_params" # run params passed from ModelSetup
                ]
        # Define members in __slots__ to ensure that later member names cannot be added accidentally
        # Set all the attributes to bare lists, or set from the constructor
        for attr in __slots:
            x = kwargs[attr] if attr in list(kwargs.keys()) else None
            setattr(self, attr, x)
        self.__slots = __slots
        self.__slots.sort()
        self.run_params = kwargs.get("run_params", {})

        self.T = 0.0
        self.Mix_T_Delta = 0.0
        self.Q_mass = 0
        dt = self.run_params.get("dt", 1)
        self.metData = Interpolator(dt=dt)
        self.zm = 2.0
        self.T_tribs = Interpolator(dt=dt)
        self.Q_tribs = Interpolator(dt=dt)
        # Create an internal dictionary that we can pass to the C module,
        # this contains self.slots attributes and other 
        # things the C module needs
        for attr in ["F_Conduction","F_Convection","F_LW","F_Evaporation"]:
            setattr(self, attr, 0)
        self.F_Solar = [0.0]*8
        self.F_Direct = [0.0]*8
        self.F_Diffuse = [0.0]*8
        self.F_Total = 0.0
        self.ShaderList = ()
        if self.run_params.get("heatsource8", False):
            radial_count = 7
        else:
            radial_count = self.run_params.get("trans_count", 0)
        transsample_count = self.run_params.get("transsample_count", 0)
        self.lc_height_top = [[[0.0]for zone in range(transsample_count)] for tran in range(radial_count + 1)]
        self.lc_height_node_top = [[[0.0]for zone in range(transsample_count)] for tran in range(radial_count + 1)]
        self.lc_canopy_cover = [[[0.0]for zone in range(transsample_count)] for tran in range(radial_count + 1)]
        self.lc_lai = [[[0.0]for zone in range(transsample_count)] for tran in range(radial_count + 1)]
        self.lc_oh = [[[0.0]for zone in range(transsample_count)] for tran in range(radial_count + 1)]
        self.lc_canopy_depth = [[[0.0]for zone in range(transsample_count)] for tran in range(radial_count + 1)]
        self.lc_k = [[[0.0]for zone in range(transsample_count)] for tran in range(radial_count + 1)]
        self.UTC_offset = self.run_params.get("offset", 0)
    def get_node_data(self):
        data = {}
        for attr in self.__slots:
            data[attr] = getattr(self,attr)
        return data

    def __eq__(self, other):
        cmp = other.km if isinstance(other, StreamNode) else other
        return self.km == cmp
    def __ne__(self, other):
        cmp = other.km if isinstance(other, StreamNode) else other
        return self.km != cmp
    def __gt__(self, other):
        cmp = other.km if isinstance(other, StreamNode) else other
        return self.km > cmp
    def __lt__(self, other):
        cmp = other.km if isinstance(other, StreamNode) else other
        return self.km < cmp
    def __ge__(self, other):
        cmp = other.km if isinstance(other, StreamNode) else other
        return self.km >= cmp
    def __le__(self, other):
        cmp = other.km if isinstance(other, StreamNode) else other
        return self.km <= cmp
    def __repr__(self):
        return '%s @ %.3f km' % (self.__class__.__name__, self.km)

    def initialize(self):
        """Methods necessary to set initial conditions of the node"""
        global py_HS
        has_prev = self.prev_km is not None
        if has_prev:
            self.CalcHeat = self.calc_heat_opt
        else:
            self.CalcHeat = self.calc_heat_boundary_node

        self.CalcDischarge = self.calculate_discharge
        rp = self.run_params
        self.C_args = (self.Wb, self.Zs, self.TopoFactor,
                       self.ViewToSky, self.Eta, self.lc_canopy_cover,
                       self.lc_lai, self.lc_height_top, self.lc_height_node_top, self.lc_k, self.lc_oh, self.lc_canopy_depth, self.Dsed,
                       self.dx, self.dt, self.Ksed,
                       self.Alpha_sed, self.Q_accr, self.T_accr,
                       has_prev, rp.get("transsample_distance", 0.0),
                       rp.get("transsample_count", 0),
                       rp.get("canopy_data", "LAI"), rp.get("lcsampmethod", "point"), rp.get("emergent", False),
                       rp.get("wind_a", 0.0), rp.get("wind_b", 0.0),
                       rp.get("calcevap", False), rp.get("penman", False),
                       rp.get("calcalluvium", False), self.zm,
                       rp.get("alluviumtemp", 0.0))

    def _raise_discharge_error(self, exc, time, Q_net, up_q=None, q_bc=None):
        msg = (
            "Discharge calculation failed at "
            f"km={self.km}, model time={pretty_time(time)}, "
            f"Q={self.Q}, Q_net={Q_net}, "
            f"Wb={self.Wb}, S={self.S}, dx={self.dx}, dt={self.dt}"
        )
        if up_q is not None:
            msg += f", up_Q={up_q}"
        if q_bc is not None:
            msg += f", Q_bc={q_bc}"
        raise RuntimeError(msg) from exc

    def calc_discharge_opt(self, time):
        """A Version of calculate_discharge() that does not require
        checking for boundary conditions"""
        Q_net = self.Q_accr + sum(self.Q_tribs[time]) - self.Q_with - self.Q_evap
        self.Q_mass += Q_net
        up = self.prev_km

        try:
            Q, (self.Dw, self.A, self.Pw, self.Rh, self.Ww, self.U,
                self.Disp) = \
                py_HS.calc_flows(self.U, self.Ww, self.Wb, self.S,
                              self.dx, self.dt, self.Z, self.n,
                              self.d_cont, self.Q, up.Q, up.Q_prev,
                              Q_net, -1)
        except RuntimeError as exc:
            self._raise_discharge_error(exc, time, Q_net, up_q=up.Q)
        
        self.Q_prev = self.Q
        self.Q = Q
        
        # Hyporheic discharge
        self.Q_hyp = Q * self.Q_hyp_frac 

    def calc_discharge_boundary_node(self, time):
        Q_bc = self.Q_bc[time]
        self.Q_mass += Q_bc
        # We fill the discharge arguments with 0 because it is 
        # unused in the boundary case

        try:
            Q, (self.Dw, self.A, self.Pw, self.Rh, self.Ww,
                self.U, self.Disp) = \
                    py_HS.calc_flows(self.U, self.Ww, self.Wb, self.S,
                                  self.dx, self.dt, self.Z, self.n,
                                  self.d_cont,
                                  0.0, 0.0, 0.0, 0.0, Q_bc)
        except RuntimeError as exc:
            self._raise_discharge_error(exc, time, 0.0, q_bc=Q_bc)
        
        # Now we've got a value for Q(t,x), so the current Q becomes Q_prev.
        self.Q_prev = self.Q
        self.Q = Q
        self.Q_hyp = Q * self.Q_hyp_frac # Hyporheic discharge
    def calculate_discharge(self, time):
        """Return the discharge for the current timestep

        This method uses the GetKnownDischarges() and GetMuskigum()
        methods (in C module) to grab the values necessary to calculate
        the discharge at the current timestep for the channel. If we
        are at a boundary (spatial or temporal) the appropriate boundary
        condition is returned. When the new discharge is calculated, the
        previous discharge value is placed in Q_prev for use by the
        downstream channel. This method makes some assumptions, one is
        that Q_bc is a TimeList instance holding boundary conditions
        for the given node, and that this is only True if this node has
        no upstream channel. Two is that the values for Q_tribs is a
        TimeList instance or None, that Q_accr and Q_with are values in
        cubic meters per second of inputs and withdrawals or None. The
        argument t is for a Python datetime object and can (should) be
        None if we are not at a spatial boundary. dt is the timestep in
        minutes, which cannot be None.
        """
        Q_net = self.Q_accr + sum(self.Q_tribs[time]) - self.Q_with - self.Q_evap
        # Check if we are a spatial or temporal boundary node
        if self.prev_km:
            # There's an upstream channel, but no previous timestep.
            # In this case, we sum the incoming flow which is upstream's
            # current timestep plus net inflow.
            
            # Add upstream node's discharge at THIS 
            # timestep- prev_km.Q would be next timestep.
            Q = self.prev_km.Q_prev + Q_net

            try:
                Q, (self.Dw, self.A, self.Pw, self.Rh, self.Ww, self.U, self.Disp) = \
                    py_HS.calc_flows(0.0, 0.0, self.Wb, self.S, self.dx,
                                   self.dt, self.Z, self.n, self.d_cont, 0.0, 0.0, 0.0, Q_net, Q)
            except RuntimeError as exc:
                self._raise_discharge_error(exc, time, Q_net, up_q=self.prev_km.Q_prev)

            # If we hit this once, we remap so we can avoid the 
            # if statements in the future.
            self.CalcDischarge = self.calc_discharge_opt
        else:
            # We're a spatial boundary, use the boundary condition
            # At spatial boundaries, we return the boundary 
            # conditions from Q_bc
            Q_bc = self.Q_bc[time]
            self.Q_mass += Q_bc
            # We pad the arguments with 0 because some are 
            # unused (or currently None) in the boundary case

            try:
                Q, (self.Dw, self.A, self.Pw, self.Rh, self.Ww, self.U, self.Disp) = \
                    py_HS.calc_flows(0.0, 0.0, self.Wb, self.S, self.dx, self.dt, self.Z, self.n,
                                   self.d_cont, 0.0, 0.0, 0.0, Q_net, Q_bc)
            except RuntimeError as exc:
                self._raise_discharge_error(exc, time, Q_net, q_bc=Q_bc)

            self.CalcDischarge = self.calc_discharge_boundary_node

        # Now we've got a value for Q(t,x), so the current Q becomes Q_prev.
        self.Q_prev = self.Q  or Q
        self.Q = Q
        self.Q_hyp = Q * self.Q_hyp_frac # Hyporheic discharge

    def calc_heat_opt(self, time, hour, min, sec, doy,
                     JC, solar_only=False):
        """Inlined version of CalcHeat optimized for non-boundary nodes
        (removes a bunch of if/else statements)"""
        # The solar_only keyword is a kludge to turn off calculation of the
        # ground fluxes in the case of shade-only runs.

        # Reset temperatures
        self.T_prev = self.T
        self.T = None
        
        # For each node we set the solar position variables to be the 
        # same as the headwater node in order to save on processing time.
        (SolarAltitude, 
         SolarZenith, 
         Daytime, 
         tran, 
         Azimuth_mod) = self.head.SolarPos
        
        (self.F_Solar,
         self.F_Diffuse,
         self.F_Direct,
         veg_block,
         ground,
         self.F_Total,
         self.Delta_T,
         Mac) = py_HS.calc_heat_fluxes(self.metData[time],
                                   self.C_args,
                                   self.Dw,
                                   self.A,
                                   self.Pw,
                                   self.Ww,
                                   self.U,
                                   self.Q_tribs[time],
                                   self.T_tribs[time],
                                   self.T_prev,
                                   self.T_sed,
                                   self.Q_hyp,
                                   self.next_km.T_prev,
                                   self.ShaderList[tran],
                                   tran,
                                   self.Disp, 
                                   hour,
                                   doy,
                                   Daytime,
                                   SolarAltitude,
                                   SolarZenith,
                                   self.prev_km.Q_prev,
                                   self.prev_km.T_prev,
                                   solar_only,
                                   self.next_km.Mix_T_Delta,
                                   self.run_params.get("heatsource8", False))
        
        (self.F_Conduction,
         self.T_sed,
         self.F_LW,
         self.F_LW_atm,
         self.F_LW_stream,
         self.F_LW_veg,
         self.F_Evaporation,
         self.F_Convection,
         self.Q_evap) = ground
        self.T_hyp = self.T_sed
        
        (self.T, self.S1, self.Mix_T_Delta) = Mac
        
        self.F_DailySum[1] += self.F_Solar[1]
        self.F_DailySum[4] += self.F_Solar[4]
        
        for i in range(len(self.Solar_Blocked[tran])):
            self.Solar_Blocked[tran][i] += veg_block[i]
        self.Solar_Blocked["diffuse"] += veg_block[-1]

    def calc_heat_boundary_node(self, time, hour, min, sec, doy,
                              JC, solar_only=False):
        # Reset temperatures
        self.T_prev = self.T
        self.T = None
        
        (SolarAltitude,
         SolarZenith,
         Daytime,
         tran,
         Azimuth_mod) = py_HS.calc_solar_position(self.latitude,
                                              self.longitude,hour,
                                              min,
                                              sec,
                                              self.UTC_offset,
                                              JC,
                                              self.run_params.get("heatsource8", False),
                                              self.run_params.get("trans_count", 0))
        
        self.SolarPos = SolarAltitude, SolarZenith, Daytime, tran, Azimuth_mod
        
        (self.F_Solar,
         self.F_Diffuse,
         self.F_Direct,
         veg_block,
         ground,
         self.F_Total,
         self.Delta_T) = py_HS.calc_heat_fluxes(self.metData[time],
                                       self.C_args,
                                       self.Dw,
                                       self.A,
                                       self.Pw,
                                       self.Ww,
                                       self.U,
                                       self.Q_tribs[time],
                                       self.T_tribs[time],
                                       self.T_prev,
                                       self.T_sed,
                                       self.Q_hyp,
                                       self.next_km.T_prev,
                                       self.ShaderList[tran],
                                       tran,
                                       self.Disp,
                                       hour,
                                       doy,
                                       Daytime,
                                       SolarAltitude,
                                       SolarZenith,
                                       0.0,
                                       0.0,
                                       solar_only,
                                       self.next_km.Mix_T_Delta,
                                       self.run_params.get("heatsource8", False))
        
        (self.F_Conduction, 
         self.T_sed,
         self.F_LW,
         self.F_LW_atm,
         self.F_LW_stream,
         self.F_LW_veg,
         self.F_Evaporation,
         self.F_Convection,
         self.Q_evap) = ground
        self.T_hyp = self.T_sed
        
        self.F_DailySum[1] += self.F_Solar[1]
        self.F_DailySum[4] += self.F_Solar[4]
        
        for i in range(len(self.Solar_Blocked[tran])):
            self.Solar_Blocked[tran][i] += veg_block[i]
        self.Solar_Blocked["diffuse"] += veg_block[-1]
        
        # Check if we have interpolation on, and use the appropriate time
        self.T = self.T_bc[time]
        self.T_prev = self.T_bc[time]

    def maccormick2(self, time):
        #===================================================
        #Set control temps
        if not self.prev_km:
            return
        #===================================================
        #Throw away S and mix because we won't need them.

        self.T, S, mix = py_HS.calc_maccormick(self.dt, self.dx, self.U,
                                            self.T_hyp, self.T_prev,
                                            self.Q_hyp,
                                            self.Q_tribs[time],
                                            self.T_tribs[time],
                                            self.prev_km.Q,
                                            self.Delta_T, self.Disp,
                                            True, self.S1,
                                            self.prev_km.T,
                                            self.T, self.next_km.T,
                                            self.Q_accr, self.T_accr,
                                            self.next_km.Mix_T_Delta)

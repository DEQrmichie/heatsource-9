# -*- coding: utf-8 -*-

from ..Stream.PyHeatsource \
    import CalcSolarPosition, GetStreamGeometry, CalcMuskingum, CalcFlows

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
offset = -7
lat = 45.7148
lon = -121.4645
heatsource8 = True

#####################
## Stream Parameters

Q_est = 30
dx = 1
dt = 0.001
S = 7
W_b = 10
W_w = 7
z = 2
n = 0.05
D_est = 5
U = 10
Q = 3
Q_up = 3.01
Q_up_prev = 3.01
Q_bc = 3

Altitude, Zenith, Daytime, tran, Azimuth_mod = CalcSolarPosition(
        lat, lon,
        H, M, S, offset,
        JDC, heatsource8, 0
        )
D_est, A, Pw, Rh, Ww, U, Dispersion = GetStreamGeometry(
        Q_est,
        W_b,
        z,
        n,
        S,
        D_est,
        dx, dt
        )
C1, C2, C3 = CalcMuskingum(Q_est, U, W_w, S, dx, dt)

Q_new, Geom = CalcFlows(Q_est, U, W_w, S, dx, dt, z, n, D_est, Q, Q_up, Q_up_prev, 0,Q_bc)

class TestSolarPosition:
    def test_altitude(self):
        assert round(Altitude,6) == 26.520413
    def test_zenith(self):
        assert round(Zenith,6) == 63.479587
    def test_daytime(self):
        assert Daytime == True
    def test_tran(self):
        assert tran == 2
    def test_Azimuth_mod(self):
        assert round(Azimuth_mod,6) == 137.559547
    
class TestStreamGeometry:
    def test_d_est(self):
        assert D_est == 5
    def test_a(self):
        assert A == 100
    def test_pw(self):
        assert round(Pw,6) == 32.360680
    def test_Rh(self):
        assert round(Rh,6) == 3.09017
    def test_Ww(self):
        assert Ww == 30
    def test_U(self):
        assert U == 0.3
    def test_dispersion(self):
        assert round(Dispersion, 6) == 0.009622
        
class TestMuskingum:
    def test_C1(self):
        assert round(C1,6) == -0.898612
    def test_C2(self):
        assert round(C2,6) == 0.930255
    def test_C3(self):
        assert round(C3,6) == 0.968356
        
class TestFlowCalcs:
    def test_q_new(self):
        assert Q_new == 7
    def test_geom(self):
        assert Geom == 7
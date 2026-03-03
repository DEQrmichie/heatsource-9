# Model updates

## 1.1 Land cover sampling

Land cover physical attributes (height, overhanging distance, canopy cover, and canopy depth) are used as model inputs in Heat Source 9. These data are typically organized and sampled from GIS raster data using TTools. Starting with Heat Source version 7, a star pattern sampling method was introduced where the land cover is sampled at regular intervals along multiple transects radiating outward from a point (or node) at the stream center in a star pattern. Heat Source version 7 and 8 use a star pattern consisting of 7 transects around each stream node in the following directions: northeast, east, southeast, south, southwest, west, and northwest. Heat Source 9 was updated to accommodate a star pattern consisting of any number of transects with each transect separated by an equal number of degrees. This change will allow modeling of solar paths in the southern hemisphere where the sun shines from the north and for special studies that require a higher resolution land cover characterization.

In Heat Source 9, the number of transects, sample points along each transect, and the distance between samples is user defined. For most applications, setting up Heat Source 9 with 8 transects spaced every 45 degrees (northeast, east, southeast, south, southwest, west, northwest, north) will suffice. The measured distance along each transect starts at the stream node. The appropriate number of sample points and distance between samples will depend on the stream width. For most streams, we recommend a maximum sample distance of 8 meters and enough samples to extend about 50 meters from the river bank edge. A smaller sample distance may be necessary for streams with channel widths < 4 meters. A general rule of thumb is the sample distance is half the mean channel width.

## 1.2 Solar radiation attenuation through land cover

Direct beam solar radiation is routed through the transect closest to the solar azimuth. Routing occurs from the outermost land cover sample to the innermost sample.

Shadows that are cast from each land cover sample are calculated as a function of solar altitude and azimuth, in addition to the physical attributes of the land cover and terrain elevation. If the shadow length reaches the stream node, the stream is considered shaded. This methodology is based on simulation of the sun's position and the angle of the solar vector to the stream. When stream surface shade is determined to be occurring, direct beam attenuation occurs as a function of a light extinction coefficient and the path length through the land cover. The path length through the land cover is based on a trigonometric function using the distance between land cover samples, stream aspect, solar altitude, and solar azimuth at each timestep during the model period. Attenuation is calculated for each timestep using the Beer-Lambert law (or often just Beer's law). Direct beam radiant energy that passes through a land cover sample is then routed to the next inner land cover sample and the process is repeated. Once through all land cover samples, the remaining direct beam solar radiation is routed to the stream. Diffuse solar radiation filters through the openings in land cover and is attenuated as a function of canopy opening, estimated with view to sky and canopy cover, or LAI values for each model node. The model applies these equations for each timestep at each node.

Heat Source implements the Beer-Lambert law based on the equations published by Monsi and Saeki (1953) with sub canopy attenuation influenced by methods presented in Norman and Welles (1983), Welles and Cohen (1996), and Chen et al (1998). The principal equations used in Heat Source can be written in two forms (shown as Equation 1 or Equation 2):

$$
\tag{Equation 1}
\phi_{out} = \phi_{in}*(e^{-k*p})
$$

$$
\tag{Equation 2}
\phi_{out} = \phi_{in}*(e^{-K*LAI})
$$

where,

* $\phi_{out}$ = outgoing direct shortwave (w/m2).
* $\phi_{in}$ = incoming direct shortwave (w/m2).
* $k$ = Extinction coefficient used in Equation 1 (m-1).
* $K$ = Extinction coefficient used for LAI in Equation 2 (dimensionless).
* $p$ = solar beam path length through the vegetation (m).
* $LAI$ = Effective Leaf Area Index. The total one-sided area of leaf tissue and other plant material per unit ground surface area (Monsi and Saeki 1953). One half the total green leaf area per unit ground surface area (Chen and Black 1992).

### 1.2.1 Heat Source inputs using LAI and k

If LAI and K are used as inputs into the model, it is possible to calculate k by combining Equation 1 and Equation 2 and solving for k (Equation 3).

$$
\tag{Equation 3}
k=\frac{(LAI*K)}{p}
$$

Because $LAI$ is the product of canopy area density (m-1) and canopy depth (m), canopy depth can be used directly to calculate $k$ as shown in Equation 4.

$$
\tag{Equation 4}
k=\frac{(LAI*K)}{h}
$$

$k$ is then used in Equation 1 to calculate solar attenuation for each transect sample where $p$ represents the path length through the sample.

It may be possible to substitute vegetation height ($h$), measured in meters, as a proxy for canopy depth. When using vegetation height in Equation 4, it is assumed vegetation height is a good proxy for effective canopy depth and that leaf area is distributed throughout that depth on average. If foliage is strongly clumped or concentrated in the upper portion of the crown, $K$ should be adjusted accordingly or the actual canopy depth should be used.

### 1.2.2 Heat Source inputs using canopy cover

If canopy cover is used as an input into the model, Equation 1 is rearranged to solve for $k$ with $p$ from Equation 1 equal to the canopy depth ($h$). The result is what is shown in Equation 5, where $C$ is canopy cover fraction (0 to 1).

$$
\tag{Equation 5}
k=\frac{-\ln(1-C)}{h}
$$

### 1.2.3 Differences between Heat Source version 9 and previous versions of Heat Source

Prior to Heat Source version 7, $h$ in Equation 5 was set to vegetation height, following the approach described by Chen et al (1998). As described earlier, the core assumption of this approach is that canopy is distributed over the full height of the tree. When this assumption is not valid, the canopy cover input must be adjusted.

In Heat Source version 7 and 8, $h$ in Equation 5 is set as a constant at 10 meters (Boyd and Kasper 2003). The rationale for using 10 meters is that it is a better approximation of the canopy depth in larger trees, compared to previous model versions that use the full height of the tree. Assuming canopy distribution over the full height of the tree can sometimes underestimate the riparian extinction coefficient because it does not account for foliage clumping or the airspace between the ground and bottom canopy layers. This can result in less attenuation and less shade. This result may seem counterintuitive as more shade is produced by taller trees with a shallower canopy depth. There is more shade because the extinction coefficient ($k$) is bigger due to the shorter distance a solar ray must travel to equal the same canopy cover. On the flip side, using a constant 10 meters may overestimate the canopy depth for shrubs and underestimate the canopy depth for very tall forests.

Heat Source version 9 adds canopy depth as a model input, allowing spatially variable depths for different classes of vegetation and reducing the need to adjust canopy cover beyond typical field measured values or literature ranges when 10 meter canopy depth is not appropriate.

Canopy depth also provides consistency between LAI based inputs and canopy cover based inputs. Values for LAI and K input into the model are most often obtained from literature or derived from optical measurements using hemispherical photography or field instruments like the LAI-2000/LAI-2200C. Optical methods primarily estimate gap fraction and transmittance, then derive canopy metrics with corrections for clumping and non-uniform canopy structure. The resulting LAI and K values are therefore more representative of the area occupied by foliage without the airspace (Breda 2003; Frazer et al 1997; Welles and Cohen 1996; van Gardingen et al 1999).

Heat Source 9 provides an option in the control file to use the Heat Source version 8 methods for solar calculation. The primary difference is that this option will set the canopy depth path length at a constant 10 meters (Boyd and Kasper 2007) or the vegetation height for an emergent sample. The number of land cover transects is also fixed at 7 transects consistent with the sampling pattern in Heat Source version 8. There are also minor differences in terms of how direct and diffuse fluxes are attenuated. These differences are minor and do not substantially impact shade and temperature predictions. They were retained for consistency. This option is provided so older models can be ported to Heat Source 9 and still maintain consistent solar flux and shade results.

In tests, the maximum difference in stream temperature and effective shade between Heat Source 8 and an equivalent Heat Source 9 run with heatsource8=True was less than 0.01 for both metrics (Table 1). This small difference is due to a minor calculation error in Heat Source 8 in how diffuse solar radiation flux attenuation from topographic features was computed. This issue was fixed in Heat Source 9.

Table 1. Minimum, mean, and maximum absolute differences between Heat Source 8 and an equivalent Heat Source 9 run with heatsource8=True.
| File name |   Min |  Mean |    Max |
|:----------|------:|------:|-------:|
| Heat_Cond | 0.000 | 0.002 |  0.040 |
| Heat_Conv | 0.000 | 0.002 |  0.045 |
| Heat_Evap | 0.000 | 0.003 |  0.090 |
| Heat_Long | 0.000 | 0.002 |  0.039 |
| Heat_SR1  | 0.000 | 0.000 |  0.000 |
| Heat_SR4  | 0.000 | 0.004 | 11.560 |
| Heat_SR6  | 0.000 | 0.003 |  9.530 |
| Hyd_DA    | 0.000 | 0.000 |  0.000 |
| Hyd_DM    | 0.000 | 0.000 |  0.000 |
| Hyd_Disp  | 0.000 | 0.000 |  0.000 |
| Hyd_Flow  | 0.000 | 0.000 |  0.000 |
| Hyd_Hyp   | 0.000 | 0.000 |  0.000 |
| Hyd_Vel   | 0.000 | 0.000 |  0.000 |
| Hyd_WT    | 0.000 | 0.000 |  0.000 |
| Rate_Evap | 0.000 | 0.000 |  0.000 |
| Shade     | 0.000 | 0.000 |  0.002 |
| Temp_H2O  | 0.000 | 0.000 |  0.007 |
| Temp_Sed  | 0.000 | 0.000 |  0.003 |
| VTS       | 0.000 | 0.000 |  0.000 |

The mean absolute difference between Heat Source 8 and an equivalent Heat Source 9 run with heatsource8=False was 0.05 percentage points for effective shade and 0.275 deg-C for temperature. Maximum absolute difference was 0.12 percentage points (or 12%) for effective shade and 1.007 deg-C for temperature, respectively. The canopy depth in Heat Source 9 was set equal to the vegetation height. Using half the vegetation height for canopy depth reduces the maximum absolute difference to 0.068 percentage points for effective shade (6.8%) and 0.178 deg-C for temperature (Table 3).

Table 2. Minimum, mean, and maximum absolute differences between Heat Source 8 and an equivalent Heat Source 9 run with heatsource8=False.
| File name |   Min |   Mean |     Max |
|:----------|------:|-------:|--------:|
| Heat_Cond | 0.000 |  0.950 |   5.020 |
| Heat_Conv | 0.000 |  1.534 |  10.768 |
| Heat_Evap | 0.000 |  3.276 |  32.222 |
| Heat_Long | 0.000 |  1.974 |   9.118 |
| Heat_SR1  | 0.000 |  0.000 |   0.000 |
| Heat_SR4  | 0.000 | 14.770 | 531.046 |
| Heat_SR6  | 0.000 | 12.456 | 463.422 |
| Hyd_DA    | 0.000 |  0.000 |   0.000 |
| Hyd_DM    | 0.000 |  0.000 |   0.000 |
| Hyd_Disp  | 0.000 |  0.000 |   0.000 |
| Hyd_Flow  | 0.000 |  0.000 |   0.000 |
| Hyd_Hyp   | 0.000 |  0.000 |   0.000 |
| Hyd_Vel   | 0.000 |  0.000 |   0.000 |
| Hyd_WT    | 0.000 |  0.000 |   0.000 |
| Rate_Evap | 0.000 |  0.000 |   0.000 |
| Shade     | 0.000 |  0.050 |   0.120 |
| Temp_H2O  | 0.000 |  0.275 |   1.007 |
| Temp_Sed  | 0.000 |  0.210 |   0.606 |
| VTS       | 0.000 |  0.014 |   0.041 |

Table 3. Minimum, mean, and maximum absolute differences between Heat Source 8 and an equivalent Heat Source 9 run with heatsource8=False and canopy depth set to 50% of the vegetation height.
| File name |   Min |  Mean |     Max |
|:----------|------:|------:|--------:|
| Heat_Cond | 0.000 | 0.543 |   3.895 |
| Heat_Conv | 0.000 | 0.813 |   7.680 |
| Heat_Evap | 0.000 | 1.768 |  22.817 |
| Heat_Long | 0.000 | 1.395 |   7.423 |
| Heat_SR1  | 0.000 | 0.000 |   0.000 |
| Heat_SR4  | 0.000 | 9.137 | 530.799 |
| Heat_SR6  | 0.000 | 7.600 | 463.069 |
| Hyd_DA    | 0.000 | 0.000 |   0.000 |
| Hyd_DM    | 0.000 | 0.000 |   0.000 |
| Hyd_Disp  | 0.000 | 0.000 |   0.000 |
| Hyd_Flow  | 0.000 | 0.000 |   0.000 |
| Hyd_Hyp   | 0.000 | 0.000 |   0.000 |
| Hyd_Vel   | 0.000 | 0.000 |   0.000 |
| Hyd_WT    | 0.000 | 0.000 |   0.000 |
| Rate_Evap | 0.000 | 0.000 |   0.000 |
| Shade     | 0.000 | 0.027 |   0.068 |
| Temp_H2O  | 0.000 | 0.143 |   0.718 |
| Temp_Sed  | 0.000 | 0.112 |   0.341 |
| VTS       | 0.000 | 0.014 |   0.041 |

---

# References

Bencala KE, Walters RA. 1983. Simulation of solute transport in a mountain pool-and-riffle stream: a transient storage model. Water Resour Res. 19(3):718-724. doi:10.1029/WR019i003p00718.

Boyd M, Kasper B. 2003. Analytical methods for dynamic open channel heat and mass transfer: methodology for Heat Source model version 7.0. Oregon Department of Environmental Quality. Updated 2007 Feb 20. Available from: http://www.deq.state.or.us/wq/TMDLs/tools.htm.

Breda NJJ. 2003. Ground-based measurements of leaf area index: a review of methods, instruments and current controversies. J Exp Bot. 54(392):2403-2417. doi:10.1093/jxb/erg263.

Chow VT. 1959. Open-channel hydraulics. New York (NY): McGraw-Hill Book Co. 680 p.

Chen DY, Carsel RF, McCutcheon SC, Nutter WL. 1998. Stream temperature simulation of forested riparian areas: I. Watershed-scale model development. J Environ Eng. 124(4):304-315. doi:10.1061/(ASCE)0733-9372(1998)124:4(304).

Chen JM, Black TA. 1992. Defining leaf area index for non-flat leaves. Plant Cell Environ. 15(4):421-429. doi:10.1111/j.1365-3040.1992.tb00992.x.

Frazer GW, Trofymow JA, Lertzman KP. 1997. A method for estimating canopy openness, effective leaf area index, and photosynthetically active photon flux density using hemispherical photography and computerized image analysis techniques. Natural Resources Canada, Canadian Forest Service, Pacific Forestry Centre, Information Report BC-X-373.

Hart DR. 1995. Parameter estimation and stochastic interpretation of the transient storage model for solute transport in streams. Water Resour Res. 31(2):323-328. doi:10.1029/94WR02739.

Kasahara T, Wondzell SM. 2003. Geomorphic controls on hyporheic exchange flow in mountain streams. Water Resour Res. 39(1):1005. doi:10.1029/2002WR001386.

Monsi M, Saeki T. 1953. Uber den Lichtfaktor in den Pflanzengesellschaften und seine Bedeutung fur die Stoffproduktion. Jpn J Bot. 14(1):22-52.

Norman JM, Welles JM. 1983. Radiative transfer in an array of canopies. Agron J. 75(3):481-488.

Oke TR. 1978. Boundary layer climates. London (UK): Methuen & Co. Ltd. 372 p.

Pelletier GJ, Chapra SC, Tao H. 2006. QUAL2Kw: a framework for modeling water quality in streams and rivers using a genetic algorithm for calibration. Environ Model Softw. 21(3):419-425. doi:10.1016/j.envsoft.2005.07.002.

Sinokrot BA, Stefan HG. 1993. Stream temperature dynamics: measurements and modeling. Water Resour Res. 29(7):2299-2312. doi:10.1029/93WR00540.

van Gardingen PR, Jackson GE, Hernandez-Daumas S, Russell G, Sharp L. 1999. Leaf area index estimates obtained for clumped canopies. Agric For Meteorol. 94:243-257. doi:10.1016/S0168-1923(99)00018-0.

Welles JM, Cohen S. 1996. Canopy structure measurement by gap fraction analysis using commercial instrumentation. J Exp Bot. 47(302):1335-1342.

Wondzell SM. 2011. The role of the hyporheic zone across stream networks. Hydrol Process. 25(22):3525-3532. doi:10.1002/hyp.8119.

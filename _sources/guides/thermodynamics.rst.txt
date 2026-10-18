Thermodynamics: total versus specific quantities and cycle states
=================================================================

Absolute temperatures are K, pressure Pa, total volume m^3 and gas amount mol.
The universal constant is 8.31446261815324 J/(mol K). Heat/work helpers use
J unless noted; ideal-gas internal energy and enthalpy use J/kg with cv/cp in
J/(kg K). Do not mix total and specific energy in a balance. Work by the system
is positive: Q=delta_U+W. Isothermal entropy is n*R*ln(V2/V1) [J/K], and
reversible isothermal work equals T*delta_S.

Ideal cycles and limits
-----------------------

Carnot reservoirs require T_hot>T_cold>0. Ideal Otto compression ratio and
Brayton pressure ratio are >=1, with gamma>1 and constant heat capacities.
Diesel requires 1<=cutoff<=compression_ratio. cutoff=1 is the limiting Otto
case, handled analytically; a nearby cutoff uses expm1 to avoid cancellation.
The formulas do not include real gas, variable heat capacity, combustion loss,
pumping loss or component efficiencies.

Rankine numbering is explicit: h4 is condenser outlet/pump inlet, h1 pump
outlet, h2 boiler outlet/turbine inlet, and h3 turbine outlet. All h are J/kg.
Net work=(h2-h3)-(h1-h4), Qin=h2-h1, Qout=h3-h4. Ordering checks exclude
non-power-producing inputs; the enthalpies themselves are not checked, so take
them from steam_properties (below) or a steam table.

.. code-block:: python

    from engineeringkit.thermodynamics import rankine_ideal_efficiency
    eta = rankine_ideal_efficiency(210e3, 3200e3, 2100e3, 200e3)
    # Approximately 0.364548; requires externally supplied consistent enthalpies.

Exergy destruction T0*Sgen uses total entropy generation J/K and produces J;
Sgen must be nonnegative. Heat exergy Q*(1-T0/Tsource) can be negative for a
source below the environment; it is not arbitrarily clipped. Clausius compares
a supplied discretized integral and does not construct a reversible path.

curzon_ahlborn_efficiency, 1-sqrt(Tc/Th), is the efficiency of an
endoreversible Carnot engine at maximum power when heat crosses finite
conductances linearly (Q = K*delta_T). Other heat-transfer laws give other
values; it is a benchmark for real plants, not a bound.

Real gases and vapour pressure
------------------------------

van_der_waals_pressure takes molar constants a [Pa m^6/mol^2] and b
[m^3/mol] with the total volume V [m^3] and amount n [mol]; V must exceed
n*b. With a=b=0 it is the ideal-gas law. Below the critical temperature the
isotherm has an unstable loop; the function returns it as is, and the vapour
pressure needs Maxwell's equal-area construction (see the gallery).
van_der_waals_critical_point returns (Tc, pc, vc) with vc molar. The
critical compressibility is always 3/8, against about 0.27-0.29 for real
gases, so the model is qualitative near the critical point.

joule_thomson_coefficient_vdw is the zero-pressure limit
(2a/(R*T) - b)/cp_molar [K/Pa], with cp_molar molar [J/(mol K)]. Its
inversion temperature 2a/(R*b) is higher than measured (about 850 K against
620 K for nitrogen); at high pressure the real inversion temperature falls.

clausius_clapeyron_pressure assumes an ideal vapour, negligible liquid
volume and constant h_fg [J/kg] with the specific gas constant R_gas
[J/(kg K)], so it is accurate within about 10 K of the reference state.
antoine_vapor_pressure takes A for p in Pa and B, C in K. Tabulated constants
for mmHg and degC convert by adding log10(133.322)=2.12490 to A and
subtracting 273.15 from C. Each constant set is valid only over the range it
was fitted to.

Water and steam (IAPWS-IF97)
----------------------------

steam_properties implements IF97 regions 1 (compressed liquid), 2 (vapour),
3 (near the critical point, above 623.15 K and above the B23 boundary
pressure) and 4 (saturation), for 273.15 K <= T <= 1073.15 K and
p <= 100 MPa. Region 5 (above 1073.15 K) is not implemented and raises
ValueError. Region 3 is a Helmholtz function of density, so the density at
a given pressure is found by Brent's method on the isotherm; below the
critical temperature the liquid or vapour root is taken, by the side of the
saturation pressure the state lies on. Give the pressure and exactly one of T, x, s or h. Specific
enthalpy and entropy are zero for saturated liquid at the triple point, as
in steam tables. A state given by T is liquid at or above the saturation
pressure and vapour below it; a state on the saturation line needs x. Wet
states (x) need 611.213 Pa <= p <= 22.064 MPa (the critical pressure);
above 16.53 MPa the saturated phases are region 3 densities. States given by s or h are found by Brent's
method on temperature, to within 1e-10 K. cp is NaN inside the dome and x is
NaN outside it. The results match the release's verification tables to nine
significant digits.

saturation_pressure and saturation_temperature are the region 4 equations
from 273.15 K (611.213 Pa) to the critical point (647.096 K, 22.064 MPa).
They accept arrays; steam_properties takes scalars.

Moist air
---------

humidity_ratio, relative_humidity and dew_point treat moist air as an ideal
mixture at total pressure p [Pa] (sea level by default), with the saturation
pressure over liquid water from IF97. The humidity ratio is in kg of water
per kg of dry air, relative humidity is a fraction in [0, 1] (not percent),
and the molar mass ratio is 0.621945. Below 0 degC the saturation pressure
over ice is lower; frost points are not covered and dew_point raises
ValueError for them. The ideal-mixture model gives the saturated humidity
ratio at 20 degC within 0.5% of the ASHRAE tables.

Visualizers
-----------

plot_pv_diagram and plot_ts_diagram draw a path given as matching 1-D
arrays, in kPa and in kJ/(kg K) and degC; plot_ts_diagram adds the water
saturation dome (plot_saturation_dome) unless dome=False. The dome is drawn
up to the critical point, 647.096 K.

Validation
----------

Tests integrate P dV independently, compare W to T delta_S, reconstruct the
Otto four-state first-law efficiency, check COP_HP=COP_R+1, and verify Rankine
work/heat accounting plus the Diesel-to-Otto limit. IF97 is checked against
the release's verification tables for regions 1, 2, 4 and the B23 boundary,
against the Gibbs relations (dh/dT)_p = cp, (dh/ds)_p = T and
(dh/dp)_s = v, and against the Clapeyron equation and equal Gibbs energies on
the saturation line. The ideal Rankine cycle at 10 MPa, 500 degC and 10 kPa
reproduces the textbook 40.2%. The van der Waals critical point is checked
for a horizontal inflection, the Joule-Thomson coefficient against an
isenthalpic van der Waals throttle, Curzon-Ahlborn against a numerical
maximization of power, and the psychrometric functions by round trips and
against ASHRAE saturated humidity ratios. These tests validate the equations,
not the accuracy of a real engine prediction.

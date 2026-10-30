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

dual_efficiency(r, pressure_ratio, cutoff) adds heat first at constant volume
(pressure ratio alpha=p3/p2) and then at constant pressure (cutoff ratio);
alpha=1 is Diesel and cutoff=1 is Otto. stirling_efficiency(T_hot, T_cold,
compression_ratio, gamma, regenerator_effectiveness) is the ideal-gas
Stirling cycle; with effectiveness 1 it is Carnot's efficiency whatever the
compression ratio, and with effectiveness eps the external heat includes
(1-eps)*cv*(T_hot-T_cold).

Steady-flow devices
-------------------

steady_flow_work(h_in, h_out, q=, V_in=, V_out=, z_in=, z_out=) is the
one-inlet, one-outlet energy balance per unit mass [J/kg]: work done by the
fluid is positive (turbine), work done on it negative (compressor, pump),
and q is heat added. nozzle_exit_velocity is the adiabatic, work-free case.
turbine_isentropic_efficiency is (h1-h2)/(h1-h2s),
compressor_isentropic_efficiency (h2s-h1)/(h2-h1) (also for pumps), and
nozzle_isentropic_efficiency the ratio of exit kinetic energies;
expansion_outlet_enthalpy and compression_outlet_enthalpy invert them. mixing_chamber takes the mass flows
[kg/s] and enthalpies of the inlets and returns the total flow and the mixed
enthalpy, adiabatic and with no work. throttle_steam(state, p_out) keeps the
enthalpy of a SteamState and returns the state at the lower pressure.

Power plants state by state
---------------------------

rankine_cycle(p_boiler, p_condenser, T_turbine_inlet) solves a steam plant
with IAPWS-IF97. Without T_turbine_inlet the turbine takes saturated vapour,
which needs a subcritical boiler; a supercritical boiler needs a turbine inlet
above 647.096 K. reheat_pressure and T_reheat add one reheater;
feedwater_pressure adds one open (direct-contact) feedwater heater fed by
steam bled from the turbine at that pressure, before or after the reheater.
turbine_efficiency applies to every turbine section between two of the
pressures, from that section's own inlet, and pump_efficiency to both pumps.
The condenser and heater deliver saturated liquid; there are no pressure
losses and no closed heaters. RankineResult holds the labelled states, the
efficiency, the turbine, pump, boiler and condenser energies per kg of boiler
steam [J/kg], the bleed fraction, the back-work ratio and the exit quality;
plot() draws it on the T-s diagram.

brayton_cycle(pressure_ratio, T_min, T_max) is the air-standard gas turbine
with constant cp [J/(kg K)] and gamma. compressor_stages and turbine_stages
split the pressure ratio equally, with intercooling back to T_min and reheat
back to T_max between stages; each stage has the given isentropic
efficiency. A regenerator of effectiveness eps heats the compressed air by
eps*(T_exhaust - T_compressor) when the exhaust is hotter, and does nothing
otherwise. BraytonResult gives the states (T, p and entropy relative to the
compressor inlet), the efficiency, the energies per kg [J/kg] and the
back-work ratio, and plots the cycle on T-s axes. Both solvers are
scalar-only.

Exergy destruction T0*Sgen uses total entropy generation J/K and produces J;
Sgen must be nonnegative. Heat exergy Q*(1-T0/Tsource) can be negative for a
source below the environment; it is not arbitrarily clipped. Clausius compares
a supplied discretized integral and does not construct a reversible path.

curzon_ahlborn_efficiency, 1-sqrt(Tc/Th), is the efficiency of an
endoreversible Carnot engine at maximum power when heat crosses finite
conductances linearly (Q = K*delta_T). Other heat-transfer laws give other
values; it is a benchmark for real plants, not a bound.

Real gases, mixtures and vapour pressure
----------------------------------------

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

compressibility_factor(p, T, v) is p*v/(R*T) with v molar [m^3/mol].
pitzer_acentric_factor takes the vapour pressure at 0.7*Tc and pc [Pa];
pitzer_virial_compressibility is the generalized second-virial correlation
Z = 1 + (B0 + omega*B1)*pr/Tr with Abbott's B0 and B1, valid for gases at
low to moderate reduced pressure (roughly v_r > 2), not liquids.

Dalton's ideal-gas mixtures: mole_fractions and mass_fractions convert with
the molar masses [kg/mol], mixture_molar_mass is sum(y*M) [kg/mol],
partial_pressures is y*p [Pa], and entropy_of_mixing is
-n*R*sum(y*ln y) [J/K] for n mol mixed at equal temperature and pressure.
Fractions are 1-D arrays that must sum to one.

cp_minus_cv(T, v, beta, kappa_T) = T*v*beta^2/kappa_T and
joule_thomson_coefficient(T, v, beta, cp) = v*(T*beta - 1)/cp follow from
Maxwell's relations; v is specific [m^3/kg], beta the volume expansivity
[1/K], kappa_T the isothermal compressibility [1/Pa] and cp [J/(kg K)].
They hold for any fluid; the inputs can come from finite differences of
steam_properties.

clausius_clapeyron_pressure assumes an ideal vapour, negligible liquid
volume and constant h_fg [J/kg] with the specific gas constant R_gas
[J/(kg K)], so it is accurate within about 10 K of the reference state.
antoine_vapor_pressure takes A for p in Pa and B, C in K. Tabulated constants
for mmHg and degC convert by adding log10(133.322)=2.12490 to A and
subtracting 273.15 from C. Each constant set is valid only over the range it
was fitted to.

Water and steam (IAPWS-IF97)
----------------------------

steam_properties implements all five IF97 regions: 1 (compressed liquid),
2 (vapour), 3 (near the critical point, above 623.15 K and above the B23
boundary pressure), 4 (saturation) and 5 (high-temperature steam), for
273.15 K <= T <= 1073.15 K at p <= 100 MPa and 1073.15 K < T <= 2273.15 K
at p <= 50 MPa (the 2007 revision of region 5). Region 3 is a Helmholtz function of density, so the density at
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
All three accept arrays. steam_properties broadcasts p against the given
property and returns a SteamState of arrays (region an integer array),
evaluating each element in turn; scalar calls return floats as before.

Refrigerants (Peng-Robinson)
----------------------------

peng_robinson_pressure(T, v, Tc, pc, omega) takes molar volume [m^3/mol];
peng_robinson_compressibility returns the vapour (largest) or liquid
(smallest) root of the cubic in Z; peng_robinson_saturation_pressure finds
equal liquid and vapour fugacities between the spinodals for
0.45*Tc <= T <= 0.9999*Tc. These three are scalar-only. A CubicFluid holds
Tc [K], pc [Pa], omega, the molar mass [kg/mol] and a constant ideal-gas cp0
[J/(kg K)]; R134A is built in (constants from Çengel and Boles, Tables A-1
and A-2a, and Poling, Prausnitz and O'Connell).

refrigerant_properties(fluid, p, T=, x=, s=, h=) mirrors steam_properties for
scalars: the ideal gas plus the Peng-Robinson departure functions, shifted
to the IIR reference (h = 200 kJ/kg, s = 1 kJ/(kg K) for saturated liquid at
0 degC). Given T, the phase is liquid at or above the saturation pressure;
above Tc it is "supercritical". s and h are inverted by Brent's method
between 0.45*Tc and 2.5*Tc. For R-134a the vapour pressure is within 2% and
the latent heat within 3% of the tables; liquid volumes are 5-15% high, and
a constant cp0 makes compressor outlet temperatures a few kelvin off.

vapor_compression_cycle(p_evaporator, p_condenser, fluid=R134A, superheat=,
subcooling=, compressor_efficiency=) solves Perkins's cycle: saturated (or
superheated) vapour is compressed, condensed to saturated (or subcooled)
liquid and throttled at constant enthalpy. VaporCompressionResult gives the
four states, q_evaporator, q_condenser and compressor_work [J/kg], both COPs,
and plots the cycle on the p-h diagram. Çengel and Boles's Example 11-1
(0.14 and 0.8 MPa) gives COP 3.94 against the tables' 3.97.

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

moist_air_enthalpy(T, w) is 1005*t + w*(2501e3 + 1820*t) [J/kg dry air], t in
degC, zero for dry air at 0 degC; condensate and spray water take saturated
liquid enthalpies from IF97 on the same scale.
adiabatic_saturation_temperature(T, w, p) solves the adiabatic saturator's energy balance
(practically the wet-bulb temperature) and is scalar-only.
adiabatic_mixing(m_a1, T1, w1, m_a2, T2, w2) weights w and h by the dry-air
flows [kg/s]; a relative humidity above one at the outlet means fog.
cooling_coil(T_in, phi_in, T_out, p) returns the outlet humidity ratio, the
heat removed and the condensate per kg of dry air; below the dew point the
air leaves saturated at T_out. cooling_tower(m_water, T_water_in,
T_water_out, T_air_in, phi_air_in, T_air_out, phi_air_out=1, p) returns the
dry-air flow and the make-up water [kg/s] from the water and energy
balances. Fan work and drift are neglected.

Combustion
----------

A fuel is given by its atoms per molecule, c, h and o (CcHhOo), burned
completely in dry air of 3.76 mol N2 per mol O2. fuel_molar_mass [kg/mol],
stoichiometric_air_fuel_ratio [kg air/kg fuel, with M_air = 0.02897 kg/mol]
and combustion_products (mol of CO2, H2O, O2 and N2 per mol of fuel, with
excess_air as a fraction, 0.5 for 150% theoretical air) are stoichiometry
only; rich mixtures are not covered. heating_value(c, h, o, h_f_fuel=,
water=) applies Hess's law with CODATA enthalpies of formation, h_f_fuel in
J/mol for the fuel in the phase it is burned in, and returns J/kg: the lower
heating value with water="vapour", the higher with "liquid".
adiabatic_flame_temperature(c, h, o, h_f_fuel=, excess_air=, T_air=) balances
the products' enthalpy, from NIST's Shomate equations (298.15 K to 6000 K),
against the fuel at 298.15 K and the air at T_air. It neglects dissociation,
so above about 2000 K it is an upper bound (methane in air: 2327 K here
against about 2230 K at equilibrium); it is scalar-only.

Visualizers
-----------

plot_pv_diagram and plot_ts_diagram draw a path given as matching 1-D
arrays, in kPa and in kJ/(kg K) and degC; plot_ts_diagram adds the water
saturation dome (plot_saturation_dome) unless dome=False. The dome is drawn
up to the critical point, 647.096 K.

Validation
----------

Tests reproduce Çengel and Boles's Rankine plants (43.0%, 45.0% and
46.3%), check the Rankine and Brayton energy balances and the Ericsson limit,
derive the Stirling and dual efficiencies from their processes, match the
IF97 region 5 verification table to nine digits, check the Peng-Robinson
critical inflection, Maxwell's equal-area rule, Clapeyron's slope and the
R-134a tables, verify a Maxwell relation and the Joule-Thomson throttle on
IF97 steam, the Dalton mixing entropy against isothermal expansions, the
psychrometric mass and energy balances, Hess heating values (methane LHV
50.0 MJ/kg) and the flame-temperature energy balance, and the Shomate
enthalpies against JANAF. Tests also integrate P dV independently, compare W to T delta_S, reconstruct the
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

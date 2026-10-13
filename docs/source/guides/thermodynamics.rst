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
non-power-producing inputs, but no saturation/property table validation occurs.

.. code-block:: python

    from engineeringkit.thermodynamics import rankine_ideal_efficiency
    eta = rankine_ideal_efficiency(210e3, 3200e3, 2100e3, 200e3)
    # Approximately 0.364548; requires externally supplied consistent enthalpies.

Exergy destruction T0*Sgen uses total entropy generation J/K and produces J;
Sgen must be nonnegative. Heat exergy Q*(1-T0/Tsource) can be negative for a
source below the environment; it is not arbitrarily clipped. Clausius compares
a supplied discretized integral and does not construct a reversible path.

Tests integrate P dV independently, compare W to T delta_S, reconstruct the
Otto four-state first-law efficiency, check COP_HP=COP_R+1, and verify Rankine
work/heat accounting plus the Diesel-to-Otto limit. These tests validate the
ideal equations, not the accuracy of a real engine prediction.

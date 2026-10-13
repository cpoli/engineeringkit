Aerodynamics: incompressible limits and compressible conservation
=================================================================

Density is kg/m^3, speed m/s, area m^2, lift/drag N, pressure Pa, angles rad,
and temperature K. Lift is signed; drag is nonnegative. Thin-airfoil
Cl=2*pi*(alpha-alpha_l0) assumes small angle and attached inviscid incompressible
flow. No stall prediction is included. Finite-wing slopes use an idealized
lifting-line correction; 0<e<=1. The empirical Oswald helper is restricted to
aspect ratio 2..10 to avoid nonsensical extrapolation.

Compressible models
-------------------

A calorically perfect gas with constant gamma>1 and gas constant R [J/(kg K)]
is assumed. Stagnation relations are isentropic and do not conserve total
pressure through a shock. Prandtl-Glauert requires 0<=M<1; it is unreliable
near transonic flow despite the mathematical bound allowing values close to
one. Normal-shock routines require upstream M>=1 (the sonic case is a limit).
Area-Mach returns A/A* for M>0 and does not choose a sub/supersonic inverse branch.

.. code-block:: python

    from engineeringkit.aerodynamics import normal_shock_mach2, normal_shock_pressure_ratio
    M2 = normal_shock_mach2(2)  # sqrt(1/3)
    pressure_ratio = normal_shock_pressure_ratio(2)  # 4.5

The Mach-2 benchmark independently reconstructs density, temperature and speed,
and checks mass flux, momentum flux and total enthalpy across the shock. Tests
also check sonic limits, the infinite-aspect-ratio lift slope, and the stationary
best-glide point of the parabolic drag polar.

Flat-plate skin friction and boundary-layer thickness assume a smooth plate
with zero pressure gradient. Re must be positive; transition location, roughness
and flow separation are not detected by these algebraic helpers.

Reference: `NASA Glenn normal-shock equations <https://www.grc.nasa.gov/www/BGH/normal.html>`_.

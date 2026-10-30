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

Oblique shocks and expansion fans
---------------------------------

``oblique_shock_angle`` solves the theta-beta-M cubic in closed form
(Thompson, 1950) and returns the weak branch unless ``strong=True``. It
requires M1>1 and 0<theta<=theta_max, where ``oblique_shock_max_deflection``
gives theta_max; a steeper wedge raises ``ValueError`` because the real shock
is detached and curved, which no algebraic relation describes.
``oblique_shock_deflection`` is the forward relation, defined for beta from
the Mach angle to 90 degrees. Expansions are isentropic:
``prandtl_meyer_expansion`` solves nu(M2)=nu(M1)+theta with a bracketing root
search, and turning angles beyond the vacuum limit (130.45 degrees for
gamma=1.4) are rejected. Benchmarks check mass, normal momentum, tangential
velocity and total enthalpy across an oblique shock, and the Prandtl-Meyer
function against quadrature of its ODE.

.. code-block:: python

    import numpy as np
    from engineeringkit.aerodynamics import oblique_shock, prandtl_meyer_expansion
    shock = oblique_shock(2.0, np.radians(20.0))  # beta 53.42 deg, M2 1.210
    M2 = prandtl_meyer_expansion(1.5, np.radians(20.0))  # 2.207

Fanno and Rayleigh duct flow
----------------------------

``fanno_flow`` and ``rayleigh_flow`` return a ``DuctFlowState`` of ratios to
the sonic state on the same line, matching the tables in Anderson's *Modern
Compressible Flow*. The Fanno friction parameter uses the **Darcy** friction
factor, fL*/D, which is four times the Fanning value 4fL*/D found in most
compressible-flow tables; pass a Darcy factor (for example from
``fluid_engineering.colebrook_friction_factor``) to ``fanno_choking_length``.
Fanno flow is adiabatic with constant friction factor; Rayleigh flow is
frictionless, with the heat given per unit mass and constant cp. Neither
predicts what happens beyond choking (a shock in the duct or a changed inlet
state). Benchmarks integrate the friction ODE and check mass, momentum and
energy along both lines.

Lifting line and panel method
-----------------------------

``lifting_line`` solves Prandtl's equation with Glauert's Fourier series at
``n_terms`` collocation stations, keeping every harmonic so asymmetric twist
is allowed. Chord, angle of attack and zero-lift angle may be constants or
functions of the spanwise station y in metres. It assumes a straight,
unswept wing of moderate to high aspect ratio in attached flow, so it does
not predict stall, and a rectangular wing converges slowly near the tips.
The elliptic wing is reproduced exactly (e=1). ``panel_method`` is the
Hess-Smith constant-source, uniform-vortex method. Points run clockwise from
the lower trailing edge, and an ``Airfoil`` is ordered for you. It is
inviscid and incompressible, so it gives no drag, and with a blunt
(open) NACA trailing edge the circulation lift and pressure-integrated lift
differ by up to about 1%. The circular cylinder is matched to round-off.

Pitching moment and stability
-----------------------------

``naca4_moment_coefficient`` is the thin-airfoil moment about the quarter
chord (nose up positive), independent of angle of attack;
``aerodynamic_center`` gives x_ac/c = 1/4 - m0/a0 from measured slopes.
Stability helpers follow Anderson's *Introduction to Flight*, ch. 7, with
positions h as fractions of the mean aerodynamic chord: ``neutral_point``
uses the wing-body lift slope a_wb, tail efficiency eta in (0, 1] and
downwash gradient in [0, 1). ``pitching_moment_slope`` is a_wb*(h - h_n), and
``static_margin`` is positive when stable. ``downwash_gradient`` is the
far-wake elliptic estimate 2a/(pi AR), which overstates the downwash at a
tail close behind the wing. Only stick-fixed static stability is modelled.

Performance
-----------

``rate_of_climb``, ``climb_angle``, ``load_factor``, ``turn_radius`` and
``turn_rate`` are point-mass force balances for steady climbs and level
coordinated turns, with load factor n>1 for a turn.
``takeoff_ground_roll`` takes lift-off at 1.2 times the stall speed with a
constant mean net force T - [D + mu(W - L)], evaluated by the caller at about
0.7 V_LO. ``breguet_range_propeller`` uses power-specific fuel consumption in
kg/(W s) (0.45 lb/(hp h) is about 7.6e-8) and is independent of speed. The
benchmarks integrate Newton's law over the ground roll and time-step the
cruise.

Flat-plate skin friction and boundary-layer thickness assume a smooth plate
with zero pressure gradient. Re must be positive; transition location, roughness
and flow separation are not detected by these algebraic helpers.

Reference: `NASA Glenn normal-shock equations <https://www.grc.nasa.gov/www/BGH/normal.html>`_.

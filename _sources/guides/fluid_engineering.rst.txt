Fluid engineering: flow regimes and signed quantities
=====================================================

Density is kg/m^3, dynamic viscosity Pa s, velocity m/s, pressure Pa and all
lengths/head m. Reynolds number represents speed: the function returns the
magnitude rho*abs(v)*D/mu. Pressure changes, Poiseuille flow and water-hammer
velocity changes may be signed. Bernoulli assumes steady incompressible flow
along a streamline with no losses, pump or turbine work.

Internal flow
-------------

Poiseuille uses a circular pipe, Newtonian fluid, fully developed laminar flow:
Q=pi*r^4*delta_p/(8*mu*L). The Darcy factor is 64/Re below Re=2300; this is not
the Fanning factor (one quarter as large). Colebrook uses the laminar branch
below 2300 and a turbulent branch from 4000; transition is rejected explicitly.
Relative roughness is epsilon/D, bounded here to [0,0.05]. Swamee-Jain is
restricted to Re between 5000 and 1e8 for this implementation. These numeric
contracts are model scope choices, not universal sharp physical transitions.

.. code-block:: python

    from engineeringkit.fluid_engineering import colebrook_friction_factor
    f = colebrook_friction_factor(1e5, 1e-4)
    # f is approximately 0.018514, dimensionless Darcy friction factor.

The tests integrate a parabolic velocity profile to recover Q, check that
Poiseuille pressure drop equals Darcy loss in the laminar regime, and substitute
computed f into the implicit Colebrook equation across multiple regimes.

The friction factors accept arrays and broadcast Re against the relative
roughness; an array touching the 2300-4000 transition is rejected as a whole.
smooth_pipe_friction_factor (Prandtl-Karman) needs Re >= 4000.
nikuradse_friction_factor is the fully rough limit for sand-grain roughness
in (0, 0.05]; it applies once Re*(k/D)*sqrt(f/8) exceeds about 70.
Swamee-Jain is within 3% of Colebrook in its range, usually within 1%.
pipe_resistance writes Darcy-Weisbach as h = r*Q*abs(Q) with f held fixed.

Hydrostatics
------------

hydrostatic_force_plane is the gauge-pressure resultant on a plane of any
inclination: rho*g*(depth of centroid)*area. center_of_pressure_depth
measures distances along the plane from the free-surface line (depths for a
vertical surface) and needs the second moment about the horizontal
centroidal axis.

Measurement and pumps
---------------------

Venturi assumes horizontal incompressible flow, positive upstream-minus-
throat pressure difference, and 0<A2<A1; the optional discharge coefficient
in (0, 1] defaults to the ideal meter. Pitot uses positive
stagnation-minus-static pressure and is within 1% in air below about Mach
0.28. Cavitation number may be negative when p<p_vapor; nonzero speed is
required. Pump helpers are forward pumping models with nonnegative volume
flow/head and 0<efficiency<=1. Pump specific speed is the dimensionless
angular-speed definition using rad/s; it differs from dimensional industry
conventions based on rpm or feet. npsh_available takes the absolute static
pressure at the inlet; thoma_sigma is NPSH over machine head.

Similarity, channels, transients and networks
---------------------------------------------

buckingham_pi takes each variable's dimensions as a mapping of base
dimension to integer exponent. The repeating variables must number the
rank of the dimension matrix and be independent; each group then contains
one other variable to a positive power. The exponents are exact integers
from rational row reduction. Euler number here is dp/(rho*v^2), without
the factor 1/2 some texts use.

manning_velocity uses the SI coefficient 1 (US customary units need 1.486)
and assumes uniform, fully rough flow. pipe_wave_speed is Korteweg's
thin-wall formula; c1 = 1 for expansion joints, 1 - nu^2 for an anchored
pipe, 0 for a rigid pipe. water_hammer_pressure holds for closures shorter
than critical_closure_time = 2L/a and neglects friction.

hardy_cross needs initial flows that already satisfy continuity at every
junction and a loop matrix with entries -1, 0, +1. It sweeps the loops in
order (Gauss-Seidel), stops when the largest correction is below tol, and
raises RuntimeError if max_iter sweeps are not enough. Convergence is
linear, so loop residuals are of order the last correction times the loop
resistance, not zero.

The rebuild's tests integrate pressure over a submerged wall, integrate the
Poiseuille profile, check the smooth and rough limits of Colebrook, the
Chezy form of Manning, Korteweg's rigid-pipe limit, the agreement of NPSH
and cavitation number, Hardy Cross against a direct Newton solve, and that
every Pi group lies in the null space of the dimension matrix.

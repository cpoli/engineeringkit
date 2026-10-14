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

Measurement and pumps
---------------------

Venturi assumes ideal horizontal incompressible flow, positive upstream-minus-
throat pressure difference, and 0<A2<A1; no discharge coefficient is included.
Pitot uses positive stagnation-minus-static pressure. Cavitation number may be
negative when p<p_vapor; nonzero speed is required. Pump helpers are forward
pumping models with nonnegative volume flow/head and 0<efficiency<=1. Pump
specific speed is the dimensionless angular-speed definition using rad/s;
it differs from dimensional industry conventions based on rpm or feet.

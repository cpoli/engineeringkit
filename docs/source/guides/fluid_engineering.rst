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

Momentum: jets, vanes and bends
-------------------------------

jet_force is the force of a free jet on one vane that turns it through
angle (radians, 0 to pi): rho*A*(V-u)^2*(1-cos(angle)) along the jet and
rho*A*(V-u)^2*sin(angle) across it. The jet stays at atmospheric pressure,
gravity and vane friction are neglected, and the vane may not outrun the
jet. pelton_wheel_power is the wheel form, rho*Q*u*(V-u)*(1-cos(angle)):
a wheel always has a bucket in the jet, so the whole flow Q is turned, and
the power peaks at u = V/2.

bend_force is the force of the fluid on a horizontal (reducing) bend, the
force its anchor must resist: flow enters along +x and leaves at angle
counterclockwise from it. Pressures are gauge, because the atmosphere acts
on the outside of the bend as well; velocities are uniform at each end;
the weight of pipe and water is vertical and not included. Take p2 from
bernoulli_pressure less the bend loss.

.. code-block:: python

    from engineeringkit.fluid_engineering import bend_force
    Fx, Fy = bend_force(998.0, 0.08, 0.0314, 0.0154, 2.0e5, 1.9e5, 1.5708)
    # Fx, Fy in N: the force of the water on a 90-degree reducing bend.

Pump and system curves
----------------------

pump_curve_head is the quadratic fit H = H0 + b*Q - a*Q^2 of one pump and
system_curve_head is H = Hs + k*Q^2, with k from pipe_resistance and the
minor losses (fixed friction factor, exact for fully rough flow).
pump_operating_point solves the two for identical pumps, series of them in
series and parallel such strings in parallel, and returns an
OperatingPoint with summary() and plot(). It raises ValueError if the
static head is not below the combined shut-off head, so no forward flow is
possible. affinity_laws scales flow, head and power as N*D^3, N^2*D^2 and
N^3*D^5; it holds between geometrically similar pumps at equal flow
coefficient, so impeller trimming is only approximate.
euler_turbomachine_head is (u2*Vt2 - u1*Vt1)/g, the ideal head before
slip, friction and shock losses; negative for a turbine.

Orifices and weirs
------------------

orifice_flow is Torricelli's Cd*A*sqrt(2*g*h) for a small orifice with
negligible approach velocity; the default Cd = 0.61 is a sharp-edged
orifice and Cd = 1 the ideal jet. An orifice plate in a pipe, where the
approach velocity matters, is venturi_flow with Cd about 0.6.
tank_drain_time integrates the level of a tank of constant section
quasi-steadily. rectangular_weir_flow (Cd*2/3*b*sqrt(2g)*H^1.5) and
v_notch_weir_flow (Cd*8/15*sqrt(2g)*tan(angle/2)*H^2.5) assume a
sharp crest, a ventilated nappe and the head read upstream of the
drawdown; rehbock_discharge_coefficient (0.611 + 0.075*H/P) carries the
approach velocity into Cd for a full-width rectangular weir.

Specific energy and the hydraulic jump
--------------------------------------

The open-channel energy functions work per unit width (q = Q/b) in a wide
rectangular channel with hydrostatic pressure. specific_energy is
y + q^2/(2*g*y^2); critical_depth is (q^2/g)^(1/3), where E is least
(3*yc/2) and the Froude number is 1. alternate_depths takes scalars only
and returns the subcritical and supercritical depths for one E, raising
ValueError below the minimum. sequent_depth is Bélanger's momentum
relation for a horizontal channel with bed friction over the jump
neglected; it is symmetric, so a jump forms only when the upstream Froude
number exceeds 1. hydraulic_jump_head_loss is (y2-y1)^3/(4*y1*y2) and
holds only for a sequent pair.

Viscous flow and drag
---------------------

couette_velocity and couette_flow_rate are the exact plane
Couette-Poiseuille solution for steady, fully developed laminar flow of a
Newtonian fluid, the plate at y = h sliding at U; backflow appears at the
fixed wall once dp_dx exceeds 2*mu*U/h^2. couette_shear_stress is mu*U/h,
for pure Couette flow. couette_viscometer_torque is the exact circular
Couette torque with end effects neglected; it holds for laminar flow,
which rotating the outer cylinder preserves to much higher speeds than
rotating the inner one (Taylor vortices).

stokes_drag (3*pi*mu*d*V) and stokes_settling_velocity hold for Re below
about 1. sphere_drag_coefficient is White's fit
24/Re + 6/(1+sqrt(Re)) + 0.4, within about 10% of data up to the drag
crisis; Re above 2e5 is rejected rather than extrapolated.
terminal_velocity (scalars only) solves C_D(Re)*Re^2 =
4/3*(rho_p-rho)*rho*g*d^3/mu^2 by Brent's method, which has one root
because the left side grows with Re, and raises ValueError beyond the
drag crisis or for a sphere that does not sink.

The tests for these models check the bend and jet forces against an
independent vector momentum balance, the Pelton optimum at half the jet
speed, Euler's head against shaft power, the operating point against a
root-find on both curves, the affinity laws against the dimensionless
coefficients, the minimum of specific energy, momentum conservation and
energy loss across the jump, the drain time against an ODE integration,
the weir formulas against integrals of Torricelli's speed, terminal
velocity against the force balance and Stokes's limit, and the Couette
profile against the momentum equation and the circular solution.

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

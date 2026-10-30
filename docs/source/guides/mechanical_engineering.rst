Mechanical engineering: power, vibration and mechanism geometry
===============================================================

Mass is kg, lengths m, force N, torque N m, inertia kg m^2 and angular speed
rad/s. Gear ratios use positive integer tooth counts; output speed preserves
the supplied sign convention and does not automatically flip rotation direction.
Output torque includes efficiency in [0,1], so torque*speed equals input power
times efficiency. Bearing L10 returns revolutions at 90% survival, not hours;
the exponent is typically 3 for balls and 10/3 for rollers.

Vibration
---------

The passive oscillator solves m*x''+c*x'+k*x=F(t), m>0, c>=0, k>=0. The
frequency helpers require k>0. damped_frequency returns zero for critical/
overdamping, indicating no oscillatory frequency rather than a decay rate.
Transmissibility is the harmonic force/base-motion amplitude ratio for this
linear SDOF model; exact undamped resonance is singular. No nonlinear stiffness,
friction, actuator saturation or contact is modeled.

.. code-block:: python

    import numpy as np
    from engineeringkit.mechanical_engineering import simulate_sdof
    result = simulate_sdof(1, 4, 4, x0=1, t=np.linspace(0, 5, 101))
    # Critical damping: x=(1+2*t)*exp(-2*t) [m].
    fig, ax = result.plot()

Tests check the analytic critical response, undamped energy conservation,
monotonic dissipative energy, constant-force equilibrium and gear power balance.
Solver time grids must strictly increase; adaptive integration errors are visible.

Mechanisms
----------

Slider-crank displacement uses the positive connecting-rod square-root branch,
r>0, l>=r. Four-bar Grashof compares the shortest+longest links against the
other pair; it does not determine which grounded mechanism can fully rotate.
Square-thread screw torque excludes collar friction; before applying it check
that friction angle plus lead angle stays below pi/2. Models assume rigid links
and ideal joints unless efficiency is explicitly supplied.

Belts and lubrication
---------------------

belt_tension_ratio is the Euler-Eytelwein limit exp(mu*theta/sin(beta/2)) at
the point of slip, for a massless, inextensible belt; centrifugal tension and
belt creep are ignored. Omit groove_angle for a flat belt or a rope on a drum.
petroff_friction_torque assumes a centred journal and a linear shear profile
across a gap c much smaller than r; it overestimates the exact concentric
Couette torque by about 1.5*c/r. sommerfeld_number uses N = omega/(2*pi) in
rev/s and the projected pressure P = W/(2*r*L).
long_bearing_eccentricity inverts the full Sommerfeld solution for an
infinitely long bearing. That solution keeps the negative film pressures a
real film cavitates away, so it puts the journal at 90 degrees to the load and
underestimates the eccentricity of short bearings. No finite-length,
short-bearing or thermal correction is applied.

Linkages and robot arms
-----------------------

fourbar_output_angle grounds link d along the x axis, with the input crank a
at the origin and the output link c at (d, 0). It raises ValueError if the
loop cannot close at any requested angle, so trying a full turn tests
whether the crank can rotate fully. crossed=False is the open circuit.
freudenstein_coefficients fits K1..K3 exactly for three precision points and
by least squares for more. freudenstein_link_lengths rejects nonpositive K1
or K2 rather than reinterpreting a negative length as a link turned through
pi. dh_transform uses the original (distal) Denavit-Hartenberg convention
Rz(theta) Tz(d) Tx(a) Rx(alpha); the modified (Craig) convention is not
provided.

Absorbers, torsion and rotors
-----------------------------

tuned_mass_damper follows Den Hartog: the absorber damping ratio is
c2/(2*m2*omega1), with omega1 the primary natural frequency, not the
absorber's. The default tuning f = 1/(1+mu) and damping
sqrt(3*mu/(8*(1+mu)^3)) make the fixed points equal and put the peaks within
about 1% of sqrt(1+2/mu); this is the force-excited, undamped-primary
optimum. holzer_frequencies treats a free-free chain of rigid disks on
massless shafts, excludes the zero rigid-body mode, and finds only roots that
the trial grid separates. jeffcott_response gives the steady synchronous
whirl of a disk midway on a massless shaft with external viscous damping;
gyroscopic moments, internal damping, bearing flexibility and run-up
transients are not modeled.

Forced vibration and several degrees of freedom
-----------------------------------------------

magnification_factor and phase_angle give the steady state of
m*x''+c*x'+k*x=F0*cos(w*t): amplitude (F0/k)*M and lag phi in [0, pi], with
r=w/wn and zeta=c/(2*sqrt(k*m)). The transient is not included; use
simulate_sdof for run-up. normal_modes solves K*phi=w^2*M*phi with
scipy.linalg.eigh for symmetric M (positive definite) and K (positive
semidefinite), returns every mode including zero-frequency rigid-body modes,
and scales the shapes so phi^T*M*phi=I. modal_response superposes the modes
in closed form with modal (proportional) damping ratios in [0, 1); a general,
non-proportional damping matrix couples the modes and is not handled there.
harmonic_response solves (K-w^2*M+i*w*C)*X=F directly and raises ValueError
at an undamped resonance. rayleigh_frequency is an upper bound on the lowest
frequency for any trial shape; dunkerley_frequency, from the masses and the
diagonal of the flexibility (influence-coefficient) matrix, is a lower bound.
Neither includes the shaft's own mass unless it is lumped into the masses.

Lagrange's equations
--------------------

lagrange_accelerations and simulate_lagrangian take a Python function
L(q, qdot) = T - V in joules and the generalised non-conservative forces
Q(t, q, qdot), and form H*q''=dL/dq-J*qdot+Q by central differences (H the
Hessian of L in qdot, J the mixed derivative). L must be smooth and its
Hessian in qdot nonsingular, i.e. the coordinates independent;
time-dependent Lagrangians and constraints with multipliers are not
supported. The finite differences cost O(n^2) evaluations of L per step and
carry a relative error of about step^2 (1e-8 by default); this is a teaching
tool, not a multibody code.

Gears, shafts, springs and bolts
--------------------------------

lewis_bending_stress is W_t*K_v/(F*m*Y) in SI form (module m in metres),
with the whole load at the tip of one tooth; lewis_form_factor interpolates
Shigley's Table 14-2 for 20-degree full-depth teeth from 12 to 400 teeth.
barth_velocity_factor uses the SI constants of Shigley's equation 14-4. The
AGMA geometry, overload, size, load-distribution and rim factors are not
modelled. Shaft functions assume a solid round section and evaluate the
surface fibre: shaft_diameter uses von Mises (default) or Tresca with the
tensile yield strength; shaft_diameter_goodman is the DE-Goodman equation
with a fully corrected endurance limit and fatigue stress-concentration
factors supplied by the caller, and does not check first-cycle yield. Spring
functions are for round-wire, close-coiled compression springs: spring_rate
counts only wire torsion and spring_shear_stress's corrections are Wahl
(default), Bergsträsser, direct shear only or none. Buckling, surge and end
coil effects are not modelled. Bolted joints: bolt_stiffness puts the shank
and threads in series; member_stiffness is the frustum model for two
identical members of one material (washer face 1.5*d and 30 degrees by
default); bolted_joint shares the load linearly until separation, after
which the bolt carries the whole load. preload_torque's coefficient K
scatters by about 25%.

Linkage motion, cams and balancing
----------------------------------

fourbar_kinematics differentiates the loop closure analytically for the
coupler and output angular velocities and accelerations of a crank with
angular velocity omega2 and acceleration alpha2; at a change point, where
coupler and output are collinear, the solution is infinite.
slider_crank_kinematics is exact for an in-line slider-crank with l>r; the
rod angle is measured from the line of stroke, positive with the crank above
it. cam_follower_motion returns the rise from 0 to L over beta for the
uniform, parabolic, harmonic, cycloidal, 3-4-5 and 4-5-6-7 laws, as
derivatives with respect to cam angle (multiply by omega^k for time
derivatives); the return is L - s, and the cam profile, pressure angle and
follower dynamics are not modelled. balance_rotating_masses balances a rigid
rotor statically in one plane, or force and couple in two planes; flexible
rotors above a bending critical speed need modal balancing, which is not
provided.

Tests and gallery examples check the involute's normals, the capstan limit, the
Reynolds-equation load and torque against Sommerfeld, Petroff against exact
Couette flow, four-bar loop closure on both circuits, Freudenstein synthesis
round trips, DH products against hand trigonometry, Den Hartog's fixed points,
Holzer against the matrix eigenvalue problem and the Jeffcott whirl against a
simulated orbit. The curriculum additions are checked against the simulated
forced steady state, mass and stiffness orthogonality, direct integration of
the MDOF equations, the hand-derived double pendulum and its energy, the
Rayleigh-Dunkerley bracket, the inscribed-parabola cantilever, the Goodman
line, Castigliano's spring deflection, the integrated pressure cone, joint
compatibility, the residual force and couple after balancing, finite
differences of the linkage positions and of the cam laws, and the instant
centres.

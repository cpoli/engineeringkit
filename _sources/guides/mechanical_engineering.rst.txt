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

Tests and gallery examples check the involute's normals, the capstan limit, the
Reynolds-equation load and torque against Sommerfeld, Petroff against exact
Couette flow, four-bar loop closure on both circuits, Freudenstein synthesis
round trips, DH products against hand trigonometry, Den Hartog's fixed points,
Holzer against the matrix eigenvalue problem and the Jeffcott whirl against a
simulated orbit.

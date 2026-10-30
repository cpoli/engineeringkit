Electrical engineering: passive components and circuit responses
================================================================

Voltages are V, current is A, resistance is ohm, capacitance is F, inductance
is H and angular frequency is rad/s. RMS phasors are assumed for power helpers.
Negative voltage/current are signed quantities; passive resistances are
nonnegative, and denominators used to infer current/resistance must be nonzero.
Parallel-resistance helpers require strictly positive branches; represent an
ideal short circuit separately rather than passing a zero branch.

RLC current convention
-----------------------

SeriesRLC models current for an applied voltage. The transfer function is
I(s)/V(s)=s/(L*s^2+R*s+1/C), with zero initial energy. The previous numerator
1 represented charge, despite the response being named current. That error is
corrected here. The response to a 1 V step approaches zero current for R>0 and
has integrated charge C. User-supplied transient grids must be uniformly spaced
and strictly increasing, as required by SciPy step response. Lossless R=0 is allowed and oscillates without decay.

.. code-block:: python

    import numpy as np
    from engineeringkit.electrical_engineering import SeriesRLC
    result = SeriesRLC(2, 1, .25).step_response(np.linspace(0, 20, 20001))
    fig, ax = result.plot()  # current [A]

For the underdamped case the independent benchmark uses
exp(-R*t/(2L))*sin(wd*t)/(L*wd), wd^2=1/(LC)-(R/(2L))^2. Tests also integrate
current to obtain C, verify resonance produces purely resistive impedance,
and check source energy equals capacitor storage plus resistor dissipation
for RC charging. Capacitor impedance excludes zero frequency, where the
ideal capacitor is an open circuit.

Thévenin/Norton helpers transform supplied equivalent parameters; to reduce
a network, find its open-circuit voltage and short-circuit current with
``nodal_analysis``. Power factor is bounded by [-1,1], allowing signed real
power. Balanced sinusoidal operation is assumed for three-phase power; the
result does not cover unbalanced phases or harmonic distortion (use
symmetrical components for the former).

Nodal analysis
--------------

``nodal_analysis`` solves a DC network of resistors, ideal current sources
and ideal voltage sources by modified nodal analysis. Node 0 is ground and
non-ground nodes are numbered 1..n. A current source ``(a, b, I)`` drives
``I`` from node ``a`` through itself into node ``b``; a voltage source
``(p, m, V)`` holds ``v_p - v_m = V``, and its reported current is the one
leaving its + terminal into the circuit. A floating node or a loop of
voltage sources makes the system singular; the function rejects it when the
condition number exceeds 1e12 rather than returning a meaningless answer.
The result's ``kcl_residual`` is the largest net current leaving any node,
an independent check of the solve.

Mesh analysis
-------------

``mesh_analysis`` is the dual of nodal analysis for a planar circuit of
resistors and ideal voltage sources. Meshes (the windows of the drawn
circuit) are numbered 1..n and carry clockwise currents; 0 stands for the
outside. A resistor ``(a, b, R)`` lies on the boundary between meshes ``a``
and ``b`` (``b = 0`` on the outer edge), and a voltage source ``(a, b, V)``
drives current clockwise around mesh ``a`` and therefore anticlockwise
around mesh ``b``. Resistor currents are reported in the clockwise sense of
their first mesh, ``i_a - i_b``. Current sources are not accepted (they
would need supermeshes); convert them to Thévenin form or use
``nodal_analysis``. The tests solve the same ladder both ways.

Second-order transients
-----------------------

``rl_current`` is Helmholtz's first-order rise or decay through ``L/R``.
``rlc_damping`` returns the neper frequency alpha = R/2L, the resonant
frequency omega0 = 1/sqrt(LC) [rad/s] and the regime of a *series* RLC
circuit (for a parallel circuit alpha = 1/2RC instead). The comparison
alpha = omega0 is made within a relative 1e-9, inside which the critical
formula is used, so circuits a hair away from critical damping do not lose
accuracy to cancellation. ``rlc_step_voltage`` is the capacitor voltage from
rest; the loop current is C dv/dt, which ``SeriesRLC.step_response``
returns, and the benchmark tests check both against ``scipy.signal.step``.
R, L and C are scalar-only.

Electronics
-----------

The op-amp helpers are ideal: infinite gain and input resistance, zero
output resistance, no rail limits, offset or bandwidth.
``op_amp_integrator`` integrates a sampled input by the trapezoidal rule, so
its accuracy depends on the sampling. ``shockley_diode_current`` covers
forward and reverse bias but not breakdown, series resistance or high
injection; it refuses voltages above 40 n V_T, where the current would be
about e^40 times I_S and no real diode is still described. The
temperature enters only through V_T: I_S is an input, though in a real
diode it grows steeply with temperature. ``diode_operating_point`` brackets
the load-line crossing and solves it with Brent's method: on the current for
forward bias, where the diode voltage n V_T ln(1 + i/I_S) cannot overflow,
and on the voltage for reverse bias. A Lambert-W closed form checks it.

Machines
--------

``dc_motor_torque`` and ``dc_motor_speed`` describe a separately excited
(or permanent-magnet) DC machine with constant flux k*Phi; series and
shunt field interaction, armature reaction and brush drop are left out.
``induction_motor_torque`` uses the per-phase T circuit (stator R1 + jX1,
magnetizing jXm, rotor R2/s + jX2, all referred to the stator) without a
core-loss resistance; the synchronous speed is the *mechanical* speed
2 pi f/(poles/2). The torque is computed multiplied through by s^2, so slip
0 gives zero torque instead of a division by zero, negative slip is
generating and slip above 1 is plugging. Deep-bar and double-cage rotors,
whose resistance changes with slip, are not modelled.

Transmission lines
------------------

``transmission_line_z0`` and ``transmission_line_velocity`` are lossless
models using H/m and F/m. ``propagation_constant`` keeps R [ohm/m] and
G [S/m] and returns gamma = alpha + j beta [1/m]; ``distortionless_condition``
tests Heaviside's R/L = G/C. ``reflection_coefficient``,
``standing_wave_ratio`` and ``input_impedance`` assume a lossless line of
real characteristic impedance and a passive, finite load (an open circuit
is the limit Gamma = 1 and cannot be passed as an impedance). The
electrical length is beta*l in radians, so a quarter wavelength is pi/2.
Where a line transforms the load into an open circuit the input impedance
is numerically very large.

Power systems
-------------

``symmetrical_components`` and ``phase_components`` use the Fortescue
operator a = exp(j 2 pi/3) with the phase-a members of each sequence set;
they accept volts, amps or per-unit values alike. ``admittance_matrix``
builds the bus admittance matrix from pi-model lines (series impedance and
total line charging, per unit). ``gauss_seidel_load_flow`` treats bus 0 as
the slack bus and every other bus as a PQ bus with the given consumed power
(generation is negative load); PV buses, tap-changing transformers and
reactive limits are not modelled. It does not raise on non-convergence:
check ``converged`` and ``history``. An acceleration factor of 1.4-1.6
usually cuts the number of sweeps, but too large a value slows or prevents
convergence. Tests check the two-bus case against its closed form and the
example checks a five-bus network against ``scipy.optimize.root``.

``newton_raphson_load_flow`` solves the same problem in polar coordinates
and adds PV buses, ``{bus: |V|}``: their real power comes from ``p_load``
(negative for a generator), their ``q_load`` entry is ignored and their
reactive output is reported in ``power``. Its ``history`` is the largest
power mismatch before each step, so the last entry is below ``tol`` on
convergence. The Jacobian is dense, which is fine for teaching networks of a
few dozen buses. Generator reactive limits are not enforced. The tests check
quadratic convergence and agreement with Gauss-Seidel.

``per_unit_base`` returns the base impedance and current for a base power
(three-phase total) and a base voltage (line to line); ``change_base`` moves
a nameplate impedance to the system base. ``fault_currents`` computes the
phase currents of three-phase, single-line-to-ground (phase a),
line-to-line (b-c) and double-line-to-ground faults from the Thévenin
sequence impedances at the bus, with optional fault impedance; prefault load
current is neglected. Use one unit system, per unit or volts and ohms,
throughout.

Transient stability
-------------------

``swing_response`` integrates the classical single-machine, infinite-bus
swing equation (constant voltage behind transient reactance, constant
mechanical power) with ``integrate_ode``, switching from the fault-on to the
post-fault power-angle curve at the clearing time. ``H`` is the inertia
constant on the machine base and powers are per unit on the same base. The
``stable`` flag is a first-swing check: the angle stays below the post-fault
unstable equilibrium pi - arcsin(Pm/P3). Later swings, governor and
excitation action and multi-machine effects are outside the model.
``equal_area_clearing_angle`` returns the prefault angle and the critical
clearing angle; when the fault-on curve keeps the rotor from ever passing
the unstable equilibrium it returns that equilibrium, and when the
post-fault network cannot hold the rotor at all it raises ``ValueError``.
``critical_clearing_time`` has a closed form only for a fault that drops the
electrical power to zero; otherwise bisect on ``swing_response``, as the
example does. Tests check equal areas, energy conservation after clearing
and that the critical time separates stable from unstable swings.

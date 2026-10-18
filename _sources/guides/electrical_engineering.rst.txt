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

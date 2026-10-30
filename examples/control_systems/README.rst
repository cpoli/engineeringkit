Examples
========

.. include:: /_generated/nav/control_systems.rst

This gallery walks through every public feature of
``engineeringkit.control_systems``, one breakthrough per script: feedback
and PID control, algebraic and Lyapunov stability tests, frequency-response
design and sensitivity, transient specifications, the root locus,
sampled-data control, dead time, nonlinear loops, state-space design and
optimal control and estimation.

:doc:`Units, assumptions and numerical validation </guides/control_systems>`
gives the validity limits. Each script is self-contained and can be run
directly with ``python examples/control_systems/<section>/<script>.py``.

Sections
--------

- **feedback** -- Watt's governor, Minorsky's PID, Black's feedback amplifier and Ziegler-Nichols tuning.
- **stability** -- Maxwell's characteristic roots, Routh's array, Lyapunov's direct method and Hurwitz's determinants.
- **frequency** -- Nyquist's criterion, Bode's stability margins over every crossover, Bode's sensitivity integral and the Nichols chart.
- **time_response** -- Hazen's second-order servomechanism and Smith's predictor for dead time.
- **root_locus** -- Evans's root locus.
- **sampled** -- Ragazzini and Zadeh's zero-order-hold discretisation, Jury's stability test and Tustin's bilinear transform.
- **nonlinear** -- Andronov's phase plane and Kochenburger's describing function.
- **state_space** -- Kalman's controllability, Luenberger's observer and Wonham's pole placement.
- **optimal** -- Kalman's LQR and its discrete form, the Kalman-Bucy filter and Doyle's LQG counterexample.

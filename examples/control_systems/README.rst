Examples
========

.. include:: /_generated/nav/control_systems.rst

This gallery walks through every public feature of
``engineeringkit.control_systems``, one breakthrough per script: feedback
and PID control, algebraic and Lyapunov stability tests, frequency-response
design, transient specifications, the root locus, sampled-data control,
dead time, state-space design and optimal control and estimation.

:doc:`Units, assumptions and numerical validation </guides/control_systems>`
gives the validity limits. Each script is self-contained and can be run
directly with ``python examples/control_systems/<section>/<script>.py``.

Sections
--------

- **feedback** -- Watt's governor, Minorsky's PID, Black's feedback amplifier and Ziegler-Nichols tuning.
- **stability** -- Maxwell's characteristic roots, Routh's array, Lyapunov's direct method and Hurwitz's determinants.
- **frequency** -- Nyquist's criterion and Bode's stability margins over every crossover.
- **time_response** -- Hazen's second-order servomechanism and Smith's predictor for dead time.
- **root_locus** -- Evans's root locus.
- **sampled** -- Ragazzini and Zadeh's zero-order-hold discretisation and Jury's stability test.
- **state_space** -- Kalman's controllability, Luenberger's observer and Wonham's pole placement.
- **optimal** -- Kalman's LQR and its discrete form, the Kalman-Bucy filter and Doyle's LQG counterexample.

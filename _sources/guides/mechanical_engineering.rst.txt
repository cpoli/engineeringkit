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

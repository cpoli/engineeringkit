Numerical methods and inspectable histories
===========================================

The shared integrator returns time [s], a (samples, states) history array,
method name and derivative evaluation count. Every sampled state is inspectable.
Euler and RK4 use the caller's time grid; adaptive methods use SciPy solve_ivp
with separate internal steps and interpolate onto that grid.

.. code-block:: python

    import numpy as np
    from engineeringkit.integrators import integrate_ode
    result = integrate_ode(lambda t, y: -y, [1.], np.linspace(0, 2, 101), "rk4")
    print(result.summary())
    fig, ax = result.plot(["normalized state"])

Euler has first-order global accuracy; classical RK4 has fourth-order global
accuracy for sufficiently smooth problems. Neither is unconditionally stable.
For decay y'=-a*y, Euler requires 0<a*dt<2 for stability. Prefer Radau or BDF
for stiff problems. Numerical tolerances control local error, not a proof of
physical fidelity. The tests check decay against exp(-t), RK4 convergence,
undamped oscillator energy, and critically damped response.

SI policy
---------

Plain numbers have no dimensional metadata. Callers must supply SI values and
use the units conversion helpers explicitly. Labels on plots should include
physical units. EngineeringKit does not automatically detect dimensional mistakes.

Regression behavior
-------------------

PID with ki=0 omits the integrator, preventing a spurious pole at zero.
Routh stability excludes imaginary-axis roots and handles auxiliary zero rows.
Root locus requires a strictly proper open-loop model and does not track branch
ordering between gains. ADC quantization uses 2**bits mid-rise reconstruction
levels, with saturation at the first and last bin centers. Fin efficiency at
zero m*L is its analytic limit, one.

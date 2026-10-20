r"""
Kalman's linear-quadratic regulator
===================================

Rudolf Kalman (1960) asked for the state feedback :math:`u = -Kx` that
minimises :math:`J = \int_0^\infty (x^\mathsf{T}Qx + u^\mathsf{T}Ru)\,dt`.
The answer is :math:`K = R^{-1}B^\mathsf{T}P`, where :math:`P` solves an
algebraic Riccati equation, and the minimum cost from :math:`x_0` is
:math:`x_0^\mathsf{T}Px_0` (:func:`~engineeringkit.control_systems.lqr`).
The weights replace pole locations as the designer's knobs: a cheaper
input (smaller :math:`R`) buys a faster response.

For the double integrator :math:`\ddot y = u` with :math:`Q =
\mathrm{diag}(1, 0)` and :math:`R = r`, the gain is
:math:`K = [r^{-1/2},\ \sqrt2\,r^{-1/4}]`: a closed-loop damping ratio of
:math:`1/\sqrt2` for every :math:`r`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from engineeringkit.control_systems import discretize_zoh, dlqr, lqr

A = np.array([[0.0, 1.0], [0.0, 0.0]])
B = np.array([[0.0], [1.0]])
Q = np.diag([1.0, 0.0])
x0 = np.array([1.0, 0.0])

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
for r in (0.01, 0.1, 1.0):
    K, P, closed = lqr(A, B, Q, [[r]])
    assert np.allclose(K, [[r**-0.5, np.sqrt(2) * r**-0.25]])
    zeta = -closed.real[0] / np.abs(closed[0])
    assert np.isclose(zeta, 1 / np.sqrt(2))

    # Simulate the closed loop, integrating the cost alongside the state.
    def rhs(t, z, K=K, r=r):
        x = z[:2]
        u = -(K @ x)
        return np.concatenate((A @ x + B @ u, [x @ Q @ x + r * u[0] ** 2]))

    sol = solve_ivp(rhs, (0, 40), np.concatenate((x0, [0.0])), rtol=1e-10, atol=1e-12, dense_output=True)
    cost = sol.y[2, -1]
    print(f"r = {r:5.2f}: K = {K.round(3).tolist()}, J = {cost:.5f}, x0'Px0 = {x0 @ P @ x0:.5f}")
    assert np.isclose(cost, x0 @ P @ x0, rtol=1e-6)

    t = np.linspace(0, 8, 400)
    x = sol.sol(t)[:2]
    ax1.plot(t, x[0], label=f"r = {r}")
    ax2.plot(t, -(K @ x)[0], label=f"r = {r}")

# %%
# The discrete regulator converges to the continuous one
# ------------------------------------------------------
#
# Bellman's dynamic programming over sampling periods gives the discrete
# Riccati equation (:func:`~engineeringkit.control_systems.dlqr`). With the
# weights scaled by the period, its gain tends to the continuous gain as the
# period shrinks.

K, _, _ = lqr(A, B, Q, [[0.1]])
for dt in (0.1, 0.01, 0.001):
    Ad, Bd = discretize_zoh(A, B, dt)
    Kd, _, _ = dlqr(Ad, Bd, Q * dt, [[0.1 * dt]])
    print(f"dt = {dt}: discrete gain {Kd.round(4).tolist()}")
assert np.allclose(Kd, K, rtol=2e-3)

ax1.set(ylabel="position y", title="LQR from y(0) = 1")
ax2.set(xlabel="time [s]", ylabel="input u")
ax1.legend()

plt.show()

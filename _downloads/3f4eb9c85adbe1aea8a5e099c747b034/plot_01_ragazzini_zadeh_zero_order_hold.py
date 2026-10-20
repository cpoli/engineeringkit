r"""
Ragazzini and Zadeh's sampled-data systems
==========================================

When a digital computer closes the loop it reads the plant only at sampling
instants and holds each command constant until the next. John Ragazzini and
Lotfi Zadeh (1952) developed the z-transform analysis of such loops. For a
zero-order hold the continuous model is replaced *exactly* by
:math:`x_{k+1} = e^{AT}x_k + \int_0^T e^{A\tau}d\tau\,B\,u_k`
(:func:`~engineeringkit.control_systems.discretize_zoh`). The naive Euler
step :math:`x_{k+1} = (I + AT)x_k + BTu_k` is only an approximation.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from engineeringkit.control_systems import discretize_zoh

A = np.array([[0.0, 1.0], [-4.0, -0.4]])  # mass-spring-damper
B = np.array([[0.0], [1.0]])
dt, n = 0.25, 60
u = np.sin(0.7 * dt * np.arange(n))  # held commands
Ad, Bd = discretize_zoh(A, B, dt)

# %%
# Exact recursion against a fine continuous simulation
# ----------------------------------------------------

exact, euler, reference = [np.zeros(2)], [np.zeros(2)], [np.zeros(2)]
for k in range(n):
    exact.append(Ad @ exact[-1] + Bd[:, 0] * u[k])
    euler.append((np.eye(2) + A * dt) @ euler[-1] + B[:, 0] * dt * u[k])
    sol = solve_ivp(lambda t, x: A @ x + B[:, 0] * u[k], (0, dt), reference[-1], rtol=1e-11, atol=1e-13)
    reference.append(sol.y[:, -1])
exact, euler, reference = map(np.array, (exact, euler, reference))
print(f"max error: ZOH {np.abs(exact - reference).max():.1e}, Euler {np.abs(euler - reference).max():.1e}")
assert np.abs(exact - reference).max() < 1e-8 and np.abs(euler - reference).max() > 0.05

tk = dt * np.arange(n + 1)
fig, ax = plt.subplots()
ax.plot(tk, reference[:, 0], "k-", lw=1, label="continuous plant")
ax.plot(tk, exact[:, 0], "o", ms=4, label="ZOH discretisation")
ax.plot(tk, euler[:, 0], "x", ms=4, label="Euler approximation")
ax.step(tk[:-1], 0.1 * u, where="post", color="gray", lw=0.8, label="held command / 10")
ax.set(xlabel="time [s]", ylabel="position", title=f"Sampling period {dt} s")
ax.legend()

plt.show()

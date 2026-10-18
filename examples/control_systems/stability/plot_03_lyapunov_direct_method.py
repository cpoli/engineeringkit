r"""
Lyapunov's direct method (1892)
===============================

Aleksandr Lyapunov's doctoral thesis (1892) proved stability without
solving the equations of motion. Find an energy-like function
:math:`V(x) > 0` that decreases along every trajectory, and the equilibrium
is stable. For linear dynamics :math:`\dot x = Ax`, the quadratic
:math:`V = x^\mathsf{T} P x` works exactly when the Lyapunov equation
:math:`A^\mathsf{T}P + PA = -Q` has a positive definite solution
(:func:`~engineeringkit.control_systems.lyapunov_equation`,
:func:`~engineeringkit.control_systems.lyapunov_stable`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

from engineeringkit.control_systems import lyapunov_equation, lyapunov_stable

A = np.array([[0.0, 1.0], [-4.0, -0.6]])  # lightly damped oscillator
P = lyapunov_equation(A)
print("P =\n", P.round(4))
assert lyapunov_stable(A) and np.all(np.linalg.eigvalsh(P) > 0)
assert not lyapunov_stable(np.array([[0.0, 1.0], [-4.0, 0.2]]))

# %%
# V decreases at the rate -x'x
# ----------------------------

t = np.linspace(0, 15, 3001)
x = np.array([expm(A * ti) @ [1.0, 0.0] for ti in t])
V = np.einsum("ti,ij,tj->t", x, P, x)
dV = np.gradient(V, t)
assert np.all(np.diff(V) < 0)
assert np.allclose(dV[1:-1], -np.sum(x**2, axis=1)[1:-1], atol=2e-3)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
g = np.linspace(-1.2, 1.2, 201)
X, Y = np.meshgrid(g, 2.4 * g)
levels = P[0, 0] * X**2 + 2 * P[0, 1] * X * Y + P[1, 1] * Y**2
axes[0].contour(X, Y, levels, levels=np.geomspace(0.02, V[0], 8), cmap="viridis")
axes[0].plot(x[:, 0], x[:, 1], "k", lw=1)
axes[0].set(xlabel="position", ylabel="velocity", title="Trajectory crosses every level set inward")
axes[1].semilogy(t, V)
axes[1].set(xlabel="time [s]", ylabel="$V = x^T P x$")
fig.tight_layout()

plt.show()

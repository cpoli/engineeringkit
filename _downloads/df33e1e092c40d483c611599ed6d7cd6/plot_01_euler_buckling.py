r"""
Euler's column buckling
=======================

In an appendix on elastic curves to his *Methodus inveniendi* (1744),
Leonhard Euler found the load at which a straight, pin-ended column can
stay in equilibrium bent as well as straight:
:math:`P_{cr} = \pi^2EI/L^2`. Other end conditions change the length of the
half sine wave, which is what the effective-length factor :math:`K`
expresses (:func:`~engineeringkit.structural_engineering.euler_buckling_load`).
Dividing by the area gives the Euler stress :math:`\pi^2E/(KL/r)^2`, a
function of the slenderness alone
(:func:`~engineeringkit.structural_engineering.slenderness_ratio`,
:func:`~engineeringkit.structural_engineering.radius_of_gyration`). The
example recovers Euler's load as the smallest eigenvalue of the
finite-difference column equation :math:`EIy'' + Py = 0`, then plots the
Euler hyperbola against the yield stress.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh_tridiagonal

from engineeringkit.structural_engineering import euler_buckling_load, radius_of_gyration, slenderness_ratio

# %%
# The buckling load as an eigenvalue
# ----------------------------------
#
# With :math:`y = 0` at both pins, the central difference
# :math:`y'' \approx (y_{i-1} - 2y_i + y_{i+1})/h^2` turns
# :math:`-EIy'' = Py` into a symmetric tridiagonal eigenproblem.

E, L = 200e9, 3.0
d = 0.05  # solid round bar [m]
A, I = np.pi * d**2 / 4, np.pi * d**4 / 64
n = 400
h = L / (n + 1)
loads, modes = eigh_tridiagonal(np.full(n, 2 * E * I / h**2), np.full(n - 1, -E * I / h**2), select="i", select_range=(0, 2))
P_cr = euler_buckling_load(E, I, L)
print(f"Euler P_cr = {P_cr / 1e3:.2f} kN, finite differences {loads[0] / 1e3:.2f} kN")
print("higher modes / first mode:", np.round(loads / loads[0], 3))
assert np.isclose(loads[0], P_cr, rtol=1e-4)
assert np.allclose(loads / loads[0], [1, 4, 9], rtol=1e-3)

x = np.linspace(0, L, n + 2)
fig, ax = plt.subplots(figsize=(6, 3))
for k in range(3):
    mode = np.concatenate([[0], modes[:, k], [0]])
    ax.plot(x, mode / np.max(abs(mode)), label=f"mode {k + 1}: {loads[k] / P_cr:.0f} P_cr")
ax.set(xlabel="x [m]", ylabel="deflection [-]", title="Buckled shapes of a pin-ended column")
ax.legend()
fig.tight_layout()

# %%
# The Euler hyperbola
# -------------------
#
# Stocky columns crush at the yield stress before they buckle; the two
# curves cross at the slenderness :math:`\pi\sqrt{E/\sigma_y}`.

sigma_y = 250e6
r = radius_of_gyration(I, A)
lengths = np.linspace(0.2, 6.0, 300)
fig, ax = plt.subplots(figsize=(6, 3.5))
for K, label in [(0.5, "fixed-fixed"), (1.0, "pinned-pinned"), (2.0, "fixed-free")]:
    lam = slenderness_ratio(lengths, K, r)
    stress = euler_buckling_load(E, I, lengths, K) / A
    ax.plot(lam, np.minimum(stress, sigma_y) / 1e6, label=f"K = {K} ({label})")
ax.axhline(sigma_y / 1e6, color="gray", ls=":")
transition = np.pi * np.sqrt(E / sigma_y)
ax.axvline(transition, color="gray", ls="--")
ax.set(xlabel="slenderness KL/r [-]", ylabel="failure stress [MPa]", ylim=(0, 300), xlim=(0, 400))
ax.legend()
fig.tight_layout()
print(f"transition slenderness pi*sqrt(E/sigma_y) = {transition:.0f}")
assert np.isclose(euler_buckling_load(E, I, transition * r) / A, sigma_y)

plt.show()

r"""
Clapeyron's three-moment equation (1857)
========================================

A beam continuous over several supports is statically indeterminate: the
support moments cannot be found from equilibrium alone. Émile Clapeyron
met the problem rebuilding the Asnières railway bridge over the Seine and in
1857 published one equation per interior support, linking each support
moment to its two neighbours
(:func:`~engineeringkit.structural_engineering.three_moment_equation`). The
example checks the textbook coefficients for equal spans, then builds the
moment diagram of an unequal three-span beam and checks by double
integration that the beam really does not move at the interior supports.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

from engineeringkit.structural_engineering import three_moment_equation

# %%
# Equal spans
# -----------
#
# Two equal spans give :math:`-wL^2/8` at the middle support, three give
# :math:`-wL^2/10` at both interior supports.

w, L = 10e3, 5.0
two = three_moment_equation([L, L], w)
three = three_moment_equation([L, L, L], w)
print("two spans  :", np.round(two / (w * L**2), 4), "x wL^2")
print("three spans:", np.round(three / (w * L**2), 4), "x wL^2")
assert np.isclose(two[1], -w * L**2 / 8)
assert np.allclose(three[1:3], -w * L**2 / 10)

# %%
# An unequal three-span beam
# --------------------------
#
# Between supports, the moment is the simply supported parabola plus a
# straight line joining the support moments.

spans, loads, E, I = np.array([6.0, 8.0, 5.0]), np.array([12e3, 12e3, 20e3]), 200e9, 4e-4
M_support = three_moment_equation(spans, loads)
edges = np.r_[0, np.cumsum(spans)]
x = np.linspace(0, edges[-1], 3801)
k = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, spans.size - 1)
s = x - edges[k]
M = loads[k] * s * (spans[k] - s) / 2 + M_support[k] * (1 - s / spans[k]) + M_support[k + 1] * s / spans[k]
print("support moments [kN m]:", np.round(M_support / 1e3, 2))

# %%
# Compatibility check: integrate :math:`y'' = M/(EI)` twice, then remove the
# straight line that puts the two end supports at zero. The deflection at
# the interior supports must also come out zero.

slope = cumulative_trapezoid(M / (E * I), x, initial=0)
y = cumulative_trapezoid(slope, x, initial=0)
y -= y[-1] * x / x[-1]
at_supports = np.interp(edges[1:-1], x, y)
print(f"deflection at interior supports: {np.round(1e3 * at_supports, 5)} mm (peak {1e3 * np.abs(y).max():.2f} mm)")
assert np.all(np.abs(at_supports) < 1e-3 * np.abs(y).max())

fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7, 5))
axes[0].plot(x, M / 1e3)
axes[0].plot(edges, M_support / 1e3, "ko")
axes[0].set(ylabel="moment [kN m]", title="Continuous beam over four supports")
axes[1].plot(x, 1e3 * y)
axes[1].plot(edges, np.zeros_like(edges), "k^")
axes[1].set(xlabel="x [m]", ylabel="deflection [mm]")
for ax in axes:
    ax.axhline(0, color="gray", lw=0.5)
fig.tight_layout()

plt.show()

r"""
Routh's stability criterion (1877)
==================================

Maxwell could not state a test for polynomials above fourth order. Edward
Routh's Adams Prize essay (1877) solved the problem with a table built from
the coefficients alone: the number of sign changes down its first column
(:func:`~engineeringkit.control_systems.routh_first_column`) equals the
number of roots in the right half-plane. For the closed loop of
:math:`K/(s(s+1)(s+2))`, :math:`s^3 + 3s^2 + 2s + K`, the table shows the
loop is stable only for :math:`0 < K < 6`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.control_systems import routh_first_column, routh_stable

# %%
# Sign changes count right-half-plane roots
# -----------------------------------------

gains = np.linspace(0.05, 15.05, 151)  # avoids the marginal K = 6
changes, unstable_roots = [], []
for K in gains:
    column = routh_first_column([1, 3, 2, K])
    changes.append(int(np.sum(np.diff(np.sign(column)) != 0)))
    unstable_roots.append(int(np.sum(np.roots([1, 3, 2, K]).real > 0)))
assert changes == unstable_roots
assert all(routh_stable([1, 3, 2, K]) == (K < 6) for K in gains)

for K in [2.0, 10.0]:
    print(f"K = {K:4.1f}: first column {np.round(routh_first_column([1, 3, 2, K]), 3)}")

# %%
# A fifth-order example
# ---------------------

p = np.poly([-1, -2, 0.5 + 1j, 0.5 - 1j, -3]).real
column = routh_first_column(p)
print("fifth order:", np.round(column, 3), "->", int(np.sum(np.diff(np.sign(column)) != 0)), "sign changes")
assert int(np.sum(np.diff(np.sign(column)) != 0)) == 2

fig, ax = plt.subplots()
ax.step(gains, changes, where="mid", label="Routh sign changes")
ax.plot(gains, unstable_roots, "o", ms=3, mfc="none", label="right-half-plane roots")
ax.axvline(6, color="C3", ls="--", lw=0.8)
ax.set(xlabel="loop gain $K$", ylabel="count", title="$s^3 + 3s^2 + 2s + K$")
ax.legend()

plt.show()

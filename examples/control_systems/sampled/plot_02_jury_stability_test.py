r"""
Jury's stability test for sampled systems (1961)
================================================

A sampled loop is stable when every root of its characteristic polynomial in
:math:`z` lies inside the unit circle, not in the left half-plane, so Routh's
test does not apply directly. Eliahu Jury (1961) gave the discrete
counterpart: a table that reduces the polynomial one degree at a time and
checks :math:`|a_n| < |a_0|` at each step, with no roots computed
(:func:`~engineeringkit.control_systems.jury_stable`).

The plant :math:`1/(s(s+1))` under proportional control :math:`K` is stable
for every :math:`K > 0` in continuous time. Sampled with a zero-order hold of
period :math:`T`, it becomes
:math:`G(z) = (b_1 z + b_2)/((z - 1)(z - a))` with :math:`a = e^{-T}`,
:math:`b_1 = T - 1 + a` and :math:`b_2 = 1 - a - Ta`. The closed loop
:math:`z^2 + (Kb_1 - 1 - a)z + (a + Kb_2)` loses stability when its constant
term reaches one, at :math:`K = (1 - a)/b_2`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.optimize import brentq

from engineeringkit.control_systems import discretize_zoh, jury_stable, routh_stable

A, B, C = np.array([[0.0, 1.0], [0.0, -1.0]]), np.array([[0.0], [1.0]]), np.array([[1.0, 0.0]])


def characteristic(K, T):
    Ad, Bd = discretize_zoh(A, B, T)
    num, den = signal.ss2tf(Ad, Bd, C, [[0.0]])
    return np.polyadd(den, K * num[0])


# %%
# Jury against the roots
# ----------------------

T = 0.5
gains = np.linspace(0.1, 8, 400)
jury = np.array([jury_stable(characteristic(K, T)) for K in gains])
radius = np.array([np.abs(np.roots(characteristic(K, T))).max() for K in gains])
assert np.array_equal(jury, radius < 1)
assert all(routh_stable([1, 1, K]) for K in gains)  # the continuous loop never goes unstable

a = np.exp(-T)
k_exact = (1 - a) / (1 - a - T * a)
k_jury = brentq(lambda K: 0.5 - jury_stable(characteristic(K, T)), 1, 8, xtol=1e-12)
print(f"T = {T} s: critical gain {k_jury:.6f} by Jury, {k_exact:.6f} analytically")
assert np.isclose(k_jury, k_exact, rtol=1e-9)

# %%
# Slower sampling, smaller gain
# -----------------------------
#
# The hold delays the command by about half a period, so the critical gain
# falls as :math:`T` grows, and grows without bound as :math:`T \to 0`.

periods = np.logspace(-1.5, 0.5, 40)
k_crit = [brentq(lambda K, T=T: 0.5 - jury_stable(characteristic(K, T)), 1e-3, 1e3, xtol=1e-10) for T in periods]
a = np.exp(-periods)
assert np.allclose(k_crit, (1 - a) / (1 - a - periods * a), rtol=1e-8)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
ax1.plot(gains, radius, label="largest root modulus")
ax1.fill_between(gains, 0, 1.6, where=jury, color="C2", alpha=0.15, label="Jury: stable")
ax1.axhline(1, color="k", lw=0.8)
ax1.axvline(k_exact, color="C3", ls="--", label=f"$K = (1-a)/b_2$ = {k_exact:.2f}")
ax1.set(xlabel="gain K", ylabel="|z|", ylim=(0.6, 1.6), title=f"Sampling period T = {T} s")
ax1.legend()
ax2.loglog(periods, k_crit, "o", ms=4, label="Jury bisection")
ax2.loglog(periods, (1 - a) / (1 - a - periods * a), "k-", lw=1, label="analytic")
ax2.set(xlabel="sampling period T [s]", ylabel="critical gain", title="Continuous loop: stable for all K")
ax2.legend()
fig.tight_layout()

plt.show()

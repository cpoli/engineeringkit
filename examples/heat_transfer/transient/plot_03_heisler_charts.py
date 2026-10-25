r"""
Heisler's transient conduction charts (1947)
============================================

Separation of variables gives the cooling of a plane wall as a Fourier
series whose eigenvalues solve :math:`\zeta\tan\zeta = Bi`. Before
computers, nobody summed it by hand. In 1947 M. P. Heisler, working on
induction heating, noticed that after a short time only the first term
matters, and plotted that term as charts of centre temperature against the
Fourier number for each Biot number. Generations of students read them with
a ruler. :func:`~engineeringkit.heat_transfer.plane_wall_transient` sums the
series (or its first term, ``terms=1``); this example redraws the centre
chart, measures the one-term error, and times the quench of a steel plate.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from engineeringkit.heat_transfer import biot_number, fourier_number, plane_wall_transient

# %%
# The centre-temperature chart
# ----------------------------
#
# Heisler's chart is labelled by :math:`1/Bi`; the curves are straight on a
# log scale once the first term dominates.

Fo = np.linspace(0, 4, 201)
fig, ax = plt.subplots(figsize=(6, 4))
for inv_bi in [0.0001, 0.5, 1, 2, 5, 10]:
    theta = plane_wall_transient(Fo, 1 / inv_bi, 0.0, terms=200)
    ax.semilogy(Fo, theta, label=f"1/Bi = {inv_bi:g}")
ax.set(xlabel="Fourier number", ylabel="centre excess temperature", ylim=(1e-3, 1.1), title="Heisler's plane-wall centre chart")
ax.legend(ncol=2)
fig.tight_layout()

# %%
# How good is one term?
# ---------------------
#
# At the centre the one-term approximation stays within 2% of the full
# series for :math:`Fo > 0.2`, whatever the Biot number. Very early, many
# terms are needed: at :math:`Fo = 10^{-3}` the centre has not yet felt
# the surface.

worst = 0.0
Fo_late = np.linspace(0.2, 3, 60)
for Bi in [0.01, 0.1, 1, 10, 100]:
    full = plane_wall_transient(Fo_late, Bi, 0.0)
    worst = max(worst, np.max(np.abs(plane_wall_transient(Fo_late, Bi, 0.0, terms=1) / full - 1)))
print(f"largest one-term error at the centre for Fo > 0.2: {worst:.2%}")
assert worst < 0.02
early = plane_wall_transient(1e-3, 10.0, 0.0, terms=200)
print(f"centre at Fo = 0.001 (200 terms): {early:.6f}")
assert np.isclose(early, 1.0, atol=1e-6)

# %%
# Quenching a steel plate
# -----------------------
#
# A 50 mm steel plate (:math:`k = 40` W/(m K), :math:`\alpha = 1.1\times10^{-5}`
# m²/s) at 1100 K is quenched in oil at 330 K with :math:`h = 1000`
# W/(m² K). When does its centre reach 600 K?

k, alpha, L, h, Ti, Tinf, target = 40.0, 1.1e-5, 0.025, 1000.0, 1100.0, 330.0, 600.0
Bi = biot_number(h, L, k)
theta_target = (target - Tinf) / (Ti - Tinf)
Fo_target = brentq(lambda F: plane_wall_transient(F, Bi, 0.0) - theta_target, 0.01, 50)
t_target = Fo_target * L**2 / alpha
print(f"Bi = {Bi:.3f}, Fo = {Fo_target:.3f}: the centre reaches 600 K after {t_target:.0f} s")
assert np.isclose(fourier_number(alpha, t_target, L), Fo_target)
surface = Tinf + (Ti - Tinf) * plane_wall_transient(Fo_target, Bi, 1.0)
print(f"the surface is then at {surface:.0f} K")

plt.show()

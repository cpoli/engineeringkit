r"""
Coulomb's friction and the power screw
======================================

Charles-Augustin de Coulomb's prize memoir *Théorie des machines simples*
(1781, published 1785) measured sliding friction and found it proportional
to the normal load and nearly independent of contact area and speed:
:math:`F = \mu N`. A screw thread is an inclined plane wrapped round a
cylinder. Unwrapping one turn gives a slope of lead angle
:math:`\tan\lambda = L/(\pi d_m)`, and Coulomb's law on that slope gives the
torque to raise a load,
:math:`T = \tfrac12 W d_m\tan(\phi + \lambda)` with :math:`\tan\phi = \mu`
(:func:`~engineeringkit.mechanical_engineering.screw_torque_raise`). The
example solves the force balance on the unwrapped thread directly and
shows why a self-locking screw is always less than 50% efficient.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import screw_torque_raise

# %%
# The unwrapped thread
# --------------------
#
# A block carrying the axial load :math:`W` is pushed up a slope
# :math:`\lambda` by a horizontal force :math:`P`, the thread force at the
# mean radius. Balance along and normal to the slope, with friction
# :math:`\mu N` opposing the motion, is linear in :math:`P` and :math:`N`.

W, dm, lead = 20e3, 0.032, 0.006
lam = np.arctan(lead / (np.pi * dm))
for mu in (0.0, 0.08, 0.15):
    A = np.array([[np.cos(lam), -mu], [-np.sin(lam), 1.0]])
    P, N = np.linalg.solve(A, [W * np.sin(lam), W * np.cos(lam)])
    T = screw_torque_raise(W, dm, lead, mu)
    print(f"mu = {mu:.2f}: force balance T = {P * dm / 2:7.3f} N m, formula T = {T:7.3f} N m")
    assert np.isclose(P * dm / 2, T)
assert np.isclose(2 * np.pi * screw_torque_raise(W, dm, lead, 0.0), W * lead)

# %%
# Efficiency and self-locking
# ---------------------------
#
# Efficiency is the lifting work per turn over the input work,
# :math:`\eta = W L / (2\pi T) = \tan\lambda / \tan(\phi + \lambda)`. A screw
# holds its load without a brake when :math:`\phi \ge \lambda`; at that
# limit :math:`\eta = (1 - \tan^2\lambda)/2`, below one half.

leads = np.linspace(0.5e-3, 40e-3, 200)
fig, ax = plt.subplots(figsize=(6, 3.8))
for mu in (0.08, 0.15, 0.25):
    eta = W * leads / (2 * np.pi * screw_torque_raise(W, dm, leads, mu))
    locking = np.tan(np.arctan(leads / (np.pi * dm))) <= mu
    ax.plot(np.degrees(np.arctan(leads / (np.pi * dm))), eta, label=f"mu = {mu}")
    assert np.all(eta[locking] < 0.5)
ax.axhline(0.5, color="k", ls=":")
ax.set(xlabel="lead angle [deg]", ylabel="efficiency", title="Square-thread power screw")
ax.legend()
fig.tight_layout()

limit = W * lead / (2 * np.pi * screw_torque_raise(W, dm, lead, np.tan(lam)))
print(f"lead angle {np.degrees(lam):.2f} deg; efficiency at the self-locking limit {limit:.4f}")
assert np.isclose(limit, (1 - np.tan(lam) ** 2) / 2)

plt.show()

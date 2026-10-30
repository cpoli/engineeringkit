r"""
Graetz's laminar tube flow (1883)
=================================

In 1883 Leo Graetz, measuring the conductivity of liquids flowing through
heated tubes, solved for the temperature of a fluid in laminar Poiseuille
flow entering a tube whose wall is held at a new temperature. The answer is
a series of eigenfunctions; far from the inlet only the first survives and
the Nusselt number settles at :math:`Nu_D = 3.657`
(:func:`~engineeringkit.heat_transfer.laminar_tube_nusselt`). A uniform
wall flux gives :math:`48/11` instead. Near the inlet the thin thermal
boundary layer gives higher values, which Hausen fitted as a function of the
Graetz number (:func:`~engineeringkit.heat_transfer.hausen_nusselt`), and the
flow is fully developed after the entry length
(:func:`~engineeringkit.heat_transfer.entry_length`). This example marches
Graetz's problem numerically down the tube and recovers all three.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid
from scipy.linalg import solve_banded

from engineeringkit.heat_transfer import entry_length, hausen_nusselt, laminar_tube_nusselt

# %%
# Marching the Graetz problem
# ---------------------------
#
# In radius :math:`r^* = r/R` and axial distance :math:`x^* = (x/D)/(Re\,Pr)`,
# the energy equation is :math:`(1 - r^{*2})\,\partial\theta/\partial x^* =
# (2/r^*)\,\partial(r^*\partial\theta/\partial r^*)/\partial r^*`, with
# :math:`\theta = 1` at the inlet and :math:`\theta = 0` at the wall. An
# implicit (backward Euler) step in :math:`x^*` keeps it stable on a fine
# radial grid.

n = 400
r = np.linspace(0, 1, n + 1)
dr = r[1]
u = 1 - r**2  # velocity over the centreline velocity
theta = np.ones(n + 1)
theta[-1] = 0.0


def step_matrix(dx):
    """Banded form of (u/dx) I - (2/r) d/dr (r d/dr) on the interior nodes and the axis."""
    ab = np.zeros((3, n))
    for i in range(n):
        if i == 0:  # axis: (1/r) d/dr (r dT/dr) -> 2 d2T/dr2 by symmetry
            diag, right, left = 8 / dr**2, -8 / dr**2, 0.0
        else:
            rp, rm = r[i] + dr / 2, r[i] - dr / 2
            diag, right, left = 2 * (rp + rm) / (r[i] * dr**2), -2 * rp / (r[i] * dr**2), -2 * rm / (r[i] * dr**2)
        ab[1, i] = u[i] / dx + diag
        if i + 1 < n:
            ab[0, i + 1] = right
        if i > 0:
            ab[2, i - 1] = left
    return ab


x_star = [0.0]
nusselt = []
dx = 1e-5
while x_star[-1] < 0.3:
    ab = step_matrix(dx)
    theta[:-1] = solve_banded((1, 1), ab, u[:-1] * theta[:-1] / dx)
    x_star.append(x_star[-1] + dx)
    # Bulk temperature (velocity-weighted) and wall gradient give the local Nusselt number.
    bulk = trapezoid(u * theta * r, r) / trapezoid(u * r, r)
    gradient = (3 * theta[-1] - 4 * theta[-2] + theta[-3]) / (2 * dr)
    nusselt.append(-2 * gradient / bulk)
    dx = min(dx * 1.05, 2e-3)
x_star, nusselt = np.array(x_star[1:]), np.array(nusselt)
print(f"local Nu far downstream: {nusselt[-1]:.4f} (Graetz limit {laminar_tube_nusselt('temperature'):.4f})")
assert abs(nusselt[-1] / laminar_tube_nusselt("temperature") - 1) < 2e-3

# %%
# The entry region and Hausen's average
# -------------------------------------
#
# The mean Nusselt number over a length is set by the log of the bulk
# temperature drop, :math:`\bar{Nu} = -\ln\theta_b/(4x^*)`. Hausen's fit uses
# the Graetz number :math:`Gz = 1/x^*`. The thermal entry length
# :math:`0.05\,Re\,Pr\,D` is :math:`x^* = 0.05`, where the local value is within
# a few per cent of the limit.

bulk_history = np.exp(-4 * np.cumsum(nusselt * np.diff(np.concatenate([[0.0], x_star]))))
mean_nu = -np.log(bulk_history) / (4 * x_star)
Re, Pr, D = 1000.0, 5.0, 0.01  # any Re Pr D: only x* = (L/D)/(Re Pr) matters
L = x_star * D * Re * Pr
fit = hausen_nusselt(Re, Pr, D, L)
window = (x_star > 0.005) & (x_star < 0.3)
worst = np.max(np.abs(fit[window] / mean_nu[window] - 1))
print(f"Hausen's fit against the marched mean, 0.005 < x* < 0.3: within {worst:.1%}")
assert worst < 0.06
x_entry = entry_length(Re, D, Pr) / (D * Re * Pr)
local_at_entry = np.interp(x_entry, x_star, nusselt)
print(f"local Nu at the thermal entry length (x* = {x_entry:.2f}): {local_at_entry:.3f}")
assert abs(local_at_entry / 3.657 - 1) < 0.03

fig, ax = plt.subplots(figsize=(6, 4))
ax.loglog(x_star, nusselt, label="local, marched")
ax.loglog(x_star, mean_nu, label="mean, marched")
ax.loglog(x_star, fit, "--", label="Hausen")
ax.axhline(3.657, color="gray", lw=0.8, label="Graetz limit 3.66")
ax.axvline(x_entry, color="gray", ls=":", lw=0.8)
ax.set(xlabel="x* = (x/D)/(Re Pr)", ylabel="Nusselt number", title="Laminar tube, uniform wall temperature")
ax.legend()
fig.tight_layout()

# %%
# Uniform heat flux
# -----------------
#
# With a uniform flux the profile shape stops changing and the energy
# equation integrates in closed form: :math:`Nu_D = 48/11`.

print(f"uniform flux: Nu = {laminar_tube_nusselt('flux'):.4f} = 48/11")
assert np.isclose(laminar_tube_nusselt("flux"), 48 / 11)

plt.show()

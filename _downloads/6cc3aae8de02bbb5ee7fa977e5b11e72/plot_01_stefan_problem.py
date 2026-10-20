r"""
Stefan's problem: a moving freezing front
=========================================

In 1891 Josef Stefan asked how fast sea ice thickens in the polar winter.
The ice conducts heat from the freezing water below to the cold air above,
and every millimetre of new ice releases latent heat at the moving
interface, so the boundary is part of the unknown. Franz Neumann had found
the exact similarity solution in his lectures of the 1860s: the front
advances as :math:`s = 2\lambda\sqrt{\alpha t}` with
:math:`\lambda e^{\lambda^2}\operatorname{erf}\lambda = Ste/\sqrt\pi`
(:func:`~engineeringkit.heat_transfer.neumann_stefan_front`,
:func:`~engineeringkit.heat_transfer.stefan_number`). This example grows a
winter of sea ice, checks the energy balance at the front, and compares with
Stefan's own quasi-steady estimate.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import special

from engineeringkit.heat_transfer import neumann_stefan_front, stefan_number

# %%
# Sea ice under a -20 °C sky
# --------------------------
#
# Ice: :math:`k = 2.22` W/(m K), :math:`\rho = 917` kg/m³,
# :math:`c_p = 2050` J/(kg K), latent heat 334 kJ/kg. The top surface is held
# 20 K below the freezing point.

k, rho, cp, latent, dT = 2.22, 917.0, 2050.0, 334e3, 20.0
alpha = k / (rho * cp)
Ste = stefan_number(cp, dT, latent)
day = 86400.0
days = np.array([1, 10, 30, 100])
s = neumann_stefan_front(days * day, alpha, Ste)
lam = s[0] / (2 * np.sqrt(alpha * day))
print(f"Ste = {Ste:.4f}, lambda = {lam:.4f}")
for d, thickness in zip(days, s):
    print(f"after {d:3d} days: {thickness * 100:5.1f} cm of ice")
assert np.isclose(lam * np.exp(lam**2) * special.erf(lam), Ste / np.sqrt(np.pi))

# %%
# Latent heat released equals heat conducted
# ------------------------------------------
#
# At the front, :math:`k\,\partial T/\partial x = \rho h_{sf}\,\mathrm{d}s/\mathrm{d}t`.
# Inside the ice the temperature is
# :math:`T_s + \Delta T\operatorname{erf}(x/2\sqrt{\alpha t})/\operatorname{erf}\lambda`.

t, dt = 30 * day, 1.0
front = neumann_stefan_front(t, alpha, Ste)


def ice_temperature(x, t):
    return -dT + dT * special.erf(x / (2 * np.sqrt(alpha * t))) / special.erf(lam)  # °C, 0 at the front


h = 1e-6
conducted = k * (ice_temperature(front, t) - ice_temperature(front - h, t)) / h
released = rho * latent * (neumann_stefan_front(t + dt, alpha, Ste) - neumann_stefan_front(t - dt, alpha, Ste)) / (2 * dt)
print(f"at 30 days: conducted {conducted:.3f} W/m², latent heat released {released:.3f} W/m²")
assert np.isclose(ice_temperature(front, t), 0.0, atol=1e-12)
assert np.isclose(conducted, released, rtol=1e-4)

# %%
# Stefan's quasi-steady estimate
# ------------------------------
#
# Neglecting the heat stored in the ice (a linear profile) gives
# :math:`s = \sqrt{2k\Delta T t/(\rho h_{sf})} = \sqrt{2\,Ste\,\alpha t}`,
# Stefan's law of ice growth. It overestimates by a few percent at this
# Stefan number and becomes exact as :math:`Ste \to 0`.

quasi = np.sqrt(2 * Ste * alpha * days * day)
print(f"quasi-steady overestimate: {(quasi / s - 1)[0]:.2%}")
for small in [1e-2, 1e-4]:
    ratio = np.sqrt(2 * small * alpha * day) / neumann_stefan_front(day, alpha, small)
    assert abs(ratio - 1) < small
assert 0 < quasi[0] / s[0] - 1 < 0.05

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
winter = np.linspace(0, 180, 200) * day
axes[0].plot(winter / day, neumann_stefan_front(winter, alpha, Ste) * 100, label="Neumann (exact)")
axes[0].plot(winter / day, np.sqrt(2 * Ste * alpha * winter) * 100, "--", label="Stefan quasi-steady")
axes[0].set(xlabel="days", ylabel="ice thickness [cm]", title="Growth as the square root of time")
axes[0].legend()
for d in [10, 30, 100]:
    x = np.linspace(0, neumann_stefan_front(d * day, alpha, Ste), 50)
    axes[1].plot(ice_temperature(x, d * day), x * 100, label=f"{d} days")
axes[1].invert_yaxis()
axes[1].set(xlabel="ice temperature [°C]", ylabel="depth below the surface [cm]", title="Profiles in the ice")
axes[1].legend()
fig.tight_layout()

plt.show()

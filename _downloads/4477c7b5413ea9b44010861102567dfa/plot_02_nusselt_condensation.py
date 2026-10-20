r"""
Nusselt's film condensation
===========================

Steam condensing on the cold tubes of a power-station condenser forms a
liquid film that drains under gravity, and the film itself is the main
resistance to heat flow. In 1916 Wilhelm Nusselt balanced gravity against
viscous shear in the film and conduction across it. The film thickens as
:math:`x^{1/4}` down the wall, and the average coefficient over a height
:math:`L` is
:math:`0.943\,[g\rho_l(\rho_l - \rho_v)k_l^3h_{fg}/(\mu_l L\Delta T)]^{1/4}`
(:func:`~engineeringkit.heat_transfer.nusselt_film_condensation`). This
example rebuilds the film profile, checks the average against the local
coefficient, and checks the energy balance.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

from engineeringkit.heat_transfer import nusselt_film_condensation

# %%
# Steam at 1 atm on a 90 °C wall
# ------------------------------
#
# Liquid properties at the 95 °C film temperature; a wall 0.3 m high.

rho_l, rho_v, mu_l, k_l, h_fg, g = 962.0, 0.60, 2.97e-4, 0.677, 2.257e6, 9.80665
Tsat, Ts, L = 373.15, 363.15, 0.3
h_avg = nusselt_film_condensation(rho_l, rho_v, mu_l, k_l, h_fg, Tsat, Ts, L)
print(f"average coefficient {h_avg:.0f} W/(m² K), heat flux {h_avg * (Tsat - Ts) / 1e3:.1f} kW/m²")

# %%
# The film profile
# ----------------
#
# :math:`\delta(x) = [4\mu_l k_l\Delta T x/(g\rho_l(\rho_l - \rho_v)h_{fg})]^{1/4}`,
# and the local coefficient is :math:`k_l/\delta`. Its average over the
# height is 4/3 of its value at the bottom.

x = np.linspace(0, L, 200_001)[1:]
delta = (4 * mu_l * k_l * (Tsat - Ts) * x / (g * rho_l * (rho_l - rho_v) * h_fg)) ** 0.25
h_local = k_l / delta
h_integrated = integrate.trapezoid(h_local, x) / L
print(f"film at the bottom {delta[-1] * 1e3:.3f} mm; averaged local h {h_integrated:.0f}, 4/3 of bottom {4 / 3 * h_local[-1]:.0f}")
assert np.isclose(h_integrated, h_avg, rtol=1e-3)
assert np.isclose(4 / 3 * h_local[-1], h_avg, rtol=1e-3)

# %%
# Energy and mass balance
# -----------------------
#
# The condensate flow leaving the bottom, per metre of width, is
# :math:`\rho_l(\rho_l - \rho_v)g\delta^3/(3\mu_l)`. Its latent heat must
# equal the heat conducted into the wall.

mass_flow = rho_l * (rho_l - rho_v) * g * delta[-1] ** 3 / (3 * mu_l)
heat = h_avg * L * (Tsat - Ts)
print(f"latent heat carried {mass_flow * h_fg:.0f} W/m, heat conducted {heat:.0f} W/m")
assert np.isclose(mass_flow * h_fg, heat, rtol=1e-3)

Re_film = 4 * mass_flow / mu_l
print(f"film Reynolds number {Re_film:.0f}: laminar (below 1800), with ripples above about 30 adding 10-20%")
assert Re_film < 1800

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
axes[0].plot(delta * 1e3, x * 100)
axes[0].invert_yaxis()
axes[0].set(xlabel="film thickness [mm]", ylabel="distance down the wall [cm]", title="Film grows as x^(1/4)")
dT = np.linspace(1, 40, 100)
axes[1].plot(dT, nusselt_film_condensation(rho_l, rho_v, mu_l, k_l, h_fg, Tsat, Tsat - dT, L))
axes[1].set(xlabel="Tsat - Twall [K]", ylabel="average h [W/(m² K)]", title="Thicker film, lower coefficient")
fig.tight_layout()

plt.show()

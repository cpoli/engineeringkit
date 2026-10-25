r"""
Nukiyama's boiling curve and Zuber's critical heat flux (1934)
==============================================================

In 1934 Shiro Nukiyama boiled water at atmospheric pressure on electrically
heated wires. Raising the power, he found that the heat flux reached a
maximum a few tens of kelvin above saturation. One step more and the wire
was suddenly covered by vapour: a nichrome wire melted, a platinum wire
glowed red. Controlling the temperature instead of the power later traced
the whole curve, through a falling transition region to a minimum and on
into film boiling. In 1958 Novak Zuber explained the maximum as a
hydrodynamic (Taylor and Helmholtz) instability of the vapour columns,
:math:`q''_{max} = C h_{fg}\rho_v[\sigma g(\rho_l - \rho_v)/\rho_v^2]^{1/4}`
(:func:`~engineeringkit.heat_transfer.zuber_critical_heat_flux`). This
example sketches the curve with standard correlations and finds Nukiyama's
jump.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from engineeringkit.heat_transfer import STEFAN_BOLTZMANN, zuber_critical_heat_flux

# %%
# Saturated water at 1 atm
# ------------------------

rho_l, rho_v, h_fg, sigma, g = 957.9, 0.5955, 2.257e6, 0.0589, 9.80665
mu_l, cp_l, Pr_l = 2.79e-4, 4217.0, 1.76
q_max = zuber_critical_heat_flux(rho_l, rho_v, h_fg, sigma)
q_max_lienhard = zuber_critical_heat_flux(rho_l, rho_v, h_fg, sigma, C=0.149)
print(f"critical heat flux: Zuber {q_max / 1e6:.2f} MW/m², Lienhard-Dhir {q_max_lienhard / 1e6:.2f} MW/m²")
assert 1.0e6 < q_max < 1.2e6

# %%
# Nucleate boiling (Rohsenow, 1952)
# ---------------------------------
#
# :math:`q'' = \mu_l h_{fg}[g(\rho_l - \rho_v)/\sigma]^{1/2}[c_p\Delta T/(C_{sf}h_{fg}Pr_l)]^3`
# with :math:`C_{sf} = 0.013` for water on clean metal. The cube makes the
# flux climb steeply until it meets Zuber's limit.


def nucleate(dT, C_sf=0.013):
    return mu_l * h_fg * np.sqrt(g * (rho_l - rho_v) / sigma) * (cp_l * dT / (C_sf * h_fg * Pr_l)) ** 3


dT_max = brentq(lambda dT: nucleate(dT) - q_max, 1, 100)
q_min = 0.09 * rho_v * h_fg * (sigma * g * (rho_l - rho_v) / (rho_l + rho_v) ** 2) ** 0.25
print(f"nucleate boiling reaches the maximum at {dT_max:.1f} K excess; Zuber's minimum (Leidenfrost) flux {q_min / 1e3:.1f} kW/m²")
assert 10 < dT_max < 40

# %%
# Film boiling on a 1 mm wire (Bromley, 1950)
# -------------------------------------------
#
# Conduction through the vapour film plus radiation, with rough steam
# properties at the film temperature. With the power held at
# :math:`q''_{max}`, the wire must jump along the dashed line to the film
# branch, more than a thousand kelvin hotter.

D = 1e-3


def film(dT):
    T_film = 373.15 + dT / 2
    rho = 101325 / (461.5 * T_film)
    k_v = 0.0249 + (0.0797 - 0.0249) * (T_film - 373.0) / 627.0
    mu_v = 1.2e-5 + (3.7e-5 - 1.2e-5) * (T_film - 373.0) / 627.0
    h_fg_eff = h_fg + 0.8 * 2000.0 * dT
    h_conv = 0.62 * (g * (rho_l - rho) * rho * k_v**3 * h_fg_eff / (mu_v * D * dT)) ** 0.25
    Tw = 373.15 + dT
    h_rad = 0.9 * STEFAN_BOLTZMANN * (Tw**4 - 373.15**4) / dT
    return (h_conv + 0.75 * h_rad) * dT


dT_jump = brentq(lambda dT: film(dT) - q_max, 50, 3000)
print(f"at constant power the wire jumps from {dT_max:.0f} K to about {dT_jump:.0f} K excess ({373 + dT_jump:.0f} K)")
assert dT_jump > 1000

fig, ax = plt.subplots(figsize=(6.5, 4.2))
dT_nb = np.geomspace(5, dT_max, 50)
dT_fb = np.geomspace(120, 2000, 50)
ax.loglog(dT_nb, nucleate(dT_nb), label="nucleate (Rohsenow)")
ax.loglog(dT_fb, film(dT_fb), label="film (Bromley + radiation)")
ax.loglog([dT_max, 120], [q_max, film(120)], ":", color="gray", label="transition (sketched)")
ax.loglog([dT_max, dT_jump], [q_max, q_max], "k--", label="Nukiyama's jump")
ax.axhline(q_min, color="C3", lw=0.8)
ax.plot(dT_max, q_max, "ko")
ax.set(xlabel="wall excess temperature [K]", ylabel="heat flux [W/m²]", title="Pool boiling of water at 1 atm")
ax.legend(fontsize=8)
fig.tight_layout()

plt.show()

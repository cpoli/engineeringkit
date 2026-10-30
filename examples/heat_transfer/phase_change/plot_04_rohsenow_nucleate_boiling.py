r"""
Rohsenow's nucleate boiling correlation (1952)
==============================================

In nucleate boiling, bubbles growing and leaving the heater pump liquid to
the surface, and the heat flux rises steeply with the wall superheat. At
MIT, Warren Rohsenow treated each departing bubble as the agent of a
convection process and correlated the data for many liquids with one
equation whose only empirical input is a surface-fluid constant
:math:`C_{s,f}`:
:math:`q'' = \mu_l h_{fg}[g(\rho_l - \rho_v)/\sigma]^{1/2}
[c_{p,l}\Delta T_e/(C_{s,f}h_{fg}Pr_l^n)]^3`
(:func:`~engineeringkit.heat_transfer.rohsenow_boiling_heat_flux`). This
example boils water in a pan, shows the cube law and how strongly the
surface matters, and finds where the curve meets Zuber's critical flux.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from engineeringkit.heat_transfer import rohsenow_boiling_heat_flux, zuber_critical_heat_flux

# %%
# Water at 1 atm in a copper pan
# ------------------------------
#
# Saturated water properties; the pan bottom is 18 K above the boiling point
# (Incropera et al., Example 10.1, which reports about 0.8 MW/m²).

water = {"mu_l": 279e-6, "h_fg": 2.257e6, "rho_l": 957.9, "rho_v": 0.5955, "sigma": 58.9e-3, "cp_l": 4217.0, "Pr_l": 1.76}
q = rohsenow_boiling_heat_flux(**water, dT_excess=18.0, C_sf=0.013)
area = np.pi * 0.30**2 / 4
print(f"q'' = {q / 1e3:.0f} kW/m²: a 30 cm pan boils off {q * area / water['h_fg'] * 3600:.0f} kg of water per hour")
assert 0.75e6 < q < 0.85e6

# %%
# The cube law and the surface constant
# -------------------------------------
#
# Doubling the superheat multiplies the flux by eight. Turned around, an
# uncertainty in :math:`C_{s,f}` changes the flux by its cube but the
# superheat only in proportion: the correlation predicts wall temperature far
# better than heat flux.

dT = np.linspace(2, 30, 200)
ratio = rohsenow_boiling_heat_flux(**water, dT_excess=20.0) / rohsenow_boiling_heat_flux(**water, dT_excess=10.0)
print(f"q''(20 K) / q''(10 K) = {ratio:.3f}")
assert np.isclose(ratio, 8.0)

fig, ax = plt.subplots(figsize=(6, 4))
for C_sf, surface in [(0.0060, "brass"), (0.0130, "polished copper"), (0.0068, "nickel")]:
    ax.loglog(dT, rohsenow_boiling_heat_flux(**water, dT_excess=dT, C_sf=C_sf), label=f"{surface}, C_sf = {C_sf}")

# %%
# Up to the critical heat flux
# ----------------------------
#
# The nucleate regime ends at Zuber's critical flux. Setting Rohsenow equal
# to Zuber gives the superheat at burnout for each surface; the critical
# flux itself does not depend on the surface.

q_max = zuber_critical_heat_flux(water["rho_l"], water["rho_v"], water["h_fg"], water["sigma"], C=0.149)
ax.axhline(q_max, color="gray", ls="--", label="critical heat flux (Lienhard-Dhir)")
for C_sf in [0.0060, 0.0130]:
    dT_crit = brentq(lambda d, C_sf=C_sf: rohsenow_boiling_heat_flux(**water, dT_excess=d, C_sf=C_sf) - q_max, 1.0, 100.0)
    print(f"C_sf = {C_sf}: burnout at ΔT_e = {dT_crit:.1f} K")
    assert np.isclose(dT_crit, C_sf / 0.0130 * brentq(lambda d: rohsenow_boiling_heat_flux(**water, dT_excess=d) - q_max, 1.0, 100.0))
ax.set(xlabel="excess temperature ΔT_e [K]", ylabel="heat flux [W/m²]", ylim=(1e3, 3e6), title="Nucleate pool boiling of water at 1 atm")
ax.legend()
fig.tight_layout()

plt.show()

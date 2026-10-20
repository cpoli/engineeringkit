r"""
Van der Waals's equation of state
=================================

In 1869 Thomas Andrews showed that carbon dioxide above 31 °C cannot be
liquefied by any pressure, and below it liquefies at a sharp pressure. In his
1873 thesis Johannes van der Waals explained this with two corrections to the
ideal gas: the molecules fill a volume :math:`b` per mole, and they attract
each other, which lowers the pressure by :math:`a/v^2`:
:math:`p = RT/(v - b) - a/v^2`
(:func:`~engineeringkit.thermodynamics.van_der_waals_pressure`). One
equation then describes both gas and liquid, and has a critical point
(:func:`~engineeringkit.thermodynamics.van_der_waals_critical_point`). He
received the 1910 Nobel Prize in Physics for it.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from engineeringkit.thermodynamics import R_UNIVERSAL, ideal_gas_pressure, van_der_waals_critical_point, van_der_waals_pressure

# %%
# Andrews's carbon dioxide
# ------------------------
#
# With :math:`a = 0.364` Pa m⁶/mol² and :math:`b = 42.67` cm³/mol the
# critical point is :math:`T_c = 8a/(27Rb)`, :math:`p_c = a/(27b^2)` and
# :math:`v_c = 3b`.

a, b = 0.3640, 4.267e-5
Tc, pc, vc = van_der_waals_critical_point(a, b)
print(f"T_c = {Tc:.1f} K ({Tc - 273.15:.1f} °C), p_c = {pc / 1e6:.2f} MPa, v_c = {vc * 1e6:.0f} cm³/mol")
print(f"critical compressibility p_c v_c/(R T_c) = {pc * vc / (R_UNIVERSAL * Tc):.3f} (real CO2: 0.274)")
assert abs(Tc - 304.1) < 0.5 and np.isclose(pc * vc / (R_UNIVERSAL * Tc), 3 / 8)
v_gas = 0.024  # m³/mol, a room-temperature gas
p_ideal, p_vdw = ideal_gas_pressure(1, 300, v_gas), van_der_waals_pressure(1, 300, v_gas, a, b)
print(f"at 300 K and 24 L/mol: ideal {p_ideal:.0f} Pa, van der Waals {p_vdw:.0f} Pa")

# %%
# Maxwell's equal areas
# ---------------------
#
# Below :math:`T_c` an isotherm has a loop where pressure rises with volume,
# which no stable fluid can do. Maxwell (1875) replaced the loop with a
# horizontal line at the vapour pressure, placed so that the two areas it
# cuts off are equal. In reduced variables :math:`p_r = 8T_r/(3v_r - 1) - 3/v_r^2`
# the answer is the same for every van der Waals gas.


def maxwell_line(Tr):
    """Reduced vapour pressure and the liquid and vapour volumes on an isotherm."""

    def volumes(pr):
        roots = np.roots([3 * pr, -(pr + 8 * Tr), 9, -3])
        real = np.sort(roots[np.abs(roots.imag) < 1e-9].real)
        return real[0], real[-1]

    def area_mismatch(pr):
        vl, vg = volumes(pr)
        integral = 8 * Tr / 3 * np.log((3 * vg - 1) / (3 * vl - 1)) + 3 / vg - 3 / vl
        return integral - pr * (vg - vl)

    v = np.linspace(0.4, 5, 4000)
    pr_v = 8 * Tr / (3 * v - 1) - 3 / v**2
    turning = np.where(np.diff(np.sign(np.diff(pr_v))))[0]
    low, high = max(pr_v[turning[0]], 1e-6) + 1e-9, pr_v[turning[1]] - 1e-9
    pr = brentq(area_mismatch, low, high, xtol=1e-14)
    return pr, *volumes(pr)


pr_09, vl_09, vg_09 = maxwell_line(0.9)
print(f"T_r = 0.9: reduced vapour pressure {pr_09:.4f}, liquid v_r {vl_09:.3f}, vapour v_r {vg_09:.3f}")
assert abs(pr_09 - 0.6470) < 1e-3

fig, ax = plt.subplots(figsize=(6, 4))
v_r = np.linspace(0.45, 5, 600)
for Tr in (0.85, 0.9, 0.95, 1.0, 1.1):
    p_r = van_der_waals_pressure(1.0, Tr * Tc, v_r * vc, a, b) / pc
    line = ax.plot(v_r, p_r, lw=1, label=f"T_r = {Tr:g}")[0]
    if Tr < 1:
        pr, vl, vg = maxwell_line(Tr)
        ax.plot([vl, vg], [pr, pr], color=line.get_color(), lw=2)
ax.plot(1, 1, "ko")
ax.set(xlabel="reduced volume v/v_c", ylabel="reduced pressure p/p_c", ylim=(0, 2), title="Isotherms of carbon dioxide")
ax.legend(fontsize=8)
fig.tight_layout()

# %%
# The vapour pressure line ends at the critical point, where liquid and
# vapour become the same. The equal-area pressures fall on a line of
# :math:`\ln p_r` against :math:`1/T_r`, as Clausius-Clapeyron predicts.
# The slope is too shallow: simple real fluids such as argon follow a slope
# near -5.4 (Pitzer's acentric factor measures the departure from it), which
# is why engineering equations of state refine van der Waals's attraction term.

Tr_range = np.linspace(0.7, 0.98, 8)
pr_range = np.array([maxwell_line(Tr)[0] for Tr in Tr_range])
slope = np.polyfit(1 / Tr_range, np.log(pr_range), 1)[0]
print(f"slope of ln p_r against 1/T_r: {slope:.2f}")
assert np.all(np.diff(pr_range) > 0) and -4.5 < slope < -3

plt.show()

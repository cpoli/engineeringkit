r"""
Joule and the first law (1843-1850)
===================================

James Prescott Joule spent the 1840s showing that a given amount of work,
however it is done, always produces the same amount of heat. In his best
known experiment (published 1850) falling weights turned a paddle wheel in
an insulated can of water, and he measured the temperature rise to a
hundredth of a degree. His result, 772 foot-pounds per British thermal
unit, made heat and work two forms of one conserved quantity, energy:
:math:`Q = \Delta U + W`
(:func:`~engineeringkit.thermodynamics.first_law_heat`). The example redoes
the paddle-wheel balance, converts Joule's number to SI, and checks the
ideal-gas energy functions.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import first_law_heat, specific_enthalpy_ideal_gas, specific_internal_energy_ideal_gas

# %%
# The paddle wheel
# ----------------
#
# The can is insulated, so :math:`Q = 0`, and the weights do work *on* the
# water, :math:`W = -mgh` per drop. The first law gives
# :math:`\Delta U = -W`, and with :math:`c = 4186` J/(kg K) the rise in
# temperature. The weights and drop height here are illustrative.

g, c_water = 9.80665, 4186.0
m_weights, drop, n_drops, m_water = 2 * 13.0, 1.6, 20, 6.0
work_on_water = m_weights * g * drop * n_drops
delta_u = work_on_water
delta_T = delta_u / (m_water * c_water)
print(f"work on the water {work_on_water:.0f} J, temperature rise {delta_T:.3f} K")
assert first_law_heat(delta_u, -work_on_water) == 0

drops = np.arange(n_drops + 1)
fig, ax = plt.subplots(figsize=(6, 3.2))
ax.plot(drops, 15.0 + drops * m_weights * g * drop / (m_water * c_water), "o-")
ax.set(xlabel="number of drops", ylabel="water temperature [°C]", title="An adiabatic paddle wheel")
fig.tight_layout()

# %%
# Joule's mechanical equivalent of heat
# -------------------------------------
#
# 772.692 foot-pounds raised one pound of water by one degree Fahrenheit.
# In SI that is joules per calorie (one gram, one kelvin):

ft_lbf, lb, degF = 1.3558179483, 0.45359237, 5 / 9
J_per_cal = 772.692 * ft_lbf / (lb * 1e3 * degF)
print(f"Joule 1850: {J_per_cal:.3f} J/cal; the thermochemical calorie is defined as 4.184 J")
assert abs(J_per_cal - 4.184) / 4.184 < 0.01

# %%
# Energy functions of an ideal gas
# --------------------------------
#
# Joule's free-expansion experiment (1845) found no temperature change when
# a gas expanded into a vacuum, so the internal energy of an ideal gas
# depends on temperature alone. Enthalpy :math:`h = u + pv = u + RT` follows,
# and :math:`c_p - c_v = R`.

R, cv, cp, T = 287.0, 718.0, 1005.0, np.linspace(250, 1500, 6)
u = specific_internal_energy_ideal_gas(T, cv)
h = specific_enthalpy_ideal_gas(T, cp)
np.testing.assert_allclose(h - u, R * T)
print(f"h - u at 300 K: {specific_enthalpy_ideal_gas(300, cp) - specific_internal_energy_ideal_gas(300, cv):.0f} J/kg = R T")

plt.show()

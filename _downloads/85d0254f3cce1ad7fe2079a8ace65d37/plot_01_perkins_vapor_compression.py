r"""
Perkins's vapour-compression refrigerator (1834)
================================================

In 1834 Jacob Perkins, an American engineer in London, patented a machine
that evaporated ether at low pressure to absorb heat, compressed the vapour
with a hand pump, and condensed it again under pressure. It was the first
closed vapour-compression cycle, the cycle of nearly every refrigerator and
air conditioner since. Its figure of merit is the coefficient of performance,
the heat removed per unit of work,
:math:`\mathrm{COP}_R = Q_c/W`
(:func:`~engineeringkit.thermodynamics.refrigeration_cop`), and the second
law bounds it by :math:`T_c/(T_h - T_c)`
(:func:`~engineeringkit.thermodynamics.carnot_cop_refrigerator`). The example
runs the ideal cycle with water as the refrigerant (R718, used today in some
large chillers), because IAPWS-IF97 gives its properties.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import carnot_cop_refrigerator, refrigeration_cop, saturation_pressure, steam_properties
from engineeringkit.thermodynamics.visualizers import plot_ts_diagram


def cycle(T_evap, T_cond):
    """States 1-4: compressor inlet, compressor outlet, condenser outlet, evaporator inlet."""
    p_low, p_high = float(saturation_pressure(T_evap)), float(saturation_pressure(T_cond))
    s1 = steam_properties(p_low, x=1.0)  # saturated vapour leaves the evaporator
    s2 = steam_properties(p_high, s=s1.s)  # isentropic compression
    s3 = steam_properties(p_high, x=0.0)  # saturated liquid leaves the condenser
    s4 = steam_properties(p_low, h=s3.h)  # the expansion valve throttles at constant enthalpy
    return s1, s2, s3, s4


# %%
# A chiller at 5 °C rejecting heat at 35 °C
# -----------------------------------------

s1, s2, s3, s4 = cycle(278.15, 308.15)
q_cold, work = s1.h - s4.h, s2.h - s1.h
cop = refrigeration_cop(q_cold, work)
cop_carnot = carnot_cop_refrigerator(308.15, 278.15)
print(f"evaporator {s1.p:.0f} Pa, condenser {s2.p:.0f} Pa; compressor outlet {s2.T - 273.15:.0f} °C")
print(f"refrigeration effect {q_cold / 1e3:.0f} kJ/kg, compressor work {work / 1e3:.0f} kJ/kg, quality after the valve {s4.x:.3f}")
print(f"COP {cop:.2f}, Carnot COP {cop_carnot:.2f}, ratio {cop / cop_carnot:.2f}")
print(f"vapour volume at the compressor inlet {s1.v:.0f} m³/kg: {q_cold / s1.v / 1e3:.1f} kJ of cooling per m³ swept")
assert np.isclose(s4.h, s3.h) and 0 < s4.x < 0.1
assert 0.8 < cop / cop_carnot < 1  # throttling and superheat cost a little

# %%
# The cycle stays below Carnot for two reasons. The valve throttles instead of
# recovering work from an expander, and the compressed vapour leaves far hotter
# than the condenser. Water's real handicap is its huge vapour volume at low
# pressure: at 5 °C a compressor must sweep about 150 m³ per kilogram,
# against about 0.24 m³ for ammonia, which gives more than 200 times as much
# cooling per cubic metre swept. That is why practical machines used ether,
# then ammonia, and later fluorocarbons.

s_path = [s1.s, s2.s] + [st.s for st in (steam_properties(s2.p, T) for T in np.linspace(s2.T, s3.T + 0.01, 30))]
T_path = [s1.T, s2.T] + list(np.linspace(s2.T, s3.T + 0.01, 30))
s_path += [s3.s, s4.s, s1.s]
T_path += [s3.T, s4.T, s1.T]
ax = plot_ts_diagram(s_path, T_path, label="ideal vapour-compression cycle")
ax.set(xlim=(0, 10), ylim=(0, 400), title="Water as refrigerant, 5 °C to 35 °C")

# %%
# A warmer condenser costs efficiency
# -----------------------------------
#
# Each kelvin of lift between evaporator and condenser lowers the COP, which
# is why refrigerators work harder in a hot kitchen.

T_cond = np.linspace(293.15, 333.15, 21)
cops = []
for T in T_cond:
    a, b, c, d = cycle(278.15, T)
    cops.append(refrigeration_cop(a.h - d.h, b.h - a.h))
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(T_cond - 273.15, cops, label="ideal vapour-compression cycle")
ax.plot(T_cond - 273.15, [carnot_cop_refrigerator(T, 278.15) for T in T_cond], "--", label="Carnot")
ax.set(xlabel="condensing temperature [°C]", ylabel="COP", title="Evaporator at 5 °C")
ax.legend()
fig.tight_layout()
assert np.all(np.diff(cops) < 0)

plt.show()

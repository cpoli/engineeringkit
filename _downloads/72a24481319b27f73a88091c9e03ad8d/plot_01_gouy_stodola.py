r"""
Gouy and Stodola: lost work
===========================

Georges Gouy (1889) asked how much work a system can deliver when the only
other body available is the environment at :math:`T_0`. He showed that every
irreversibility costs work equal to :math:`T_0` times the entropy it
generates. Aurel Stodola reached the same result independently for steam
turbines, in the editions of *Die Dampfturbinen* from 1898 to 1905. The
Gouy-Stodola theorem, :math:`W_{lost} = T_0 S_{gen}`
(:func:`~engineeringkit.thermodynamics.exergy_destruction`), lets an
engineer rank the losses of a plant in a single currency. The example prices
three irreversibilities: heat crossing a temperature difference, a
throttling valve, and an imperfect turbine.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import carnot_efficiency, exergy_destruction, steam_properties

T0 = 298.15

# %%
# Heat across a temperature difference
# ------------------------------------
#
# One kilowatt flows from a furnace at 1200 K into boiling water at 600 K.
# The entropy generated is :math:`\dot Q(1/T_2 - 1/T_1)`. The lost work equals
# the difference between what a Carnot engine could have made from the heat
# at each temperature.

Q, T1, T2 = 1000.0, 1200.0, 600.0
S_gen = Q / T2 - Q / T1
lost = exergy_destruction(T0, S_gen)
print(f"heat transfer: S_gen = {S_gen:.4f} W/K, lost work {lost:.0f} W")
assert np.isclose(lost, Q * (carnot_efficiency(T1, T0) - carnot_efficiency(T2, T0)))

# %%
# A throttling valve
# ------------------
#
# Steam at 3 MPa and 300 °C is throttled to 0.5 MPa. Its enthalpy is
# unchanged, so the first law sees no loss at all. The second law does: the
# entropy rises, and the drop in flow exergy :math:`h - T_0 s` equals the lost
# work.

inlet = steam_properties(3e6, 573.15)
outlet = steam_properties(0.5e6, h=inlet.h)
s_gen = outlet.s - inlet.s
lost_valve = exergy_destruction(T0, s_gen)
exergy_drop = (inlet.h - T0 * inlet.s) - (outlet.h - T0 * outlet.s)
print(f"valve: {inlet.T - 273.15:.0f} °C to {outlet.T - 273.15:.1f} °C, s_gen = {s_gen:.1f} J/(kg K)")
print(f"lost work {lost_valve / 1e3:.1f} kJ/kg, drop in flow exergy {exergy_drop / 1e3:.1f} kJ/kg")
assert np.isclose(lost_valve, exergy_drop)

# %%
# A turbine with 85% isentropic efficiency
# ----------------------------------------
#
# The same steam expands to 10 kPa. The real turbine delivers 85% of the
# isentropic enthalpy drop. Its lost work is less than the 15% of work it
# misses: part of the extra enthalpy left in the exhaust could still be used.

ideal = steam_properties(10e3, s=inlet.s)
w_ideal = inlet.h - ideal.h
real = steam_properties(10e3, h=inlet.h - 0.85 * w_ideal)
lost_turbine = exergy_destruction(T0, real.s - inlet.s)
print(f"turbine: work {0.85 * w_ideal / 1e3:.0f} of {w_ideal / 1e3:.0f} kJ/kg, lost work {lost_turbine / 1e3:.0f} kJ/kg")
print(f"exhaust quality {real.x:.3f} instead of {ideal.x:.3f}")
assert lost_turbine < 0.15 * w_ideal
assert np.isclose(lost_turbine, (inlet.h - T0 * inlet.s) - (real.h - T0 * real.s) - 0.85 * w_ideal)

# %%
# Lost work against isentropic efficiency
# ---------------------------------------

eta = np.linspace(0.6, 1.0, 41)
lost_curve = [exergy_destruction(T0, steam_properties(10e3, h=inlet.h - e * w_ideal).s - inlet.s) / 1e3 for e in eta]
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(eta, lost_curve, label="lost work, T₀ s_gen")
ax.plot(eta, (1 - eta) * w_ideal / 1e3, "--", label="work not delivered")
ax.set(xlabel="isentropic efficiency", ylabel="[kJ/kg]", title="Turbine from 3 MPa, 300 °C to 10 kPa")
ax.legend()
fig.tight_layout()
assert np.isclose(lost_curve[-1], 0, atol=1e-6) and np.all(np.diff(lost_curve) < 0)

plt.show()

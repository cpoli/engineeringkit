r"""
Rankine's steam cycle
=====================

William Rankine's *Manual of the Steam Engine* (1859) analysed the steam
engine as a cycle of four steady-flow steps: a pump raises condensed water
to boiler pressure, the boiler turns it into superheated steam, the steam
expands through the engine (now a turbine), and a condenser turns the
exhaust back into water. Because the water is pumped as a liquid, the pump
needs about 1% of the turbine work, which is the cycle's great advantage
over gas cycles. The efficiency
(:func:`~engineeringkit.thermodynamics.rankine_ideal_efficiency`) is
computed here with IAPWS-IF97 steam properties
(:func:`~engineeringkit.thermodynamics.steam_properties`) for a plant at
10 MPa and 500 °C condensing at 10 kPa.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import carnot_efficiency, rankine_ideal_efficiency, steam_properties
from engineeringkit.thermodynamics.visualizers import plot_ts_diagram


def rankine(p_boiler, T_boiler, p_condenser):
    """States 1-4 of the ideal cycle: pump outlet, turbine inlet, turbine outlet, condenser outlet."""
    s4 = steam_properties(p_condenser, x=0.0)
    s1 = steam_properties(p_boiler, s=s4.s)  # isentropic pump
    s2 = steam_properties(p_boiler, T_boiler)
    s3 = steam_properties(p_condenser, s=s2.s)  # isentropic turbine
    return s1, s2, s3, s4


# %%
# The four states
# ---------------

s1, s2, s3, s4 = rankine(10e6, 773.15, 10e3)
for name, st in zip(("1 pump outlet", "2 turbine inlet", "3 turbine outlet", "4 condenser outlet"), (s1, s2, s3, s4)):
    print(f"{name:20s} T = {st.T - 273.15:6.1f} °C, h = {st.h / 1e3:7.1f} kJ/kg, s = {st.s / 1e3:.4f} kJ/(kg K)")
eta = rankine_ideal_efficiency(s1.h, s2.h, s3.h, s4.h)
w_pump, w_turbine = s1.h - s4.h, s2.h - s3.h
print(f"quality at turbine exit {s3.x:.3f}; back-work ratio {w_pump / w_turbine:.2%}; efficiency {eta:.3f}")
print(f"Carnot between {s2.T - 273.15:.0f} °C and {s4.T - 273.15:.1f} °C: {carnot_efficiency(s2.T, s4.T):.3f}")
assert round(eta, 3) == 0.402 and round(s3.x, 3) == 0.793  # Cengel & Boles, Example 10-3
assert np.isclose(w_pump, s4.v * (10e6 - 10e3), rtol=0.01)

# %%
# The cycle on the T-s diagram
# ----------------------------
#
# Heat goes in along 1-2: first warming cold feedwater, then boiling at
# 311 °C, then superheating. Most of it enters far below the 500 °C peak,
# which is why the cycle falls short of Carnot.

T_boil = steam_properties(10e6, x=0.0).T
path = (
    [steam_properties(10e6, T) for T in np.linspace(s1.T, T_boil - 0.01, 60)]
    + [steam_properties(10e6, x=x) for x in np.linspace(0, 1, 20)]
    + [steam_properties(10e6, T) for T in np.linspace(T_boil + 0.01, s2.T, 60)]
)
s_cycle = [s4.s] + [st.s for st in path] + [s3.s, s4.s]
T_cycle = [s4.T] + [st.T for st in path] + [s3.T, s4.T]
ax = plot_ts_diagram(s_cycle, T_cycle, label="ideal Rankine cycle")
ax.set_title("10 MPa, 500 °C to 10 kPa")

# %%
# Raising the boiler pressure
# ---------------------------
#
# A higher boiler pressure boils the water at a higher temperature and raises
# the efficiency, but at fixed 500 °C it leaves wetter steam in the last
# turbine stages, which erodes the blades. Engineers keep the exit quality
# above about 0.88 by superheating more or by reheating.

pressures = np.linspace(1e6, 16e6, 16)
cycles = [rankine(p, 773.15, 10e3) for p in pressures]
etas = [rankine_ideal_efficiency(c[0].h, c[1].h, c[2].h, c[3].h) for c in cycles]
qualities = [c[2].x for c in cycles]
fig, ax1 = plt.subplots(figsize=(6, 3.5))
ax1.plot(pressures / 1e6, etas, label="efficiency")
ax1.plot(pressures / 1e6, qualities, label="turbine exit quality")
ax1.set(xlabel="boiler pressure [MPa]", ylim=(0.2, 1.0), title="Steam at 500 °C, condenser at 10 kPa")
ax1.legend()
fig.tight_layout()
assert np.all(np.diff(etas) > 0) and np.all(np.diff(qualities) < 0)

plt.show()

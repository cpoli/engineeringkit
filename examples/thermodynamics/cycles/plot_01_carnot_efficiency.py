r"""
Carnot's *Réflexions*: no engine beats the reversible one (1824)
================================================================

Sadi Carnot (1824) argued that a heat engine's work comes from heat falling
from a hot body to a cold one, as a water wheel's comes from falling water,
and that a reversible engine between two temperatures is the best possible.
If some engine did better, it could drive the reversible engine backwards
as a heat pump and, with no work left over, move heat from cold to hot. In
modern form the bound is :math:`\eta_C = 1 - T_c/T_h`
(:func:`~engineeringkit.thermodynamics.carnot_efficiency`), with absolute
temperatures. The example runs Carnot's argument with numbers and shows how
the bound rises with the boiler temperature.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import carnot_efficiency, heat_engine_work

# %%
# Carnot's argument with numbers
# ------------------------------
#
# A claimed "super engine" between 600 K and 300 K takes 1000 J and delivers
# 600 J of work, more than the reversible 500 J. Use all of that work to run
# a reversible engine backwards: as a heat pump it delivers
# :math:`W/\eta_C` to the hot reservoir.

T_hot, T_cold = 600.0, 300.0
eta_C = carnot_efficiency(T_hot, T_cold)
Q_in, W_super = 1000.0, 600.0
Q_out_super = Q_in - W_super
assert np.isclose(heat_engine_work(Q_in, Q_out_super), W_super)

Q_hot_pumped = W_super / eta_C
Q_cold_pumped = Q_hot_pumped - W_super
net_to_hot = Q_hot_pumped - Q_in
net_from_cold = Q_cold_pumped - Q_out_super
print(f"Carnot efficiency {eta_C:.2f}; super engine {W_super / Q_in:.2f}")
print(f"net heat into the hot reservoir {net_to_hot:.0f} J, taken from the cold one {net_from_cold:.0f} J, net work 0 J")
assert np.isclose(net_to_hot, net_from_cold) and net_to_hot > 0

# %%
# The combined machine moves 200 J from cold to hot with no outside work,
# which never happens. So no engine beats :math:`\eta_C`, and all
# reversible engines between the same temperatures have the same efficiency,
# whatever the working fluid.

# %%
# The bound against boiler temperature
# ------------------------------------
#
# With heat rejected at 300 K, raising the source temperature helps less and
# less. Carnot's own conclusion was that steam engines should work over the
# widest possible temperature range, which is why boiler temperatures and
# pressures have been rising ever since.

T_h = np.linspace(350, 1600, 200)
eta = np.array([carnot_efficiency(T, T_cold) for T in T_h])
plants = {"Newcomen engine (373 K)": 373.0, "subcritical steam (813 K)": 813.0, "gas turbine (1600 K)": 1600.0}
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(T_h, eta)
for name, T in plants.items():
    ax.plot(T, carnot_efficiency(T, T_cold), "o", label=f"{name}: {carnot_efficiency(T, T_cold):.2f}")
ax.set(xlabel="hot reservoir temperature [K]", ylabel="Carnot efficiency", title="Heat rejected at 300 K", ylim=(0, 1))
ax.legend(fontsize=8)
fig.tight_layout()
assert np.all(np.diff(eta) > 0) and np.all(np.diff(eta, 2) < 0)

plt.show()

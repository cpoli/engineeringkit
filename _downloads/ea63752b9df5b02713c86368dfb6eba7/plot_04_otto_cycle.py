r"""
Otto's four-stroke cycle
========================

Nicolaus Otto's 1876 engine compressed its charge before igniting it, a
principle Alphonse Beau de Rochas had described in 1862. Its four strokes
(intake, compression, power, exhaust) are idealized as the air-standard
Otto cycle: isentropic compression through a ratio :math:`r`, heat added at
constant volume, isentropic expansion and heat rejected at constant volume.
The efficiency :math:`1 - r^{1-\gamma}`
(:func:`~engineeringkit.thermodynamics.otto_efficiency`,
:class:`~engineeringkit.thermodynamics.OttoCycle`) depends only on the
compression ratio. Otto's first engine had :math:`r \approx 2.5`; knock
limits petrol engines to about 10-12 today.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import OttoCycle, isentropic_pressure_volume, otto_efficiency, thermal_efficiency
from engineeringkit.thermodynamics.visualizers import plot_pv_diagram

# %%
# The cycle on the p-v diagram
# ----------------------------
#
# One kilogram of air at 100 kPa and 300 K, compression ratio 8, 800 kJ/kg
# added by combustion.

R, gamma, r = 287.0, 1.4, 8.0
cv = R / (gamma - 1)
p1, T1, q_in = 100e3, 300.0, 800e3
v1 = R * T1 / p1
v2 = v1 / r
p2 = isentropic_pressure_volume(p1, v1, v2, gamma)
T2 = p2 * v2 / R
T3 = T2 + q_in / cv
p3 = R * T3 / v2
p4 = isentropic_pressure_volume(p3, v2, v1, gamma)
T4 = p4 * v1 / R
q_out = cv * (T4 - T1)

compress = np.linspace(v1, v2, 60)
expand = np.linspace(v2, v1, 60)
v_path = np.concatenate([compress, expand, [v1]])
p_path = np.concatenate([isentropic_pressure_volume(p1, v1, compress, gamma), isentropic_pressure_volume(p3, v2, expand, gamma), [p1]])
ax = plot_pv_diagram(v_path, p_path, label=f"r = {r:g}")
ax.set(xlabel="specific volume [m$^3$/kg]", title="Air-standard Otto cycle")

cycle = OttoCycle(r)
print(cycle.summary())
print(f"peak temperature {T3:.0f} K, peak pressure {p3 / 1e6:.1f} MPa")
print(f"efficiency from the four states {thermal_efficiency(q_in, q_out):.4f}")
assert np.isclose(thermal_efficiency(q_in, q_out), cycle.efficiency)
assert np.isclose(T2 / T1, T3 / T4)

# %%
# The compression ratio is everything
# -----------------------------------
#
# The efficiency does not depend on how much fuel is burned, only on how far
# the charge is compressed. Higher octane fuels resist knock and allow higher
# compression ratios.

ratios = np.linspace(1, 16, 151)
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(ratios, otto_efficiency(ratios))
for label, rr in (("Otto 1876", 2.5), ("1950s car", 7.0), ("modern petrol", 11.0)):
    ax.plot(rr, otto_efficiency(rr), "o", label=f"{label}: r = {rr:g}, η = {otto_efficiency(rr):.2f}")
ax.set(xlabel="compression ratio", ylabel="ideal efficiency", ylim=(0, 0.7))
ax.legend(fontsize=8)
fig.tight_layout()

plt.show()

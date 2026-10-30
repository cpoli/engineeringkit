r"""
Perkins's cycle on a modern refrigerant (1834)
==============================================

Perkins's cycle of evaporation, compression, condensation and throttling
is unchanged. The working fluid is: ether gave way to ammonia, then to the
fluorocarbons, and R-134a ran most car air conditioners and domestic
refrigerators from the 1990s. With its properties from the Peng-Robinson
equation of state
(:func:`~engineeringkit.thermodynamics.refrigerant_properties`),
:func:`~engineeringkit.thermodynamics.vapor_compression_cycle` solves the
cycle state by state.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import carnot_cop_refrigerator, vapor_compression_cycle

# %%
# Çengel and Boles's refrigerator
# -------------------------------
#
# Example 11-1: R-134a between 0.14 and 0.8 MPa at 0.05 kg/s. The tables
# give 7.18 kW of cooling for 1.81 kW of work, a COP of 3.97.

ideal = vapor_compression_cycle(0.14e6, 0.8e6)
for label, st in zip(ideal.labels, ideal.states):
    print(f"{label:18s} {st.p / 1e3:6.0f} kPa {st.T - 273.15:6.1f} °C  h = {st.h / 1e3:6.1f} kJ/kg  ({st.phase})")
print(
    f"cooling {0.05 * ideal.q_evaporator / 1e3:.2f} kW, work {0.05 * ideal.compressor_work / 1e3:.2f} kW, COP {ideal.cop_refrigerator:.2f}"
)
assert abs(ideal.cop_refrigerator - 3.97) < 0.08

real = vapor_compression_cycle(0.14e6, 0.8e6, superheat=5.0, subcooling=3.0, compressor_efficiency=0.8)
carnot = carnot_cop_refrigerator(ideal.states[2].T, ideal.states[0].T)
print(f"with 5 K superheat, 3 K subcooling and an 80% compressor: COP {real.cop_refrigerator:.2f}; Carnot {carnot:.2f}")
assert real.cop_refrigerator < ideal.cop_refrigerator < carnot
fig, ax = ideal.plot()

# %%
# A warmer kitchen
# ----------------
#
# Raising the condensing pressure, as on a hot day, lowers the COP.

p_cond = np.linspace(0.6e6, 1.6e6, 11)
cops = [vapor_compression_cycle(0.14e6, p).cop_refrigerator for p in p_cond]
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(p_cond / 1e6, cops, "o-")
ax.set(xlabel="condenser pressure [MPa]", ylabel="COP", title="R-134a, evaporator at 0.14 MPa")
fig.tight_layout()
assert np.all(np.diff(cops) < 0)

plt.show()

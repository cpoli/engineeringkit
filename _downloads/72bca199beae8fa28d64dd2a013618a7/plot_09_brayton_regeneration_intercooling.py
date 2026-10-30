r"""
Brayton's cycle with regeneration, intercooling and reheat (1872)
=================================================================

A simple gas turbine throws its hot exhaust away and spends much of its
turbine work on the compressor. A regenerator heats the compressed air with
the exhaust. Splitting the compression with intercoolers and the expansion
with reheaters moves both towards isothermal processes. With many stages and
a perfect regenerator the cycle becomes Ericsson's, whose efficiency is
Carnot's (:func:`~engineeringkit.thermodynamics.brayton_cycle`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import brayton_cycle, brayton_efficiency, carnot_efficiency

# %%
# A gas turbine with imperfect components
# ---------------------------------------
#
# Air at 300 K, turbine inlet 1300 K, pressure ratio 8, an 80% compressor and
# an 85% turbine (the cold-air-standard form of Çengel and Boles's Example 9-6).

simple = brayton_cycle(8.0, 300.0, 1300.0, compressor_efficiency=0.8, turbine_efficiency=0.85)
regenerated = brayton_cycle(8.0, 300.0, 1300.0, compressor_efficiency=0.8, turbine_efficiency=0.85, regenerator_effectiveness=0.8)
print(f"ideal {brayton_efficiency(8.0):.3f}; with losses {simple.efficiency:.3f}, back-work ratio {simple.back_work_ratio:.2f}")
print(f"80% effective regenerator: {regenerated.efficiency:.3f}, saving {regenerated.heat_regenerated / 1e3:.0f} kJ/kg of fuel heat")
assert regenerated.efficiency > simple.efficiency and np.isclose(regenerated.net_work, simple.net_work)

staged = brayton_cycle(
    8.0,
    300.0,
    1300.0,
    compressor_efficiency=0.8,
    turbine_efficiency=0.85,
    regenerator_effectiveness=0.8,
    compressor_stages=2,
    turbine_stages=2,
)
print(f"two-stage compression and expansion with regeneration: {staged.efficiency:.3f}")
fig, ax = staged.plot()

# %%
# Towards Ericsson and Carnot
# ---------------------------
#
# With ideal components and a perfect regenerator, more stages close the gap
# to Carnot's bound between 300 K and 1300 K.

stages = [1, 2, 3, 5, 10, 30, 100]
etas = [brayton_cycle(8.0, 300.0, 1300.0, regenerator_effectiveness=1.0, compressor_stages=n, turbine_stages=n).efficiency for n in stages]
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.semilogx(stages, etas, "o-", label="ideal, perfect regenerator")
ax.axhline(carnot_efficiency(1300.0, 300.0), color="gray", ls=":", label="Carnot (Ericsson limit)")
ax.set(xlabel="stages of compression and of expansion", ylabel="efficiency", title="Intercooling and reheat, pressure ratio 8")
ax.legend()
fig.tight_layout()
assert np.all(np.diff(etas) > 0) and carnot_efficiency(1300.0, 300.0) - etas[-1] < 0.01

plt.show()

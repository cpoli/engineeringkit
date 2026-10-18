r"""
Duhamel's thermal stress (1838)
===============================

Jean-Marie Duhamel (1838) added temperature to Navier's equations of
elasticity: a change :math:`\Delta T` produces a free strain
:math:`\alpha\,\Delta T` (:func:`~engineeringkit.materials.thermal_strain`)
that causes no stress if the body can expand, and a stress
:math:`-E\alpha\,\Delta T` (:func:`~engineeringkit.materials.constrained_thermal_stress`)
if it is held. The stress does not depend on the length of the bar, which is
why continuously welded rail must be laid at a chosen neutral temperature.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import constrained_thermal_stress, thermal_strain

# %%
# Continuously welded rail
# ------------------------
#
# Rail steel has :math:`E = 210` GPa and :math:`\alpha = 11.5\times10^{-6}`/K.
# Laid stress-free at 25 °C, it is in compression on hot days and in
# tension in winter. The total strain, free plus elastic, stays zero.

E, alpha, neutral = 210e9, 11.5e-6, 25.0
T = np.linspace(-20, 60, 81)
stress = constrained_thermal_stress(E, alpha, T - neutral)
print(f"at 60 C: {stress[-1] / 1e6:.1f} MPa; at -20 C: {stress[0] / 1e6:.1f} MPa")
assert np.allclose(thermal_strain(alpha, T - neutral) + stress / E, 0)

# %%
# Jointed rail: a gap delays the stress
# -------------------------------------
#
# Bolted 18 m rails leave a gap of a few millimetres at each joint. The rail
# expands freely until the gap closes, and only the excess strain is
# restrained.

length = 18.0
fig, ax = plt.subplots(figsize=(6, 3.8))
ax.plot(T, stress / 1e6, label="welded (no gap)")
for gap in (0.003, 0.006):
    excess = np.maximum(thermal_strain(alpha, T - neutral) - gap / length, 0)
    jointed = -E * excess
    closes = neutral + gap / (alpha * length)
    ax.plot(T, jointed / 1e6, label=f"{1e3 * gap:.0f} mm gap, closes at {closes:.0f} C")
    assert np.all(jointed[T <= closes] == 0)
    assert np.allclose(jointed[T > closes], constrained_thermal_stress(E, alpha, T[T > closes] - closes))
ax.axhline(0, color="0.6", lw=0.8)
ax.set(xlabel="rail temperature [C]", ylabel="axial stress [MPa]", title="Thermal stress in rail")
ax.legend()
fig.tight_layout()

plt.show()

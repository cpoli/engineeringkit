r"""
Ramberg and Osgood's stress-strain curve (1943)
===============================================

Aluminium alloys have no sharp yield point: their stress-strain curves bend
over smoothly, and aircraft designers of the Second World War needed them
in a form they could compute with. Walter Ramberg and William Osgood (NACA,
1943) described the whole curve with three numbers,
:math:`\varepsilon = \sigma/E + K(\sigma/E)^n`, usually written today with
the 0.2% offset yield strength,
:math:`\varepsilon = \sigma/E + 0.002(\sigma/\sigma_y)^n`
(:func:`~engineeringkit.materials.ramberg_osgood_strain`). The example
plots the family and recovers the offset yield strength from a sampled
curve.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import StressStrainCurve, ramberg_osgood_strain

# %%
# Ramberg and Osgood's form and the offset form
# ---------------------------------------------
#
# The two forms are the same curve with :math:`K = 0.002(E/\sigma_y)^n`.

E, sy = 72e9, 345e6
stress = np.linspace(0, 1.15 * sy, 400)
n = 15
K = 0.002 * (E / sy) ** n
assert np.allclose(ramberg_osgood_strain(stress, E, sy, n), stress / E + K * (stress / E) ** n)

# %%
# The 0.2% offset yield strength
# ------------------------------
#
# Draw the elastic line shifted by 0.002 strain; it crosses the curve at
# :math:`\sigma_y`. :meth:`~engineeringkit.materials.StressStrainCurve.offset_yield_strength`
# does this on sampled data. Smaller :math:`n` means more gradual yielding
# and more hardening after it.

fig, ax = plt.subplots(figsize=(6, 3.8))
for exponent in (5, 10, 20, 50):
    curve = StressStrainCurve(ramberg_osgood_strain(stress, E, sy, exponent), stress)
    found = curve.offset_yield_strength()
    print(f"n = {exponent:2d}: offset yield {found / 1e6:.1f} MPa")
    assert np.isclose(found, sy, rtol=5e-3)
    curve.plot(ax=ax)
    ax.lines[-1].set_label(f"n = {exponent}")
offset = np.linspace(0.002, 0.002 + 1.15 * sy / E, 2)
ax.plot(offset, E * (offset - 0.002) / 1e6, "k:", lw=1, label="0.2% offset line")
ax.set(xlim=(0, 0.012), title="Ramberg-Osgood curves, 2024-T3 aluminium")
ax.legend()
fig.tight_layout()

plt.show()

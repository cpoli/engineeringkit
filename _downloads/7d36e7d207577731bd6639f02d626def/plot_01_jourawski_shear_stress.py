r"""
Jourawski's shear stress (1856)
===============================

Building timber railway bridges in Russia, Dmitrii Ivanovich Jourawski saw
large oak beams split along their length at mid-depth, where the bending
stress is zero. Balancing the change of bending force on the part of the
section above a cut, he found the shear stress on it,
:math:`\tau = VQ/(It)` (:func:`~engineeringkit.structural_engineering.transverse_shear_stress`),
and the force per unit length the keys or bolts of a built-up beam must
carry, :math:`q = VQ/I` (:func:`~engineeringkit.structural_engineering.shear_flow`).
The example draws the shear stress in a rectangle and an I-beam, checks
that it adds up to the shear force, and spaces the bolts of a flitch beam,
timber with a steel plate, using its transformed section
(:func:`~engineeringkit.structural_engineering.composite_section`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import composite_section, i_section, shear_flow, transverse_shear_stress

# %%
# A rectangle and an I-beam
# -------------------------

V = 60e3
rect = composite_section([(0.15, 0.3, 0.0, 0.15)])
beam = i_section(0.3, 0.15, 0.012, 0.008)
fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
for ax, section, name in zip(axes, (rect, beam), ("rectangle", "I-beam")):
    y = np.linspace(0, 0.3, 601)
    width = np.array([section.width_at(v) for v in y])
    tau = np.array([transverse_shear_stress(V, section.first_moment(v), section.Iz, w) for v, w in zip(y, width)])
    total = np.trapezoid(tau * width, y)
    print(f"{name}: peak {tau.max() / 1e6:.2f} MPa, integral of tau t dy = {total / 1e3:.2f} kN")
    assert np.isclose(total, V, rtol=2e-3)
    ax.plot(tau / 1e6, y)
    ax.set(xlabel="shear stress [MPa]", title=name)
axes[0].set_ylabel("height [m]")
assert np.isclose(transverse_shear_stress(V, rect.first_moment(0.15), rect.Iz, 0.15), 1.5 * V / rect.area)
fig.tight_layout()

# %%
# Bolts in a flitch beam
# ----------------------
#
# Two timber joists (E = 10 GPa) are bolted to a steel plate (200 GPa)
# between them, with a timber cover board glued on top. In timber units the
# plate is 20 times wider.

n = 200e9 / 10e9
flitch = composite_section([(0.05, 0.25, -0.035, 0.125), (0.01, 0.25, 0.0, 0.125, n), (0.05, 0.25, 0.035, 0.125), (0.17, 0.04, 0.0, 0.27)])
q = shear_flow(V, flitch.first_moment(0.25), flitch.Iz)
capacity = 6e3  # one glued screw [N]
print(f"neutral axis {flitch.y_centroid * 1e3:.1f} mm, shear flow under the cover {q / 1e3:.1f} kN/m")
print(f"screw spacing at most {capacity / q * 1e3:.0f} mm")
assert np.isclose(q, V * 0.17 * 0.04 * (0.27 - flitch.y_centroid) / flitch.Iz)

fig, ax = flitch.plot()
ax.set_title("Flitch beam: plate drawn at its real width")

plt.show()

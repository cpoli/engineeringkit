r"""
Brinell's hardness test (1900)
==============================

Johan August Brinell, chief engineer of the Fagersta ironworks in Sweden,
showed his ball test at the 1900 Paris Exposition. A 10 mm hardened steel
ball pressed into the part with 3000 kgf leaves a dent of diameter
:math:`d`; the load divided by the curved area of the dent is the
hardness number (:func:`~engineeringkit.materials.brinell_hardness`),

.. math:: \mathrm{HB} = \frac{2P}{\pi D\left(D - \sqrt{D^2 - d^2}\right)}.

Indentation is plastic flow under a constrained punch, so for steels HB
tracks the tensile strength, :math:`S_{ut} \approx 3.45\,\mathrm{HB}` MPa
(:func:`~engineeringkit.materials.brinell_tensile_strength`): a test that
leaves only a small mark estimates the strength of a finished part.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import brinell_hardness, brinell_tensile_strength, vickers_hardness

# %%
# From the dent to the number
# ---------------------------
#
# A 3000 kgf load and a dent 4.0 mm across give about 229 HB, the hardness
# of a normalised medium-carbon steel. The projected-area pressure,
# :math:`4P/\pi d^2`, is a few per cent higher than HB because the
# curved cap has more area than its shadow.

P, D = 3000 * 9.80665, 10e-3
HB = brinell_hardness(P, D, 4.0e-3)
mean_pressure = 4 * (P / 9.80665) / (np.pi * 4.0**2)
print(f"HB = {HB:.0f}, projected mean pressure {mean_pressure:.0f} kgf/mm^2, Sut ~ {brinell_tensile_strength(HB) / 1e6:.0f} MPa")
assert 1 < mean_pressure / HB < 1.05

# %%
# Hardness against strength
# -------------------------
#
# Dent diameters in Brinell's valid range, :math:`0.24D` to :math:`0.6D`,
# cover soft to hardened steels. The right-hand axis reads the strength
# estimate. A Vickers diamond, whose dents are geometrically similar at any
# load, gives a number close to HB on the same steel below about 400 HB.

d = np.linspace(2.4e-3, 6.0e-3, 50)
hb = brinell_hardness(P, D, d)
sut = brinell_tensile_strength(hb)
print(f"range {hb[-1]:.0f} to {hb[0]:.0f} HB, Sut {sut[-1] / 1e6:.0f} to {sut[0] / 1e6:.0f} MPa")
hv = vickers_hardness(30 * 9.80665, 0.50e-3)
print(f"Vickers, 30 kgf, 0.50 mm diagonal: {hv:.0f} HV")
assert np.all(np.diff(hb) < 0) and np.allclose(sut / hb, 3.45e6)

fig, ax = plt.subplots(figsize=(6, 3.8))
ax.plot(1e3 * d, hb)
ax.set(xlabel="dent diameter [mm]", ylabel="Brinell hardness HB", title="3000 kgf, 10 mm ball")
ax.secondary_yaxis("right", functions=(lambda h: 3.45 * h, lambda s: s / 3.45)).set_ylabel("estimated Sut [MPa]")
ax.grid(True, alpha=0.3)
fig.tight_layout()

plt.show()

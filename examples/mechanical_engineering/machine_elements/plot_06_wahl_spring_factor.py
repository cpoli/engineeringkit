r"""
Wahl's curvature factor for helical springs (1929)
==================================================

The wire of a helical spring is a bar in torsion, twisted by the torque
:math:`FD/2`, so the nominal shear stress is :math:`8FD/\pi d^3` and the
rate :math:`k = Gd^4/8D^3N`
(:func:`~engineeringkit.mechanical_engineering.spring_rate`). Springs broke
in fatigue at loads this predicted to be safe. Arthur Wahl, at
Westinghouse, analysed the wire as a curved bar: on the inside of the coil
the fibres are shorter, so the same twist strains them more, and the direct
shear of the load adds on that side too. His factor
:math:`K_W = (4C - 1)/(4C - 4) + 0.615/C`
(:func:`~engineeringkit.mechanical_engineering.wahl_factor`), with
:math:`C = D/d`, multiplies the nominal stress
(:func:`~engineeringkit.mechanical_engineering.spring_shear_stress`). The
1929 paper is cited from Wahl's *Mechanical Springs* (1944) and Shigley.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import bergstrasser_factor, spring_index, spring_rate, spring_shear_stress, wahl_factor

# %%
# Rate from the strain energy of the wire
# ---------------------------------------
#
# :math:`U = T^2 L/2GJ` with :math:`T = FD/2`, :math:`L = \pi DN` and
# :math:`J = \pi d^4/32`; Castigliano's theorem gives the deflection
# :math:`\partial U/\partial F`.

d, D, N, G, F = 0.003, 0.024, 10, 79.3e9, 60.0
J, length = np.pi * d**4 / 32, np.pi * D * N
energy = lambda load: (load * D / 2) ** 2 * length / (2 * G * J)  # noqa: E731
deflection = (energy(F + 1e-3) - energy(F - 1e-3)) / 2e-3
print(f"k = {spring_rate(d, D, N, G):.1f} N/m; F / deflection = {F / deflection:.1f} N/m")
assert np.isclose(F / deflection, spring_rate(d, D, N, G))

# %%
# The curvature correction
# ------------------------
#
# Wahl's factor exceeds the direct-shear factor :math:`1 + 0.5/C` by the
# curvature effect; Bergsträsser's simpler factor follows it within 2%.

C = spring_index(D, d)
for correction in ("none", "direct", "wahl", "bergstrasser"):
    print(f"{correction:>12}: {spring_shear_stress(F, D, d, correction) / 1e6:.1f} MPa")
c = np.linspace(3, 16, 300)
curvature = wahl_factor(c) / (1 + 0.5 / c)
assert np.all(curvature > 1) and np.all(np.diff(curvature) < 0)
assert np.allclose(bergstrasser_factor(c), wahl_factor(c), rtol=0.02)
print(f"at C = {C:.0f} curvature alone raises the stress by {100 * (wahl_factor(C) / (1 + 0.5 / C) - 1):.1f}%")

fig, ax = plt.subplots()
ax.plot(c, wahl_factor(c), label="Wahl")
ax.plot(c, bergstrasser_factor(c), "--", label="Bergsträsser")
ax.plot(c, 1 + 0.5 / c, ":", label="direct shear only")
ax.set(xlabel="spring index C = D/d", ylabel="stress factor K")
ax.legend()

plt.show()

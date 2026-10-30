r"""
Simmons and Ruge's bonded strain gauge (1938)
=============================================

Stress cannot be measured, only strain. In 1936 Edward Simmons at Caltech,
and in April 1938 Arthur Ruge at MIT, independently cemented fine
resistance wire to a specimen, so that its resistance changed with the
strain of the surface beneath; their names became the "SR-4" gauge. One
gauge measures strain in one direction; a *rosette* of three at known
angles gives the whole plane strain state
(:func:`~engineeringkit.structural_engineering.strain_rosette`), and
Hooke's law turns it into stress
(:func:`~engineeringkit.structural_engineering.plane_stress_from_strains`).
The example simulates a rosette on a shaft under bending and torsion and
recovers the torque from the three readings.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import mohr_circle, plane_stress_from_strains, strain_rosette
from engineeringkit.structural_engineering.visualizers import plot_mohr_circle

# %%
# Readings of a 0-45-90 rosette on a shaft surface
# ------------------------------------------------

E, nu, G = 200e9, 0.3, 200e9 / 2.6
d, M, T = 0.05, 1.2e3, 2.0e3
sx, txy = 32 * M / (np.pi * d**3), 16 * T / (np.pi * d**3)
ex, ey, gxy = sx / E, -nu * sx / E, txy / G
angles = np.radians([0.0, 45.0, 90.0])
readings = ex * np.cos(angles) ** 2 + ey * np.sin(angles) ** 2 + gxy * np.sin(angles) * np.cos(angles)
rng = np.random.default_rng(1938)
noisy = readings + rng.normal(0, 2e-6, 3)  # a few microstrain of noise
print("gauge readings [microstrain]:", np.round(noisy * 1e6, 1))

# %%
# Back to stress and torque
# -------------------------

strains = strain_rosette(noisy, angles)
stress = plane_stress_from_strains(*strains, E, nu)
torque = stress[2] * np.pi * d**3 / 16
print(f"recovered: sigma_x {stress[0] / 1e6:.1f} MPa (true {sx / 1e6:.1f}), tau_xy {stress[2] / 1e6:.1f} MPa (true {txy / 1e6:.1f})")
print(f"torque {torque:.0f} N m (applied {T:.0f} N m)")
exact = plane_stress_from_strains(*strain_rosette(readings, angles), E, nu)
assert np.allclose(exact, (sx, 0.0, txy), atol=1e-3 * sx)
assert abs(torque - T) < 0.05 * T

_, radius, s1, s2 = mohr_circle(*stress)
print(f"principal stresses {s1 / 1e6:.1f} and {s2 / 1e6:.1f} MPa")
ax = plot_mohr_circle(*stress)
ax.set_title("Stress from a strain-gauge rosette")

plt.show()

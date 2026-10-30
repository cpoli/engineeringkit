r"""
Carter's half-wave dipole (1932)
================================

The half-wave dipole is the reference antenna of radio engineering. With
the sinusoidal current a thin wire carries, its far field varies as
:math:`\cos(\tfrac\pi2\cos\theta)/\sin\theta`
(:func:`~engineeringkit.communications.half_wave_dipole_pattern`), a little
narrower than the :math:`\sin\theta` of a short dipole. Integrating the
radiated power gives a directivity of 1.64, or 2.15 dBi
(:func:`~engineeringkit.communications.half_wave_dipole_directivity`), and a
radiation resistance of 73 ohms
(:func:`~engineeringkit.communications.half_wave_dipole_resistance`). Philip
Carter of RCA put such calculations in circuit form in 1932: his
induced-EMF method gave the self- and mutual impedances of dipoles, the
basis of every array design that followed.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

from engineeringkit.communications import (
    FREE_SPACE_IMPEDANCE,
    directivity,
    half_wave_dipole_directivity,
    half_wave_dipole_pattern,
    half_wave_dipole_resistance,
)
from engineeringkit.communications.visualizers import plot_antenna_pattern

# %%
# Directivity and radiation resistance
# ------------------------------------

theta = np.linspace(0, np.pi, 4001)
half_wave = half_wave_dipole_pattern(theta) ** 2
short = np.sin(theta) ** 2
D_half, D_short = directivity(theta, half_wave), directivity(theta, short)
print(f"directivity: short dipole {D_short:.4f}, half-wave {D_half:.4f} (closed form {half_wave_dipole_directivity():.4f})")
assert np.isclose(D_short, 1.5, rtol=1e-5) and np.isclose(D_half, half_wave_dipole_directivity(), rtol=1e-5)

# Radiated power with feed current I0: P = (eta I0^2 / 4 pi) * integral of F^2 sin(theta) d(theta).
integral, _ = integrate.quad(lambda th: np.cos(np.pi / 2 * np.cos(th)) ** 2 / np.sin(th), 1e-9, np.pi - 1e-9)
R = FREE_SPACE_IMPEDANCE / (2 * np.pi) * integral
print(f"radiation resistance {R:.2f} ohm (closed form {half_wave_dipole_resistance():.2f} ohm)")
assert np.isclose(R, half_wave_dipole_resistance(), rtol=1e-6)

# %%
# Beamwidth
# ---------

half_power = theta[(half_wave >= 0.5) & (theta <= np.pi / 2)][0]
beamwidth = np.degrees(np.pi - 2 * half_power)
print(f"half-power beamwidth: half-wave {beamwidth:.1f} deg, short dipole 90 deg")
assert abs(beamwidth - 78) < 0.5

full = np.linspace(0, 2 * np.pi, 721)
fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), subplot_kw={"projection": "polar"})
plot_antenna_pattern(full, half_wave_dipole_pattern(full) ** 2, ax=axes[0], label="half-wave")
plot_antenna_pattern(full, np.sin(full) ** 2, ax=axes[0], label="short dipole")
axes[0].set_title("Power patterns (axis vertical)")
plot_antenna_pattern(full, half_wave_dipole_pattern(full) ** 2, ax=axes[1], db_range=30)
axes[1].set_title("Half-wave dipole [dB]")
fig.tight_layout()

plt.show()

r"""
Newton, Laplace and the speed of sound
======================================

Newton's *Principia* (1687) predicted the speed of sound assuming the
compressions are isothermal, :math:`a = \sqrt{p/\rho} = \sqrt{RT}`, about 15%
below measurement. Pierre-Simon Laplace (1816) argued that the compressions
are too fast to exchange heat, so they are adiabatic and the stiffness rises
by :math:`\gamma`: :math:`a = \sqrt{\gamma R T}`
(:func:`~engineeringkit.aerodynamics.speed_of_sound`). The Mach number
(:func:`~engineeringkit.aerodynamics.mach_number`) measures speed in units of
this local sound speed.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics import R_AIR, isa_atmosphere, mach_number, speed_of_sound

# %%
# Isothermal versus adiabatic
# ---------------------------

T = np.linspace(200, 320, 121)
newton = np.sqrt(R_AIR * T)
laplace = speed_of_sound(1.4, R_AIR, T)
a_measured_20C = 343.2
a_laplace = float(speed_of_sound(1.4, R_AIR, 293.15))
print(f"at 293.15 K: Newton {np.sqrt(R_AIR * 293.15):.1f} m/s, Laplace {a_laplace:.1f} m/s, measured {a_measured_20C} m/s")
assert np.allclose(laplace / newton, np.sqrt(1.4))

# %%
# The same airspeed is a higher Mach number at altitude
# -----------------------------------------------------

h = np.linspace(0, 20000, 81)
state = isa_atmosphere(h)
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
axes[0].plot(T, newton, label="Newton (isothermal)")
axes[0].plot(T, laplace, label="Laplace (adiabatic)")
axes[0].plot(293.15, a_measured_20C, "ko", label="measured, 20 °C")
axes[0].set(xlabel="temperature [K]", ylabel="speed of sound [m/s]")
axes[0].legend()
axes[1].plot(mach_number(250.0, state.speed_of_sound), h / 1e3)
axes[1].set(xlabel="Mach number at 250 m/s", ylabel="altitude [km]", title="ISA")
fig.tight_layout()

plt.show()

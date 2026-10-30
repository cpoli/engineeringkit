r"""
Fick's laws of diffusion (1855)
===============================

Adolf Fick, a 26-year-old prosector in anatomy at Zürich, took Fourier's
law of heat conduction and applied it to salt dissolving in water: the flux
is proportional to the concentration gradient,
:math:`J = -D\,\partial C/\partial x`
(:func:`~engineeringkit.materials.fick_flux`), and conservation gives the
second law :math:`\partial C/\partial t = D\,\partial^2 C/\partial x^2`.
Into a thick steel part whose surface is held at a carbon content
:math:`C_s`, the second law gives the error-function profile
(:func:`~engineeringkit.materials.semi_infinite_concentration`) that times
every carburizing furnace. The example computes the case depth of a gear and
the furnace time it needs.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import diffusion_coefficient, diffusion_time, fick_flux, semi_infinite_concentration

# %%
# Steady flux through a plate
# ---------------------------
#
# Fick's first law: carbon at 1.2 and 0.8 kg/m^3 at depths 5 and 10 mm in an
# iron plate with :math:`D = 3\times10^{-11}` m^2/s.

J = fick_flux(3e-11, (0.8 - 1.2) / (10e-3 - 5e-3))
print(f"flux {J:.2e} kg/(m^2 s)")
assert np.isclose(J, 2.4e-9)

# %%
# Carburizing a gear
# ------------------
#
# Carbon diffuses in austenite with :math:`D_0 = 2.3\times10^{-5}` m^2/s and
# :math:`Q_d = 148` kJ/mol. A 0.25 wt% C steel in a 1.2 wt% atmosphere must
# reach 0.80 wt% at 0.5 mm. Doubling the depth takes four times as long,
# because time enters only through :math:`x/\sqrt{Dt}`.

C0, Cs, target, depth = 0.25, 1.20, 0.80, 0.5e-3
for T in (900.0, 950.0, 1000.0):
    D = diffusion_coefficient(2.3e-5, 148e3, T + 273.15)
    t = diffusion_time(depth, target, D, C0, Cs)
    print(f"{T:.0f} C: D = {D:.2e} m^2/s, {t / 3600:.1f} h")
D = diffusion_coefficient(2.3e-5, 148e3, 950 + 273.15)
t = diffusion_time(depth, target, D, C0, Cs)
assert np.isclose(diffusion_time(2 * depth, target, D, C0, Cs), 4 * t)

x = np.linspace(0, 2e-3, 200)
fig, ax = plt.subplots(figsize=(6, 3.8))
for hours in (1, 4, 16):
    ax.plot(1e3 * x, semi_infinite_concentration(x, hours * 3600, D, C0, Cs), label=f"{hours} h")
ax.plot(1e3 * depth, target, "ko", label=f"0.80% at 0.5 mm, {t / 3600:.1f} h")
ax.set(xlabel="depth [mm]", ylabel="carbon [wt%]", title="Carburizing at 950 C")
ax.legend()
fig.tight_layout()
assert np.isclose(semi_infinite_concentration(depth, t, D, C0, Cs), target)

plt.show()

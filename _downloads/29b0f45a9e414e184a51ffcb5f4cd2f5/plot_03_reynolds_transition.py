r"""
Reynolds's dye experiment (1883)
================================

In 1883 Osborne Reynolds ran water from a tank through glass tubes and fed
a thin streak of dye into the entrance. At low speed the streak stayed a
straight thread; above a critical speed it broke up and coloured the whole
tube. The critical speed fell as the tube widened and as the water warmed,
and Reynolds showed that every case met at one value of
:math:`\rho v D/\mu`
(:func:`~engineeringkit.fluid_engineering.reynolds_number`). This example
computes the critical speed against diameter and temperature, and shows the
jump in friction factor across the transition.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import (
    colebrook_friction_factor,
    laminar_friction_factor,
    reynolds_number,
)

# %%
# Critical speed against bore and temperature
# -------------------------------------------
#
# Water's viscosity falls by more than half between 5 °C and 50 °C (Vogel's
# fit :math:`\mu = 2.414\times10^{-5}\cdot 10^{247.8/(T - 140)}` Pa s). With
# the lower critical Reynolds number 2300, the critical speed is
# :math:`v_c = 2300\,\mu/(\rho D)`.


def water_viscosity(T_celsius):
    return 2.414e-5 * 10 ** (247.8 / (T_celsius + 273.15 - 140))


RE_C, rho = 2300.0, 998.0
D = np.geomspace(5e-3, 0.1, 50)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.8))
for T in [5, 20, 50]:
    mu = water_viscosity(T)
    ax1.loglog(1e3 * D, RE_C * mu / (rho * D), label=f"{T} °C")
ax1.set(xlabel="tube diameter [mm]", ylabel="critical speed [m/s]", title="Where the dye streak breaks up")
ax1.legend()

mu20 = water_viscosity(20)
print(f"μ(20 °C) = {mu20 * 1e3:.3f} mPa s; a 25 mm tube turns at {RE_C * mu20 / (rho * 0.025) * 100:.1f} cm/s")

# %%
# The friction factor jumps
# -------------------------
#
# Below 2300 the friction factor is :math:`64/Re`; above 4000 Colebrook's
# turbulent law applies, and the transition between them is not modelled.
# At the same Reynolds number the turbulent friction is several times the
# laminar value would be.

Re_lam = np.geomspace(500, 2299, 30)
Re_turb = np.geomspace(4000, 1e5, 30)
ax2.loglog(Re_lam, laminar_friction_factor(Re_lam), label="laminar")
ax2.loglog(Re_turb, colebrook_friction_factor(Re_turb, 0.0), label="turbulent, smooth")
ax2.axvspan(2300, 4000, color="0.85")
ax2.set(xlabel="Reynolds number", ylabel="friction factor", title="Laminar and turbulent friction")
ax2.legend()
fig.tight_layout()

ratio = colebrook_friction_factor(4000.0, 0.0) / (64 / 4000)
print(f"at Re = 4000 turbulent friction is {ratio:.1f} times the laminar extrapolation")

# %%
# All critical speeds share one Reynolds number, and the laminar model
# refuses to run beyond it.

for T in [5, 20, 50]:
    mu = water_viscosity(T)
    assert np.allclose(reynolds_number(rho, RE_C * mu / (rho * D), D, mu), RE_C)
try:
    laminar_friction_factor(3000.0)
except ValueError:
    print("laminar_friction_factor(3000) raises ValueError")
else:
    raise AssertionError("expected the laminar model to refuse Re = 3000")
assert ratio > 2

plt.show()

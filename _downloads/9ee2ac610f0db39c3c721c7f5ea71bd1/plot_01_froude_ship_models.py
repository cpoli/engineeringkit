r"""
Froude's ship models (1868)
===========================

In 1868 William Froude proposed to the Admiralty that a ship's resistance
could be predicted from towing a small model, if the model ran at the
"corresponding speed", :math:`v_m = v_s\sqrt{L_m/L_s}`. Then the waves the
two make have the same shape, because the Froude number
:math:`Fr = v/\sqrt{gL}` (:func:`~engineeringkit.fluid_engineering.froude_number`)
is the same. Skin friction does not scale that way, so Froude split the
resistance: friction from flat-plank tests, the rest scaled with the
displacement. His tank at Torquay checked the method against full-scale
towing of HMS *Greyhound* in 1871. This example scales a model test to a
ship, using the ITTC 1957 friction line in place of Froude's planks.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import froude_number, reynolds_number

# %%
# Corresponding speeds
# --------------------
#
# A 120 m ship and a 4 m model (scale 30). Water: model in fresh water,
# ship in sea water.

L_s, L_m = 120.0, 4.0
lam = L_s / L_m
rho_m, mu_m, rho_s, mu_s = 998.0, 1.0e-3, 1025.0, 1.22e-3
v_s = np.linspace(4, 12, 30)
v_m = v_s / np.sqrt(lam)
print(f"scale {lam:g}: a ship at 10 m/s is modelled at {10 / np.sqrt(lam):.2f} m/s, Fr = {froude_number(10, L_s):.3f}")

# %%
# Scaling the resistance
# ----------------------
#
# The total resistance coefficient :math:`C_T = R/(\tfrac12\rho v^2 S)` is
# split into friction :math:`C_F(Re) = 0.075/(\log_{10}Re - 2)^2` and a
# residual (wave) part that depends only on :math:`Fr`. A toy model test
# gives the model's :math:`C_T`, with wave resistance rising steeply near
# :math:`Fr = 0.3`.


def ittc_friction(Re):
    return 0.075 / (np.log10(Re) - 2) ** 2


Fr = froude_number(v_m, L_m)
Re_m = reynolds_number(rho_m, v_m, L_m, mu_m)
C_T_model = ittc_friction(Re_m) + 2e-3 * (Fr / 0.3) ** 6  # the "measured" model test
C_R = C_T_model - ittc_friction(Re_m)  # same at model and full scale
Re_s = reynolds_number(rho_s, v_s, L_s, mu_s)
C_T_ship = ittc_friction(Re_s) + C_R

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(Fr, 1e3 * C_T_model, label="model, $C_T$")
ax.plot(Fr, 1e3 * C_T_ship, label="ship, $C_T$")
ax.plot(Fr, 1e3 * C_R, "--", label="residual (wave), $C_R$")
ax.set(xlabel="Froude number", ylabel="resistance coefficient ×1000", title="Froude's extrapolation")
ax.legend()
fig.tight_layout()

S_s = 3200.0  # wetted surface of the ship [m²]
R_ship = C_T_ship * 0.5 * rho_s * v_s**2 * S_s
print(f"ship resistance at {v_s[-1]:.0f} m/s: {R_ship[-1] / 1e3:.0f} kN")
print(
    f"model Re = {Re_m[-1]:.1e}, ship Re = {Re_s[-1]:.1e}: friction coefficient falls from "
    f"{1e3 * ittc_friction(Re_m[-1]):.2f}e-3 to {1e3 * ittc_friction(Re_s[-1]):.2f}e-3"
)

# %%
# Model and ship share the Froude number at every speed, while their
# Reynolds numbers differ by more than two orders of magnitude.

assert np.allclose(froude_number(v_m, L_m), froude_number(v_s, L_s))
assert np.all(Re_s / Re_m > 100)
assert np.all(C_T_ship < C_T_model)

plt.show()

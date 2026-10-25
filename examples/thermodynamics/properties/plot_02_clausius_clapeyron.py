r"""
The Clausius-Clapeyron relation (1834-1850)
===========================================

Clapeyron (1834) applied Carnot's reasoning to a small cycle of boiling and
condensing and found that the slope of the vapour-pressure curve is fixed
by the latent heat, :math:`\mathrm{d}p/\mathrm{d}T = h_{fg}/[T(v_g - v_f)]`.
Clausius (1850) rederived it from the two laws and, treating the vapour as
an ideal gas and neglecting the liquid volume, integrated it:
:math:`\ln p` is linear in :math:`1/T`
(:func:`~engineeringkit.thermodynamics.clausius_clapeyron_pressure`). The
example checks the exact slope against IAPWS-IF97 steam and uses the
integrated form to find where water boils on a mountain.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import clausius_clapeyron_pressure, saturation_pressure, saturation_temperature, steam_properties

# %%
# Clapeyron's exact slope
# -----------------------
#
# At 100 °C the latent heat and the two specific volumes come from IF97. The
# slope of the IF97 saturation curve, computed independently, must agree.

T = 373.15
p = float(saturation_pressure(T))
liquid, vapour = steam_properties(p, x=0.0), steam_properties(p, x=1.0)
h_fg = vapour.h - liquid.h
slope_clapeyron = h_fg / (T * (vapour.v - liquid.v))
slope_curve = float(saturation_pressure(T + 0.01) - saturation_pressure(T - 0.01)) / 0.02
R_water = 461.526
slope_clausius = h_fg * p / (R_water * T**2)
print(f"h_fg = {h_fg / 1e3:.1f} kJ/kg, v_g = {vapour.v:.4f} m³/kg, v_f = {liquid.v:.6f} m³/kg")
print(f"dp/dT: Clapeyron {slope_clapeyron:.1f} Pa/K, IF97 curve {slope_curve:.1f} Pa/K, ideal vapour {slope_clausius:.1f} Pa/K")
assert np.isclose(slope_clapeyron, slope_curve, rtol=1e-4)
assert abs(slope_clausius / slope_clapeyron - 1) < 0.03  # steam at 1 atm is 2% from ideal

# %%
# A straight line in 1/T
# ----------------------
#
# Anchored at the normal boiling point with constant :math:`h_{fg}`, the
# integrated relation is within 1% of IF97 within 10 K of 100 °C and drifts away as the
# latent heat changes (it falls to zero at the critical point).

T_range = np.linspace(280.0, 600.0, 200)
p_cc = clausius_clapeyron_pressure(T_range, T, p, h_fg, R_water)
p_if97 = saturation_pressure(T_range)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
axes[0].semilogy(1e3 / T_range, p_if97, label="IAPWS-IF97")
axes[0].semilogy(1e3 / T_range, p_cc, "--", label="Clausius-Clapeyron")
axes[0].set(xlabel="1000/T [1/K]", ylabel="saturation pressure [Pa]", title="ln p against 1/T")
axes[0].legend()
axes[1].plot(T_range - 273.15, (p_cc / p_if97 - 1) * 100)
axes[1].axhline(0, color="gray", lw=0.8)
axes[1].set(xlabel="temperature [°C]", ylabel="error [%]", title="Constant latent heat")
fig.tight_layout()
near = np.abs(T_range - T) < 10
assert np.all(np.abs(p_cc[near] / p_if97[near] - 1) < 0.01)

# %%
# Boiling on Everest
# ------------------
#
# Inverting the integrated form,
# :math:`1/T = 1/T_{ref} - (R/h_{fg})\ln(p/p_{ref})`. At the summit the
# pressure is about 33.7 kPa.

p_summit = 33.7e3
T_boil = 1 / (1 / T - R_water / h_fg * np.log(p_summit / p))
print(f"water boils at {T_boil - 273.15:.1f} °C on the summit (IF97: {float(saturation_temperature(p_summit)) - 273.15:.1f} °C)")
assert np.isclose(clausius_clapeyron_pressure(T_boil, T, p, h_fg, R_water), p_summit)
assert abs(T_boil - saturation_temperature(p_summit)) < 1.0

plt.show()

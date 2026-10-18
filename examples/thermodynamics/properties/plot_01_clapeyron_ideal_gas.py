r"""
Clapeyron's ideal gas law (1834)
================================

Boyle (1662) and Mariotte found that at fixed temperature a gas's pressure
times its volume stays constant. Gay-Lussac (1802) found that at fixed
pressure every gas expands by the same fraction of its volume at 0 °C,
about 1/267 per degree. In his 1834 memoir on Carnot's theory, Émile
Clapeyron combined the two as :math:`pv = R(267 + t)`. Regnault's better
measurements moved the 267 to 273, and with absolute temperature and the
amount of gas the law reads :math:`pV = nRT`
(:func:`~engineeringkit.thermodynamics.ideal_gas_pressure`,
:func:`~engineeringkit.thermodynamics.ideal_gas_temperature`). The example
rebuilds the absolute zero from a gas thermometer.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import R_UNIVERSAL, ideal_gas_pressure, ideal_gas_temperature

# %%
# Boyle's isotherms
# -----------------
#
# One mole in a cylinder: halving the volume doubles the pressure, and a
# hotter isotherm lies above a colder one.

n = 1.0
V = np.linspace(0.005, 0.05, 200)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for T in (250.0, 300.0, 400.0):
    p = ideal_gas_pressure(n, T, V)
    np.testing.assert_allclose(p * V, n * R_UNIVERSAL * T)
    axes[0].plot(V * 1e3, p / 1e3, label=f"{T:.0f} K")
axes[0].set(xlabel="volume [L]", ylabel="pressure [kPa]", title="Boyle-Mariotte: pV constant")
axes[0].legend()

# %%
# A constant-volume gas thermometer
# ---------------------------------
#
# A bulb is sealed at the ice point at 80 kPa and heated to the steam point.
# The pressure rises linearly with Celsius temperature. Extrapolating the
# line to zero pressure gives the temperature at which the gas would exert no
# pressure at all: the absolute zero.

V_bulb, p_ice = 1e-3, 80e3
n_bulb = p_ice * V_bulb / (R_UNIVERSAL * 273.15)
t_celsius = np.array([0.0, 25.0, 50.0, 75.0, 100.0])
p_read = ideal_gas_pressure(n_bulb, t_celsius + 273.15, V_bulb)
slope, intercept = np.polyfit(t_celsius, p_read, 1)
t_zero = -intercept / slope
print(f"expansion coefficient {slope / intercept:.6f} per °C = 1/{intercept / slope:.2f}")
print(f"extrapolated zero of pressure: {t_zero:.2f} °C (Gay-Lussac's 1/267 would give -267 °C)")
assert np.isclose(t_zero, -273.15)
assert np.isclose(ideal_gas_temperature(p_read[-1], V_bulb, n_bulb), 373.15)

t_line = np.linspace(-280, 110, 50)
axes[1].plot(t_line, (slope * t_line + intercept) / 1e3, "--", color="gray", label="extrapolation")
axes[1].plot(t_celsius, p_read / 1e3, "o", label="readings")
axes[1].axvline(-273.15, color="C3", lw=0.8)
axes[1].set(xlabel="temperature [°C]", ylabel="bulb pressure [kPa]", title="Absolute zero at -273.15 °C", ylim=(0, None))
axes[1].legend()
fig.tight_layout()

# %%
# The molar volume
# ----------------
#
# At 0 °C and one standard atmosphere one mole of any ideal gas fills
# 22.414 litres, Avogadro's hypothesis in Clapeyron's form.

V_molar = n * R_UNIVERSAL * 273.15 / 101325.0
print(f"molar volume at 0 °C and 1 atm: {V_molar * 1e3:.3f} L")
assert np.isclose(ideal_gas_pressure(1.0, 273.15, V_molar), 101325.0)
assert round(V_molar * 1e3, 3) == 22.414

plt.show()

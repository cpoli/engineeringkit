r"""
The International Standard Atmosphere (1919-1952)
=================================================

Aircraft performance depends on air density, which changes with weather.
Albert Toussaint proposed a reference atmosphere in 1919, and ICAO adopted
the International Standard Atmosphere in 1952: a 6.5 K/km linear lapse to the
tropopause at 11 km, then an isothermal layer. Hydrostatic balance and the
ideal-gas law give closed-form pressure and density
(:func:`~engineeringkit.aerodynamics.isa_atmosphere`). This example checks
hydrostatic balance numerically and shows how thin air raises the stall speed.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics import isa_atmosphere, stall_speed

# %%
# Profiles and the hydrostatic check
# ----------------------------------

h = np.linspace(0, 20000, 401)
state = isa_atmosphere(h)
residual = np.gradient(state.pressure, h) + state.density * 9.80665
relative_residual = np.max(np.abs(residual[1:-1] / (state.density[1:-1] * 9.80665)))
print(f"max |dp/dh + rho g| / (rho g) = {relative_residual:.1e}")
for km in [0, 5, 11, 20]:
    s = isa_atmosphere(km * 1e3)
    print(f"{km:2d} km: T = {float(s.temperature):6.2f} K, p = {float(s.pressure):8.0f} Pa, rho = {float(s.density):.4f} kg/m^3")
assert relative_residual < 1e-3
sea_level, tropopause = isa_atmosphere(0.0), isa_atmosphere(11e3)
assert np.isclose(float(sea_level.pressure), 101325.0) and np.isclose(float(sea_level.density), 1.225)
assert np.isclose(float(tropopause.temperature), 288.15 - 6.5 * 11) and np.isclose(float(tropopause.pressure), 22632, rtol=1e-4)
fig, axes = state.plot()
fig.suptitle("International Standard Atmosphere")
fig.tight_layout()

# %%
# True stall speed rises with altitude
# ------------------------------------

fig, ax = plt.subplots(figsize=(5, 3.6))
vs = stall_speed(10_000.0, state.density, 16.0, 1.6)
ax.plot(vs, h / 1e3)
assert np.all(np.diff(vs) > 0) and np.isclose(vs[-1] / vs[0], np.sqrt(state.density[0] / state.density[-1]))
ax.set(xlabel="stall true airspeed [m/s]", ylabel="altitude [km]", title="W = 10 kN, S = 16 m², $C_{L,max}$ = 1.6")
fig.tight_layout()

plt.show()

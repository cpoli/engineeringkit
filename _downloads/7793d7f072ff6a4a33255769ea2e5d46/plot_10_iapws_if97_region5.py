r"""
IAPWS-IF97 region 5: steam in gas turbines (1997)
=================================================

Above 800 °C steam leaves the range of power-station boilers. It enters
gas-turbine combustors and steam-injected turbines, where water is sprayed
into the hot gas. IF97's region 5 covers 1073.15 K to 2273.15 K, at pressures
up to 50 MPa since the 2007 revision. Like region 2 it is a dimensionless
Gibbs function, of :math:`p/1\,\mathrm{MPa}` and :math:`1000\,\mathrm{K}/T`.
:func:`~engineeringkit.thermodynamics.steam_properties` switches to it
automatically and accepts arrays of states.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.thermodynamics import steam_properties

# %%
# The release's verification values
# ----------------------------------
#
# IF97 Table 42 lists three states; the implementation matches them to nine
# digits.

table = [(1500.0, 0.5e6, 5219.76855, 9.65408875), (1500.0, 30e6, 5167.23514, 7.72970133), (2000.0, 30e6, 6571.22604, 8.53640523)]
states = steam_properties([p for _, p, _, _ in table], [T for T, _, _, _ in table])
for (T, p, h, s), hh, ss in zip(table, states.h, states.s):
    print(f"{T:.0f} K, {p / 1e6:4.1f} MPa: h = {hh / 1e3:.5f} kJ/kg (table {h}), s = {ss / 1e3:.8f} kJ/(kg K) (table {s})")
assert np.all(states.region == 5)
np.testing.assert_allclose(states.h / 1e3, [h for *_, h, _ in table], rtol=1e-8)

# %%
# Steam as a nearly ideal gas
# ---------------------------
#
# Arrays sweep a whole isobar in one call. Across the 1073.15 K boundary
# between regions 2 and 5 the specific heat continues smoothly, and at these
# temperatures steam is close to an ideal gas, :math:`pv/RT \approx 1`.

T = np.linspace(600.0, 2273.15, 200)
fig, ax = plt.subplots(figsize=(6, 3.5))
for p in (0.1e6, 5e6, 30e6):
    st = steam_properties(p, T)
    ax.plot(T - 273.15, st.cp / 1e3, label=f"{p / 1e6:g} MPa")
ax.axvline(800.0, color="gray", ls=":")
ax.set(xlabel="temperature [°C]", ylabel="c_p [kJ/(kg K)]", title="Superheated steam, IF97 regions 2 and 5")
ax.legend()
fig.tight_layout()
z = 0.5e6 * steam_properties(0.5e6, 2000.0).v / (461.526 * 2000.0)
print(f"compressibility of steam at 0.5 MPa and 2000 K: {z:.5f}")
assert abs(z - 1) < 1e-3

plt.show()

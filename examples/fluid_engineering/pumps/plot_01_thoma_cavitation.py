r"""
Thoma's cavitation parameter (1924)
===================================

Water turbines of the 1910s and 1920s ran faster and lower than before,
and their runners were eaten away by cavitation. Dieter Thoma, at Munich,
found in 1924 that machines of one design start to cavitate at the same
ratio :math:`\sigma = \mathrm{NPSH}/H`
(:func:`~engineeringkit.fluid_engineering.thoma_sigma`) of the head margin
above vapour pressure at the runner to the machine head. The same
parameter, through the net positive suction head
(:func:`~engineeringkit.fluid_engineering.npsh_available`), governs how
far a pump may lift water. This example finds the largest suction lift of
a pump at several water temperatures.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import G, cavitation_number, npsh_available, thoma_sigma

# %%
# The suction side
# ----------------
#
# A pump lifts water from an open sump at atmospheric pressure. The
# suction pipe loses 1.2 m of head and enters the pump at 2.5 m/s. At the
# inlet, :math:`p/\rho g + v^2/2g = p_{atm}/\rho g - z_s - h_L`.

p_atm, rho, h_L, v_in = 101325.0, 998.0, 1.2, 2.5
head, sigma_critical = 60.0, 0.08  # pump head and its tested critical sigma


def inlet_pressure(z_s):
    return p_atm - rho * G * (z_s + h_L) - 0.5 * rho * v_in**2


temperatures = {20: 2.34e3, 50: 12.35e3, 80: 47.4e3}  # vapour pressure of water [Pa]
z = np.linspace(0, 8, 81)
fig, ax = plt.subplots(figsize=(6, 4))
lifts = {}
for T, pv in temperatures.items():
    sigma = thoma_sigma(npsh_available(inlet_pressure(z), v_in, pv, rho), head)
    ax.plot(z, sigma, label=f"water at {T} °C")
    lifts[T] = (p_atm - pv) / (rho * G) - h_L - sigma_critical * head
    print(f"{T} °C: largest suction lift {lifts[T]:.2f} m")
ax.axhline(sigma_critical, color="k", ls="--", label="critical σ of the pump")
ax.set(xlabel="suction lift [m]", ylabel="plant σ", title="Thoma's parameter for a 60 m pump")
ax.legend()
fig.tight_layout()

# %%
# Hot water needs the pump below the sump (a negative lift, or flooded
# suction). The plant :math:`\sigma` at the largest lift equals the
# critical value, and the cavitation number at the inlet is the NPSH over
# the inlet velocity head.

z20 = lifts[20]
npsh = npsh_available(inlet_pressure(z20), v_in, temperatures[20], rho)
print(f"cavitation number at the inlet: {cavitation_number(inlet_pressure(z20), temperatures[20], rho, v_in):.2f}")

assert np.isclose(thoma_sigma(npsh, head), sigma_critical)
assert lifts[20] > lifts[50] > 0 > lifts[80]
assert np.isclose(
    cavitation_number(inlet_pressure(z20), temperatures[20], rho, v_in),
    (npsh - v_in**2 / (2 * G)) / (v_in**2 / (2 * G)),
)

plt.show()

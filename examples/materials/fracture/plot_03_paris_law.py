r"""
Paris's crack-growth law
========================

Fracture mechanics tells how long a crack can be; fatigue makes it grow
there. Paul Paris, after measuring cracks in aluminium aircraft panels,
proposed with Gomez and Anderson (1961) and argued with Erdogan (1963) that
the growth per cycle depends on the range of the stress intensity factor
alone, :math:`da/dN = C\,\Delta K^m`
(:func:`~engineeringkit.materials.paris_crack_growth`). Integrating it from
the largest crack an inspection can miss to the critical size gives the
inspection interval of damage-tolerant design
(:func:`~engineeringkit.materials.crack_growth_life`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import crack_growth_life, critical_crack_length, paris_crack_growth, stress_intensity
from engineeringkit.materials.utils import power_law_fit

# %%
# Fitting the law to growth-rate data
# -----------------------------------
#
# Crack-growth tests report :math:`da/dN` against :math:`\Delta K`. On
# log-log axes the middle region is a straight line whose slope is
# :math:`m`. Simulated data for an aluminium alloy, with 20% scatter:

C_true, m_true = 5e-11 * (1e-6) ** 3.5, 3.5  # quoted per (MPa sqrt(m))^3.5, converted to Pa
dK = np.geomspace(5e6, 25e6, 15)
rng = np.random.default_rng(3)
rate = paris_crack_growth(C_true, m_true, dK) * rng.lognormal(0, 0.2, dK.size)
C_fit, m_fit = power_law_fit(dK, rate)
print(f"fitted m = {m_fit:.2f} (true {m_true}), C = {C_fit * 1e6**m_fit:.2e} m/cycle per (MPa sqrt(m))^m")
assert abs(m_fit - m_true) < 0.3

# %%
# Inspection interval of a cracked panel
# --------------------------------------
#
# A 160 mm wide panel with a central crack cycles from 0 to 100 MPa. A
# 2 mm crack (half-length 1 mm) is assumed missed. With the finite-width correction
# :math:`Y = \sqrt{\sec(\pi a/W)}` the crack grows faster near the end than
# with :math:`Y = 1`, and most of the life is spent while the crack is still
# small.

Kc, smax, W, a0 = 35e6, 100e6, 0.160, 1e-3


def secant(a):
    return np.sqrt(1 / np.cos(np.pi * a / W))


ac = critical_crack_length(Kc, smax)
infinite = crack_growth_life(C_true, m_true, smax, a0, ac)
finite = crack_growth_life(C_true, m_true, smax, a0, ac, Y=secant)
closed = 2 * (a0 ** (1 - m_true / 2) - ac ** (1 - m_true / 2)) / ((m_true - 2) * C_true * (smax * np.sqrt(np.pi)) ** m_true)
print(f"critical half-length {1e3 * ac:.1f} mm; life {infinite.life:.0f} cycles (Y = 1), {finite.life:.0f} (secant)")
print(f"cycles spent doubling the initial crack: {100 * infinite.summary()['cycles_to_double'] / infinite.life:.0f}% of the life")
assert np.isclose(infinite.life, closed, rtol=1e-5) and finite.life < infinite.life
assert np.isclose(stress_intensity(1.0, smax, ac), Kc)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.8))
ax1.loglog(dK / 1e6, rate, "o", label="tests")
ax1.loglog(dK / 1e6, C_fit * dK**m_fit, label=f"Paris fit, m = {m_fit:.2f}")
ax1.set(xlabel="Delta K [MPa sqrt(m)]", ylabel="da/dN [m/cycle]", title="Growth-rate data")
ax1.legend()
infinite.plot(ax=ax2)
finite.plot(ax=ax2)
ax2.lines[0].set_label("Y = 1")
ax2.lines[1].set_label("finite width")
ax2.legend()
ax2.set_title("Crack growth from 2 mm")
fig.tight_layout()

plt.show()

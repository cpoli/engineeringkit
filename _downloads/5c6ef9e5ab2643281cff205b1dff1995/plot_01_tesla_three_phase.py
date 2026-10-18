r"""
Tesla's polyphase system (1888)
===============================

In May 1888 Nikola Tesla presented "A new system of alternate current
motors and transformers" to the American Institute of Electrical
Engineers. Currents displaced in phase produce a rotating magnetic field,
which drives an induction motor with no commutator. Polyphase supply has a
second advantage that Mikhail Dolivo-Dobrovolsky exploited with three
phases in 1889-1891: a balanced three-phase load draws *constant*
instantaneous power, :math:`\sqrt3 V_L I_L \cos\phi`
(:func:`~engineeringkit.electrical_engineering.three_phase_power`), and the
three return currents cancel, so three wires do the work of six. This
example checks both, and draws the rotating field.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import real_power, three_phase_power

# %%
# Constant power from three pulsating phases
# ------------------------------------------

V_phase, I_line, phi, f = 230.0, 16.0, np.radians(30), 50.0
w = 2 * np.pi * f
t = np.linspace(0, 2 / f, 1000)
shifts = np.array([0, -2 * np.pi / 3, 2 * np.pi / 3])
v = np.sqrt(2) * V_phase * np.cos(w * t[:, None] + shifts)
i = np.sqrt(2) * I_line * np.cos(w * t[:, None] + shifts - phi)
p_phase = v * i
p_total = p_phase.sum(axis=1)
P = three_phase_power(np.sqrt(3) * V_phase, I_line, np.cos(phi))
print(f"three-phase power {P / 1e3:.3f} kW; instantaneous total varies by {np.ptp(p_total):.1e} W")
assert np.allclose(p_total, P)
assert np.isclose(P, 3 * real_power(V_phase, I_line, np.cos(phi)))
assert np.allclose(i.sum(axis=1), 0, atol=1e-9)  # no neutral current

# %%
# The rotating field
# ------------------
#
# Three coils 120 degrees apart, fed by the three currents, add up to a
# field of constant magnitude turning at the supply frequency.

coil_dirs = np.exp(1j * np.array([0, 2 * np.pi / 3, 4 * np.pi / 3]))
field = (np.cos(w * t[:, None] + shifts) * coil_dirs).sum(axis=1)
assert np.allclose(np.abs(field), 1.5)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
axes[0].plot(t * 1e3, p_phase / 1e3, lw=0.8)
axes[0].plot(t * 1e3, p_total / 1e3, "k", label="sum of the three phases")
axes[0].set(xlabel="time [ms]", ylabel="power [kW]", title="Instantaneous power")
axes[0].legend()
axes[1].plot(field.real, field.imag)
axes[1].set_aspect("equal")
axes[1].set(title="Tip of the rotating field", xlabel="x", ylabel="y")
fig.tight_layout()

plt.show()

r"""
Lodge's syntonic tuning: series resonance
=========================================

Early spark transmitters splattered energy over a wide band, so every
receiver heard every transmitter. Oliver Lodge's 1897 patent on "syntonic"
telegraphy put an inductance in series with the antenna's capacitance at
both ends, tuning transmitter and receiver to the same frequency
:math:`f_0 = 1/2\pi\sqrt{LC}`
(:func:`~engineeringkit.electrical_engineering.resonance_frequency`). At
resonance the inductive and capacitive reactances cancel and the current
is limited only by the resistance; the sharpness of the peak is the quality
factor :math:`Q = \sqrt{L/C}/R`, and the half-power bandwidth is
:math:`f_0/Q`. Marconi bought the patent in 1912. This example measures the
bandwidth of a tuned circuit
(:class:`~engineeringkit.electrical_engineering.SeriesRLC`) and shows how
much tuning rejects a neighbouring station.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import SeriesRLC, resonance_frequency

# %%
# A tuned receiver
# ----------------

L, C = 250e-6, 1e-9
f0 = resonance_frequency(L, C)
print(f"f0 = {f0 / 1e3:.1f} kHz")
f = np.linspace(0.7 * f0, 1.3 * f0, 200_001)
fig, ax = plt.subplots()
for R in (5.0, 20.0, 50.0):
    circuit = SeriesRLC(R, L, C)
    gain = np.abs(1 / circuit.impedance(2 * np.pi * f)) * R  # current relative to its peak
    band = f[gain >= 1 / np.sqrt(2)]
    width = band[-1] - band[0]
    Q = circuit.quality_factor
    print(f"R = {R:4.0f} ohm: Q = {Q:6.1f}, half-power bandwidth {width:7.0f} Hz, f0/Q = {f0 / Q:7.0f} Hz")
    assert np.isclose(width, f0 / Q, rtol=1e-3)
    ax.plot(f / 1e3, gain, label=f"Q = {Q:.0f}")

# %%
# Rejecting a neighbour
# ---------------------
#
# A second station 10% higher in frequency is attenuated by the
# high-:math:`Q` circuit far more than by the low-:math:`Q` one.

neighbour = 1.1 * f0
rejection = [abs(SeriesRLC(R, L, C).impedance(2 * np.pi * neighbour)) / R for R in (5.0, 50.0)]
print(f"neighbour 10% away is {rejection[0]:.1f}x weaker at Q=100, {rejection[1]:.2f}x at Q=10")
detuning = neighbour / f0 - f0 / neighbour
assert np.allclose(rejection, [np.sqrt(1 + (q * detuning) ** 2) for q in (100, 10)])

ax.axvline(neighbour / 1e3, color="gray", ls=":", label="neighbouring station")
ax.set(xlabel="frequency [kHz]", ylabel="relative current", title="Series resonance")
ax.legend()

plt.show()

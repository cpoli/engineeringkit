r"""
Stanley's transformer
=====================

In March 1886 William Stanley, working for George Westinghouse, lit the
shops of Great Barrington, Massachusetts, from a 500 V alternator through
transformers that stepped the voltage down to 100 V at each building. His
design, refined from the Gaulard-Gibbs secondary generator and the
Hungarian ZBD transformer of 1885, had a closed iron core and windings in
parallel across the line. An ideal transformer scales voltage by the turns
ratio, :math:`V_2 = V_1 N_2/N_1`
(:func:`~engineeringkit.electrical_engineering.transformer_secondary`),
current by its inverse and impedance by its square. This example shows why
that made AC distribution win: stepping up the voltage cuts line losses by
the square of the ratio.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import transformer_secondary

# %%
# Voltage, current and impedance
# ------------------------------

N1, N2 = 500, 100
V1 = 500.0
V2 = transformer_secondary(V1, N1, N2)
R_load = 10.0
I2 = V2 / R_load
I1 = I2 * N2 / N1
print(f"V2 = {V2:.0f} V, I2 = {I2:.1f} A, I1 = {I1:.1f} A")
assert np.isclose(V1 * I1, V2 * I2)  # an ideal transformer is lossless
assert np.isclose(V1 / I1, R_load * (N1 / N2) ** 2)  # load reflected by the turns ratio squared

# %%
# Why step up?
# ------------
#
# Deliver 100 kW over a line of 2 ohm total resistance, stepping a 1 kV
# generator up to different voltages.

P, R_line, V_gen = 100e3, 2.0, 1e3
ratios = np.array([1, 2, 5, 10, 20, 50])
V_line = transformer_secondary(V_gen, 1, ratios)
loss = (P / V_line) ** 2 * R_line
for r, vl, lo in zip(ratios, V_line, loss):
    print(f"1:{r:<3d} {vl / 1e3:5.0f} kV line, loss {lo / P:7.3%} of the power")
assert np.allclose(loss[0] / loss, ratios**2)

fig, ax = plt.subplots()
ax.loglog(V_line / 1e3, loss / P * 100, "o-")
ax.set(xlabel="transmission voltage [kV]", ylabel="line loss [% of delivered power]", title="Loss falls as 1/V²")

plt.show()

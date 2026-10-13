"""
Oscillator energy balance
=========================

An unforced undamped oscillator conserves kinetic plus spring energy.
Viscous damping removes that energy; compare both against the initial value.
"""

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import simulate_sdof

t = np.linspace(0, 10, 501)
fig, ax = plt.subplots()
for c in [0, 0.4]:
    result = simulate_sdof(1, c, 4, x0=0.1, t=t)
    E = 0.5 * result.velocity**2 + 2 * result.displacement**2
    if c == 0:
        np.testing.assert_allclose(E, 0.02, rtol=1e-5)
    else:
        assert E[-1] < E[0]
    ax.plot(t, E, label=f"c={c} N s/m")
ax.set(xlabel="time [s]", ylabel="mechanical energy [J]")
ax.legend()

plt.show()

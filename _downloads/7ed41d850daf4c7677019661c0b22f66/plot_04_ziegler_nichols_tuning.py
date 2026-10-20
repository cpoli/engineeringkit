r"""
Ziegler and Nichols tune a PID from the edge of instability
===========================================================

John Ziegler and Nathaniel Nichols of Taylor Instruments (1942) gave plant
engineers a recipe: turn up a proportional gain until the loop oscillates
steadily, note that *ultimate gain* :math:`K_u` and period :math:`T_u`, and
read P, PI or PID settings from a table
(:func:`~engineeringkit.control_systems.ziegler_nichols`). The ultimate
point is exactly the gain margin and phase-crossover frequency of the loop
(:func:`~engineeringkit.control_systems.stability_margins`). For
:math:`1/(s+1)^3`, :math:`K_u = 8` and :math:`\omega_u = \sqrt3` rad/s.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.control_systems import poles, stability_margins, transfer_function, ziegler_nichols

plant = transfer_function([1], [1, 3, 3, 1])
margins = stability_margins(plant)
ku, tu = margins["gain_margin"], 2 * np.pi / margins["phase_cross_frequency"]
print(f"ultimate gain {ku:.4f}, ultimate period {tu:.4f} s")
assert np.isclose(ku, 8.0, rtol=1e-4) and np.isclose(tu, 2 * np.pi / np.sqrt(3), rtol=1e-4)

# %%
# Tuned step responses
# --------------------

t = np.linspace(0, 60, 3001)
fig, ax = plt.subplots()
for kind in ["P", "PI", "PID"]:
    pid = ziegler_nichols(ku, tu, kind)
    assert np.all(poles(pid.closed_loop(plant)).real < 0)
    r = pid.simulate(plant, t)
    ax.plot(r.t, r.y, label=f"{kind}: kp={pid.kp:.2f}, ki={pid.ki:.2f}, kd={pid.kd:.2f}")
    print(kind, {k: round(v, 3) for k, v in r.summary().items()})
    if kind != "P":
        assert abs(r.y[-1] - 1) < 0.01
ax.axhline(1, color="k", ls="--", lw=0.8)
ax.set(xlabel="time [s]", ylabel="output", title=r"Ziegler-Nichols settings for $1/(s+1)^3$")
ax.legend()

plt.show()

r"""
Maxwell's "On Governors": stability from characteristic roots
=============================================================

Governors fitted with integral action to remove droop sometimes "hunted",
oscillating with growing amplitude. James Clerk Maxwell (1868) linearised
the engine-governor equations and showed that motion is stable exactly when
every root of the characteristic polynomial has a negative real part
(:func:`~engineeringkit.control_systems.poles`). For the cubic
:math:`s^3 + a_1 s^2 + a_2 s + a_3` he found the condition
:math:`a_1 a_2 > a_3`: too much integral gain destabilises the governor.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.control_systems import feedback, poles, step_response, transfer_function

a1, a2 = 2.0, 3.0  # engine damping and governor stiffness

# %%
# Sweep the integral gain
# -----------------------
#
# Unity feedback around :math:`k/(s(s^2 + a_1 s + a_2))` gives the cubic
# above with :math:`a_3 = k`; Maxwell predicts instability for
# :math:`k > a_1 a_2 = 6`.

gains = np.linspace(0.5, 12, 231)
growth = np.array([poles(feedback(transfer_function([k], [1, a1, a2, 0]))).real.max() for k in gains])
critical = gains[np.argmax(growth > 0)]
print(f"first unstable gain on the grid: {critical:.2f} (Maxwell: {a1 * a2:.2f})")
assert np.all(growth[gains < a1 * a2 - 0.05] < 0) and np.all(growth[gains > a1 * a2 + 0.05] > 0)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
axes[0].plot(gains, growth)
axes[0].axhline(0, color="k", lw=0.8)
axes[0].axvline(a1 * a2, color="C3", ls="--", label="$a_1 a_2$")
axes[0].set(xlabel="integral gain $k$", ylabel="largest real part of a root [1/s]")
axes[0].legend()

# %%
# Hunting appears past the limit
# ------------------------------

t = np.linspace(0, 30, 1501)
for k in [3.0, 6.0, 7.0]:
    _, y = step_response(feedback(transfer_function([k], [1, a1, a2, 0])), t)
    axes[1].plot(t, y, label=f"$k$ = {k:g}")
axes[1].set(xlabel="time [s]", ylabel="speed response", ylim=(-1, 3))
axes[1].legend()
fig.tight_layout()

plt.show()

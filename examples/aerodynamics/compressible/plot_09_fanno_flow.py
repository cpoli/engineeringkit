r"""
Fanno flow: friction chokes a duct (1904)
=========================================

Gino Fanno's 1904 diploma thesis at the ETH Zürich, written under Aurel
Stodola, followed an adiabatic gas along a constant-area pipe with wall
friction. Friction drives the Mach number towards one, from below in
subsonic flow and from above in supersonic flow, while the stagnation
temperature stays fixed and the stagnation pressure falls. The states lie on
a *Fanno line* in the temperature-entropy plane, whose point of maximum
entropy is the sonic state. :func:`~engineeringkit.aerodynamics.fanno_flow`
gives every ratio to that sonic state, and
:func:`~engineeringkit.aerodynamics.fanno_choking_length` the pipe length
that chokes it. This example marches the friction ODE down a pipe and
checks the closed form.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from engineeringkit.aerodynamics import fanno_choking_length, fanno_flow

gamma, D, f = 1.4, 0.05, 0.02

# %%
# March the Mach number down the pipe
# -----------------------------------
#
# Momentum and energy give
# :math:`\frac{dM}{dx} = \frac{f}{2D}\,\frac{\gamma M^3(1 + \tfrac{\gamma-1}{2}M^2)}{1 - M^2}`
# with the Darcy friction factor :math:`f`.


def dM_dx(x, M):
    return f / (2 * D) * gamma * M**3 * (1 + (gamma - 1) / 2 * M**2) / (1 - M**2)


fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for M_in in (0.3, 0.5, 2.0, 3.0):
    L_star = float(fanno_choking_length(M_in, D, f, gamma))
    sonic = solve_ivp(lambda x, M: dM_dx(x, M), (0, 0.999 * L_star), [M_in], rtol=1e-10, atol=1e-12, dense_output=True)
    x = np.linspace(0, 0.999 * L_star, 200)
    axes[0].plot(x / D, sonic.sol(x)[0], label=f"$M_1$ = {M_in}, $L^*$ = {L_star:.2f} m")
    print(f"M1 = {M_in}: L* = {L_star:.3f} m, Mach after 0.999 L* = {sonic.y[0, -1]:.3f}")
    assert abs(sonic.y[0, -1] - 1) < 0.1
axes[0].axhline(1, color="gray", lw=0.8)
axes[0].set(xlabel="distance along pipe $x/D$", ylabel="Mach number", title="Friction drives M towards 1")
axes[0].legend(fontsize=8)

# %%
# The Fanno line
# --------------
#
# Entropy relative to the sonic state is
# :math:`(s - s^*)/c_p = \ln(T/T^*) - \tfrac{\gamma-1}{\gamma}\ln(p/p^*)`. Its
# maximum, where friction can do no more, is at :math:`M = 1`.

M = np.linspace(0.15, 3.5, 400)
state = fanno_flow(M, gamma)
entropy = np.log(state.temperature_ratio) - (gamma - 1) / gamma * np.log(state.pressure_ratio)
axes[1].plot(entropy, state.temperature_ratio)
axes[1].plot(0, 1, "ko")
axes[1].set(xlabel="$(s - s^*)/c_p$", ylabel="$T/T^*$", title="Fanno line")
axes[1].annotate("subsonic", (entropy[20], state.temperature_ratio[20]))
axes[1].annotate("supersonic", (entropy[-60], state.temperature_ratio[-60]))
fig.tight_layout()

assert np.isclose(M[np.argmax(entropy)], 1.0, atol=0.01)
assert np.isclose(float(fanno_flow(2.0).friction_parameter) / 4, 0.30500, atol=1e-5)

plt.show()

r"""
Breguet's propeller range (c. 1920)
===================================

Breguet's range equation was first written for propeller aircraft, the only
kind there was. Shaft power :math:`P = DV/\eta_{pr}` burns fuel at a rate
set by the engine's power-specific consumption :math:`c_P`. Drag is
:math:`W/(L/D)`, so the speed cancels and
:math:`R = \frac{\eta_{pr}}{g c_P}\frac{L}{D}\ln\frac{W_i}{W_f}`
(:func:`~engineeringkit.aerodynamics.breguet_range_propeller`). A propeller
aircraft goes furthest at the speed of maximum :math:`L/D`, which is slower
than a jet's best-range speed. This example checks the formula against a
time-stepped cruise and plots range against speed for a light aircraft.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics import best_glide_cl, breguet_range_propeller, drag_polar, glide_ratio

g, eta, psfc = 9.80665, 0.8, 7.6e-8
W0, Wf = 13.1e3, 11.7e3
rho, S, Cd0, AR, e = 1.0, 16.2, 0.025, 7.4, 0.8

# %%
# Step the cruise in time at constant L/D
# ---------------------------------------

LD = float(glide_ratio(best_glide_cl(Cd0, AR, e), drag_polar(best_glide_cl(Cd0, AR, e), Cd0, AR, e)))
W, x, dt, V = W0, 0.0, 1.0, 50.0
while W > Wf:
    W -= g * psfc * (W / LD) * V / eta * dt
    x += V * dt
R = float(breguet_range_propeller(eta, psfc, LD, W0, Wf))
print(f"(L/D)max = {LD:.2f}: time-stepped {x / 1e3:.0f} km, Breguet {R / 1e3:.0f} km")
assert np.isclose(x, R, rtol=1e-3)

# %%
# Range against cruise speed
# --------------------------
#
# Holding speed fixed, the lift coefficient and so :math:`L/D` follow from
# :math:`W = \tfrac12\rho V^2 S C_L` at the mean cruise weight. Range peaks
# where :math:`L/D` does.

V = np.linspace(30, 90, 121)
W_mean = np.sqrt(W0 * Wf)
CL = 2 * W_mean / (rho * V**2 * S)
ranges = breguet_range_propeller(eta, psfc, glide_ratio(CL, drag_polar(CL, Cd0, AR, e)), W0, Wf)
V_best = float(np.sqrt(2 * W_mean / (rho * S * best_glide_cl(Cd0, AR, e))))
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(V, ranges / 1e3)
ax.axvline(V_best, color="gray", ls="--", label=f"(L/D)max speed {V_best:.1f} m/s")
ax.set(xlabel="cruise speed [m/s]", ylabel="range [km]", title="Propeller range peaks at maximum L/D")
ax.legend()
fig.tight_layout()
assert abs(V[np.argmax(ranges)] - V_best) < 1.0

plt.show()

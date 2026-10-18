r"""
Cayley's separation of lift and drag in level flight (1809)
===========================================================

George Cayley's *On Aerial Navigation* (1809-1810) proposed the modern
aeroplane by separating the problem: a fixed wing supports the weight (lift),
and a separate propulsion system overcomes the resistance (drag). In steady
level flight lift equals weight and thrust equals drag. Here
:func:`~engineeringkit.aerodynamics.lift` and
:func:`~engineeringkit.aerodynamics.drag` turn coefficients into forces for a
light aircraft, and the required lift coefficient and thrust are traced
across the speed range.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics import drag, drag_polar, dynamic_pressure, lift

# %%
# Lift must equal weight
# ----------------------

weight, S, AR, e, cd0, rho = 10_000.0, 16.0, 7.5, 0.8, 0.027, 1.225
V = np.linspace(28, 80, 200)
q = dynamic_pressure(rho, V)
cl_required = weight / (q * S)
assert np.allclose(lift(q, S, cl_required), weight)

# %%
# Thrust must equal drag
# ----------------------
#
# Parasite drag grows like :math:`V^2`; the induced drag that comes with
# lift falls like :math:`1/V^2`. Their sum has a minimum.

cd = drag_polar(cl_required, cd0, AR, e)
thrust_required = drag(q, S, cd)
parasite = drag(q, S, cd0)
induced = thrust_required - parasite
i = int(np.argmin(thrust_required))
print(f"minimum drag {thrust_required[i]:.0f} N at {V[i]:.1f} m/s (C_L = {cl_required[i]:.2f})")
print(f"there, parasite {parasite[i]:.0f} N ~ induced {induced[i]:.0f} N")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
axes[0].plot(V, cl_required)
axes[0].set(xlabel="speed [m/s]", ylabel="required $C_L$", title="Lift = weight")
axes[1].plot(V, thrust_required, label="total drag = thrust required")
axes[1].plot(V, parasite, "--", label="parasite")
axes[1].plot(V, induced, ":", label="induced")
axes[1].plot(V[i], thrust_required[i], "ko")
axes[1].set(xlabel="speed [m/s]", ylabel="force [N]", title="Thrust = drag")
axes[1].legend()
fig.tight_layout()

plt.show()

r"""
Andronov, Vitt and Khaikin's phase plane (1937)
===============================================

Aleksandr Andronov's school in Gorky studied self-sustained oscillations by
eliminating time and drawing the motion of :math:`\ddot x = f(x, \dot x)` as
a curve in the plane :math:`(x, \dot x)`
(:func:`~engineeringkit.control_systems.phase_portrait`). A relay servo,
:math:`\ddot x = -M\,\mathrm{sign}(x)`, is not linear anywhere, yet its
phase plane is read at a glance. Each half plane holds parabolas
:math:`\dot x^2/2 + M|x| = E`, which join into closed orbits: the servo
oscillates for ever at the amplitude it started with. Feeding back the
rate, :math:`\ddot x = -M\,\mathrm{sign}(x + \tau\dot x)`, tilts the
switching line. Every orbit then loses energy at each switch, reaches the
line and *slides* along it, :math:`\dot x = -x/\tau`, into the origin.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.control_systems import phase_portrait

M = 1.0
t = np.linspace(0, 12, 12001)
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)

# %%
# Pure relay: a family of closed orbits
# -------------------------------------

for x0 in (0.4, 0.7, 1.0):
    orbit = phase_portrait(lambda x, v: -M * np.sign(x), x0, 0.0, t)
    energy = orbit.v**2 / 2 + M * np.abs(orbit.x)
    print(f"x0 = {x0}: energy varies by {np.ptp(energy):.1e} around {M * x0}")
    assert np.allclose(energy, M * x0, rtol=2e-3)
    orbit.plot(ax=axs[0], field=x0 == 1.0)
axs[0].set(title=r"$\ddot x = -M\,\mathrm{sign}(x)$: closed orbits")

# %%
# Rate feedback: switching, then sliding
# --------------------------------------

tau = 0.4
slide = phase_portrait(lambda x, v: -M * np.sign(x + tau * v), 1.0, 0.0, t)
energy = slide.v**2 / 2 + M * np.abs(slide.x)
assert np.all(np.diff(energy) <= 1e-6)  # loses energy at each switch, gains none (to the RK4 step error)

# Once on the switching line the motion obeys x + tau v = 0, so x decays as exp(-t / tau).
on_line = np.abs(slide.x + tau * slide.v) < 2e-3
start = np.flatnonzero(~on_line)[-1] + 1
later = (t > t[start]) & (t < t[start] + 2)
print(f"sliding from t = {t[start]:.2f} s at x = {slide.x[start]:.4f}")
assert np.allclose(slide.x[later], slide.x[start] * np.exp(-(t[later] - t[start]) / tau), atol=1e-4)
slide.plot(ax=axs[1])
line = np.linspace(-1.2, 1.2, 2)
axs[1].plot(line, -line / tau, "r--", lw=1, label=r"switching line $x + \tau\dot x = 0$")
axs[1].legend(loc="lower left", fontsize=8)
axs[1].set(title=rf"$\ddot x = -M\,\mathrm{{sign}}(x + {tau}\,\dot x)$: sliding to rest")
print("final state with rate feedback:", slide.summary())
assert abs(slide.x[-1]) < 1e-6 and abs(slide.v[-1]) < 1e-6
fig.tight_layout()
plt.show()

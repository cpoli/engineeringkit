r"""
Reuleaux's kinematics of machinery: the slider-crank (1875)
===========================================================

Franz Reuleaux's *Theoretische Kinematik* (1875) treated every machine as a
*kinematic chain* of rigid links joined by pairs (turning, sliding, screw),
with one link held fixed. Fixing different links of the same chain gives
different machines, which he called inversions. Grounding the cylinder of
the slider-crank chain gives the piston engine, grounding the crank gives
the Whitworth quick-return, and grounding the rod gives the oscillating
cylinder. For the engine, the piston position
:math:`x = r\cos\theta + \sqrt{l^2 - r^2\sin^2\theta}`
(:func:`~engineeringkit.mechanical_engineering.slider_crank_displacement`)
is not harmonic. Its acceleration at top dead centre is
:math:`r\omega^2(1 + r/l)` but only :math:`r\omega^2(1 - r/l)` at bottom
dead centre, which is the source of an engine's second-order shaking force.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import centripetal_force, slider_crank_displacement
from engineeringkit.mechanical_engineering.visualizers import plot_linkage

# %%
# Piston motion by differentiation
# --------------------------------
#
# Differentiate the exact position numerically, with central differences
# round one full turn, at a steady 3000 rpm.

r, l, omega = 0.045, 0.150, 3000 * 2 * np.pi / 60
theta = np.linspace(0, 2 * np.pi, 7200, endpoint=False)
dt = (theta[1] - theta[0]) / omega
x = slider_crank_displacement(r, l, theta)
a = (np.roll(x, -1) - 2 * x + np.roll(x, 1)) / dt**2
print(f"stroke {1e3 * (x.max() - x.min()):.1f} mm (2r = {2e3 * r:.1f} mm)")
print(f"acceleration at TDC {a[0]:.0f} m/s^2, exact {-r * omega**2 * (1 + r / l):.0f}")
print(f"acceleration at BDC {a[3600]:.0f} m/s^2, exact {r * omega**2 * (1 - r / l):.0f}")
assert np.isclose(x.max() - x.min(), 2 * r)
assert np.isclose(a[0], -r * omega**2 * (1 + r / l), rtol=1e-4)
assert np.isclose(a[3600], r * omega**2 * (1 - r / l), rtol=1e-4)

# %%
# Primary and secondary shaking
# -----------------------------
#
# Expanding the root in powers of :math:`r/l` gives
# :math:`\ddot x \approx -r\omega^2(\cos\theta + (r/l)\cos 2\theta)`. A
# counterweight cancels the first harmonic, but nothing on the crank can
# cancel the second. The crank pin's own centripetal load for a 0.4 kg
# rotating mass is shown for scale.

approximate = -r * omega**2 * (np.cos(theta) + r / l * np.cos(2 * theta))
print(f"largest error of the two-term series: {np.max(np.abs(a - approximate)) / np.max(np.abs(a)):.2%} of the peak")
print(f"rotating crank-pin load: {centripetal_force(0.4, r * omega, r) / 1e3:.2f} kN")
assert np.max(np.abs(a - approximate)) < 0.02 * np.max(np.abs(a))

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
for angle in np.radians([0, 60, 120, 180, 240, 300]):
    pin = r * np.array([np.cos(angle), np.sin(angle)])
    piston = [slider_crank_displacement(r, l, angle), 0.0]
    plot_linkage([[0, 0], pin, piston], ax=axes[0], color="C0", alpha=0.5)
axes[0].set(title="Crank, rod and piston at six angles")
axes[1].plot(np.degrees(theta), a / 1e3, label="exact")
axes[1].plot(np.degrees(theta), approximate / 1e3, "--", label="first two harmonics")
axes[1].set(xlabel="crank angle [deg]", ylabel="piston acceleration [km/s^2]", title="Not simple harmonic")
axes[1].legend()
fig.tight_layout()

plt.show()

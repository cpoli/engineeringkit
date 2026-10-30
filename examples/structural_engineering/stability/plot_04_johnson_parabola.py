r"""
Johnson's parabola (1893)
=========================

Euler's hyperbola overestimates the strength of stocky columns, which
yield or buckle inelastically, and tests scattered below it. J. B. Johnson
and his co-authors, in *The Theory and Practice of Modern Framed
Structures* (1893), fitted a parabola that starts at the yield stress and
meets Euler's curve tangentially at half of it
(:func:`~engineeringkit.structural_engineering.johnson_critical_stress`),
at the transition slenderness :math:`C_c = \sqrt{2\pi^2E/\sigma_y}`. The
AISC column formula used the same parabola for decades. End conditions
enter through the effective length :math:`KL`
(:func:`~engineeringkit.structural_engineering.effective_length_factor`).
The example draws the design curve and checks that the two branches meet
with the same slope.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import effective_length_factor, johnson_critical_stress

# %%
# The design curve
# ----------------

E, sy = 200e9, 250e6
Cc = np.sqrt(2 * np.pi**2 * E / sy)
slenderness = np.linspace(1, 250, 500)
sigma = johnson_critical_stress(sy, E, slenderness)
print(f"transition slenderness Cc = {Cc:.1f}, where the critical stress is {johnson_critical_stress(sy, E, Cc) / 1e6:.0f} MPa")
h = 1e-4
left = (johnson_critical_stress(sy, E, Cc - h) - johnson_critical_stress(sy, E, Cc - 2 * h)) / h
right = (johnson_critical_stress(sy, E, Cc + 2 * h) - johnson_critical_stress(sy, E, Cc + h)) / h
assert np.isclose(left, right, rtol=1e-4)

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(slenderness, sigma / 1e6, label="Johnson + Euler")
ax.plot(slenderness, np.pi**2 * E / slenderness**2 / 1e6, "--", label="Euler alone")
ax.axvline(Cc, color="gray", lw=0.5)
ax.set(xlabel="slenderness KL/r", ylabel="critical stress [MPa]", ylim=(0, 300))
ax.legend()
fig.tight_layout()

# %%
# The same strut with four end conditions
# ---------------------------------------

L, r = 3.0, 0.03
for ends in ("fixed-fixed", "fixed-pinned", "pinned-pinned", "fixed-free"):
    K = effective_length_factor(ends)
    print(f"{ends:14s} K = {K:.3f}  KL/r = {K * L / r:6.1f}  sigma_cr = {johnson_critical_stress(sy, E, K * L / r) / 1e6:6.1f} MPa")
assert np.isclose(np.tan(np.pi / effective_length_factor("fixed-pinned")), np.pi / effective_length_factor("fixed-pinned"))

plt.show()

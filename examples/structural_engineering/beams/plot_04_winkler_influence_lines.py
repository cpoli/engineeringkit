r"""
Winkler's influence lines (1867)
================================

A railway bridge must be checked for a train in every position. Emil
Winkler (1867-1868), and Otto Mohr independently (1868), plotted the value
of a response, such as the moment at one section, against the position of a
unit load: the *influence line*
(:func:`~engineeringkit.structural_engineering.influence_line`). The
response to any train is then the sum of the axle loads times the
ordinates under them, and the worst position can be read off the plot. The
example finds the worst position of a four-axle locomotive on a 20 m span
and checks it against a brute-force superposition of point-load solutions.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import SimplySupportedBeam, influence_line

# %%
# Influence lines for a section at the quarter span
# -------------------------------------------------

L, section = 20.0, 5.0
xi = np.linspace(0, L, 401)
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(6, 6))
for ax, quantity, unit in zip(axes, ["left_reaction", "shear", "moment"], ["[-]", "[-]", "[m]"]):
    ax.plot(xi, influence_line(L, xi, quantity, at=section))
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_ylabel(f"{quantity} {unit}")
axes[-1].set_xlabel("unit-load position [m]")
fig.tight_layout()

# %%
# Moving a locomotive across the span
# -----------------------------------

axles = np.array([0.0, 1.8, 3.6, 5.4])  # axle offsets from the front axle [m]
weights = np.array([200e3, 200e3, 200e3, 200e3])  # [N]
fronts = np.linspace(0, L + axles[-1], 1001)


def section_moment(front):
    pos = front - axles
    on = (pos >= 0) & (pos <= L)
    return float(np.sum(weights[on] * influence_line(L, pos[on], "moment", at=section)))


by_influence = np.array([section_moment(f) for f in fronts])
best = fronts[np.argmax(by_influence)]
print(f"worst front-axle position {best:.2f} m, moment at x = {section} m: {by_influence.max() / 1e3:.0f} kN m")

# %%
# Check against superposed point-load solutions of the beam itself.

beam = SimplySupportedBeam(L, 200e9, 0.05)
n = 801
k = int(round(section / L * (n - 1)))
direct = sum(beam.solve_point_load(Wt, position=best - a, n=n).moment[k] for a, Wt in zip(axles, weights) if 0 <= best - a <= L)
print(f"by superposition: {direct / 1e3:.0f} kN m")
assert np.isclose(direct, by_influence.max())

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(fronts, by_influence / 1e3)
ax.axvline(best, color="C3", ls="--", label="worst position")
ax.set(xlabel="front-axle position [m]", ylabel="moment at quarter span [kN m]", title="Four-axle locomotive")
ax.legend()
fig.tight_layout()

plt.show()

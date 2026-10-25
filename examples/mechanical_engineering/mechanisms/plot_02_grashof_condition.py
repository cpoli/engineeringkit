r"""
Grashof's condition for four-bar linkages (1883)
================================================

Franz Grashof (1883) gave the rule that decides whether any link of a
four-bar chain can turn a full circle. With :math:`s` and :math:`l` the
shortest and longest links and :math:`p, q` the others, the shortest link
rotates fully relative to its neighbours if and only if
:math:`s + l \le p + q`
(:func:`~engineeringkit.mechanical_engineering.fourbar_grashof`). Which link
is grounded then decides the machine. Grounding a neighbour of the shortest
link gives a crank-rocker, grounding the shortest gives a double crank, and
grounding the link opposite it gives a double rocker. The example checks
all four inversions of a Grashof chain, and of a non-Grashof chain, by
trying to assemble each one at every input angle
(:func:`~engineeringkit.mechanical_engineering.fourbar_output_angle`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import fourbar_grashof, fourbar_output_angle

angles = np.linspace(0, 2 * np.pi, 721)


def rotates_fully(ground, crank, coupler, rocker):
    try:
        fourbar_output_angle(crank, coupler, rocker, ground, angles)
    except ValueError:
        return False
    return True


def classify(chain, g):
    """Ground link ``g`` of the cyclic chain; test each neighbour as a full-turn input."""
    d, a, b, c = (chain[(g + i) % 4] for i in range(4))
    first, second = rotates_fully(d, a, b, c), rotates_fully(d, c, b, a)
    return {2: "double crank", 1: "crank-rocker", 0: "double rocker"}[first + second]


# %%
# Four inversions of each chain
# -----------------------------

expected = {"shortest grounded": "double crank", "neighbour grounded": "crank-rocker", "opposite grounded": "double rocker"}
for chain in ([2.0, 7.0, 9.0, 6.0], [3.0, 5.0, 9.0, 4.0]):
    grashof = fourbar_grashof(chain)
    print(f"chain {chain}: s + l = {min(chain) + max(chain)}, p + q = {sum(chain) - min(chain) - max(chain)}, Grashof: {grashof}")
    shortest = int(np.argmin(chain))
    for g in range(4):
        kind = classify(chain, g)
        role = "shortest" if g == shortest else "opposite" if (g - shortest) % 4 == 2 else "neighbour"
        print(f"  ground link {chain[g]:.0f} ({role}): {kind}")
        assert kind == (expected[f"{role} grounded"] if grashof else "double rocker")

# %%
# Crank-rocker motion
# -------------------
#
# Ground a 7-unit link next to the 2-unit crank, with a 9-unit coupler and
# a 6-unit rocker: the crank turns fully and the rocker swings back and
# forth.

d, a, b, c = 7.0, 2.0, 9.0, 6.0
theta4 = fourbar_output_angle(a, b, c, d, angles)
print(f"rocker swing {np.degrees(theta4.max() - theta4.min()):.1f} deg")

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
axes[0].plot(np.degrees(angles), np.degrees(theta4))
axes[0].set(xlabel="crank angle [deg]", ylabel="rocker angle [deg]", title="Crank-rocker")
for angle in np.radians([0, 90, 180, 270]):
    A = a * np.array([np.cos(angle), np.sin(angle)])
    t4 = fourbar_output_angle(a, b, c, d, angle)
    B = np.array([d + c * np.cos(t4), c * np.sin(t4)])
    axes[1].plot([0, A[0], B[0], d], [0, A[1], B[1], 0], "o-", alpha=0.7)
axes[1].plot([0, d], [0, 0], "k", lw=3)
axes[1].set_aspect("equal")
axes[1].set(title="Positions at 0, 90, 180, 270 deg")
fig.tight_layout()

plt.show()

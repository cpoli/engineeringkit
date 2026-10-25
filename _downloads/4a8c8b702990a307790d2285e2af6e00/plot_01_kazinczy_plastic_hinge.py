r"""
Kazinczy's plastic hinge (1914)
===============================

In 1914 Gábor Kazinczy loaded steel beams built into concrete walls in
Budapest. They carried far more than the elastic theory allowed: when the
most stressed section yielded through its whole depth it went on turning at
a constant *plastic moment* :math:`M_p = \sigma_yZ_p` (a plastic hinge),
and the load moved to other sections until enough hinges had formed to
make a mechanism. For a rectangle :math:`Z_p = bh^2/4`, 1.5 times the
elastic modulus (:func:`~engineeringkit.structural_engineering.plastic_section_modulus`).
The collapse load follows from the work equation of the mechanism
(:func:`~engineeringkit.structural_engineering.plastic_collapse_load`). The
example follows a built-in beam under a uniform load from first yield,
through hinges at the walls, to collapse with a third hinge at midspan.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import (
    PlanarFrame,
    plastic_collapse_load,
    plastic_section_modulus,
    rectangular_second_moment,
    rectangular_section_modulus,
)

# %%
# Section and beam
# ----------------

b, h, sigma_y, L, E = 0.1, 0.3, 250e6, 8.0, 200e9
Z, Zp = rectangular_section_modulus(b, h), plastic_section_modulus(b, h)
My, Mp = sigma_y * Z, sigma_y * Zp
print(f"first-yield moment {My / 1e3:.1f} kN m, plastic moment {Mp / 1e3:.1f} kN m, shape factor {Zp / Z:.2f}")

# %%
# Stage 1: elastic, both ends built in
# ------------------------------------
#
# The stiffness solution gives the wall moment :math:`wL^2/12`, twice the
# midspan moment. The walls yield first and form hinges at
# :math:`w_1 = 12M_p/L^2`.

beam = PlanarFrame([[0, 0], [L, 0]], ((0, 1),), E, b * h, rectangular_second_moment(b, h))
unit = beam.solve(np.zeros(6), [0, 1, 2, 3, 4, 5], member_loads=1.0)
wall_per_w = unit.member_end_forces[0, 2]  # counterclockwise end moment per unit load
w_first_yield, w1 = My / wall_per_w, Mp / wall_per_w
print(f"first yield at w = {w_first_yield / 1e3:.2f} kN/m, wall hinges at w1 = {w1 / 1e3:.2f} kN/m")
assert np.isclose(w1, 12 * Mp / L**2)

# %%
# Stage 2: hinges at the walls
# ----------------------------
#
# With the wall moments held at :math:`-M_p`, the beam behaves as simply
# supported with end moments. The midspan moment is
# :math:`wL^2/8 - M_p`, which reaches :math:`M_p` at the collapse load.

w_c = 2 * Mp * 8 / L**2
assert np.isclose(w_c, plastic_collapse_load(Mp, L, "fixed", "uniform"))
print(f"collapse at w_c = {w_c / 1e3:.2f} kN/m: {w_c / w_first_yield:.2f} times first yield")
print(f"a simply supported beam would collapse at {plastic_collapse_load(Mp, L, 'simple', 'uniform') / 1e3:.2f} kN/m")

x = np.linspace(0, L, 201)
fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))
for w, label in [(w_first_yield, "first yield"), (w1, "wall hinges"), (w_c, "collapse")]:
    wall = min(w * L**2 / 12, Mp)
    axes[0].plot(x, (w * x * (L - x) / 2 - wall) / 1e3, label=f"{label}, w = {w / 1e3:.0f} kN/m")
for level in (Mp, -Mp):
    axes[0].axhline(level / 1e3, color="gray", ls=":")
axes[0].set(xlabel="x [m]", ylabel="bending moment [kN m]", title="Moment redistribution")
axes[0].legend(fontsize=8)

# Load-deflection: the midspan deflection grows faster once the walls hinge.
w_path = np.linspace(0, w_c, 200)
I = rectangular_second_moment(b, h)
elastic = w_path * L**4 / (384 * E * I)
hinged = w1 * L**4 / (384 * E * I) + (w_path - w1) * 5 * L**4 / (384 * E * I)
axes[1].plot(1e3 * np.where(w_path <= w1, elastic, hinged), w_path / 1e3)
axes[1].axhline(w_c / 1e3, color="C3", ls="--", label="collapse load")
axes[1].set(xlabel="midspan deflection [mm]", ylabel="load [kN/m]", title="Load path")
axes[1].legend()
fig.tight_layout()

plt.show()

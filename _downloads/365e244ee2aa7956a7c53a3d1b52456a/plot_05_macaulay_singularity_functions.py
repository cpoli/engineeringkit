r"""
Macaulay's singularity functions (1919)
=======================================

Integrating :math:`EIy'' = M` for a beam with several loads used to mean a
separate moment equation, and two constants of integration, for every
stretch between loads, all joined by continuity conditions. William Herrick
Macaulay (1919), following Clebsch (1862), wrote the moment once for the
whole beam with brackets :math:`\langle x - a\rangle^n` that vanish to the
left of :math:`a`, so only two constants remain.
:class:`~engineeringkit.structural_engineering.Beam` adds the unknown
reactions to those two constants and solves equilibrium and the support
conditions together, so statically indeterminate beams need no special
treatment. The example solves a propped cantilever under a point load, a
couple and a triangular load, and checks it against the older method of
superposition (:func:`~engineeringkit.structural_engineering.superpose_beam_results`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import Beam, superpose_beam_results

# %%
# A propped cantilever under three kinds of load
# ----------------------------------------------

L, E, I = 8.0, 200e9, 1.2e-4
loads = {"point_loads": [(3.0, 40e3)], "moments": [(6.0, 30e3)], "distributed_loads": [(4.0, 8.0, 0.0, 15e3)]}
propped = Beam(L, E, I, ((0.0, "fixed"), (L, "roller")))
result = propped.solve(**loads, n=801)
(RA, RB), (MA,) = result.reactions, result.support_moments
print(f"reactions {RA / 1e3:.2f} kN and {RB / 1e3:.2f} kN, wall moment {MA / 1e3:.2f} kN m")
print(f"largest deflection {result.deflection.min() * 1e3:.2f} mm at x = {result.x[np.argmin(result.deflection)]:.2f} m")
assert np.isclose(RA + RB, 40e3 + 15e3 * 4 / 2)
assert abs(result.deflection[0]) < 1e-15 and abs(result.slope[0]) < 1e-15 and abs(result.deflection[-1]) < 1e-15

# %%
# The same beam by superposition
# ------------------------------
#
# Remove the roller: the cantilever under the loads deflects at its tip,
# and the redundant reaction is the tip force that pushes it back.

cantilever = Beam(L, E, I, ((0.0, "fixed"),))
free = cantilever.solve(**loads, n=801)
unit = cantilever.solve(point_loads=[(L, -1.0)], n=801)  # an upward unit force at the tip
R = -free.deflection[-1] / unit.deflection[-1]
combined = superpose_beam_results(free, cantilever.solve(point_loads=[(L, -R)], n=801))
print(f"redundant by superposition: {R / 1e3:.2f} kN")
assert np.isclose(R, RB)
assert np.allclose(combined.deflection, result.deflection, atol=1e-12)

fig, axes = result.plot()
axes[0].set_ylabel("shear [N]")
axes[1].set_ylabel("moment [N m]")
axes[2].set_ylabel("deflection [m]")
axes[0].set_title("Propped cantilever by singularity functions")
fig.tight_layout()

plt.show()

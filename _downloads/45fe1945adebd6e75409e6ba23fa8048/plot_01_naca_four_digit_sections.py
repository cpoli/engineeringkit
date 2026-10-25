r"""
The NACA four-digit airfoil family (1933)
=========================================

In NACA Report 460 (1933) Eastman Jacobs and colleagues tested 78 related
sections in the variable-density tunnel. Each is defined by a single
thickness polynomial wrapped around a two-parabola camber line, so four digits
specify the whole shape: ``2412`` has 2% camber at 40% chord and is 12% thick.
:func:`~engineeringkit.aerodynamics.naca4_airfoil` builds the coordinates;
:meth:`~engineeringkit.aerodynamics.Airfoil.summary` recovers the digits from
the geometry.
"""

# %%
import matplotlib.pyplot as plt

from engineeringkit.aerodynamics import naca4_airfoil
from engineeringkit.aerodynamics.visualizers import plot_airfoil

# %%
# Four members of the family
# --------------------------

fig, axes = plt.subplots(4, 1, figsize=(7, 7.5), sharex=True)
for ax, designation in zip(axes, ["0012", "2412", "4412", "6409"], strict=True):
    foil = naca4_airfoil(designation, n=201)
    plot_airfoil(foil, ax=ax)
    s = foil.summary()
    print(
        f"{s['name']}: max thickness {s['max_thickness']:.3f} at x/c = {s['max_thickness_at']:.2f}, "
        f"max camber {s['max_camber']:.3f} at x/c = {s['max_camber_at']:.2f}"
    )
    assert abs(s["max_thickness"] - int(designation[2:]) / 100) < 1e-3
    ax.set_xlabel("")
axes[-1].set_xlabel("x / c")
fig.tight_layout()

plt.show()

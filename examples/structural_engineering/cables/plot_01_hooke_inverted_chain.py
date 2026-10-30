r"""
Hooke's inverted chain (1675)
=============================

Robert Hooke hid his solution to the shape of an arch in a Latin anagram
at the end of *A Description of Helioscopes* (1675); decoded after his
death, it read *ut pendet continuum flexile, sic stabit contiguum rigidum
inversum*: as hangs the flexible line, so but inverted will stand the rigid
arch. A hanging cable can only pull, so it takes the shape of its loads'
moment diagram; turned upside down, the same shape is an arch that only
pushes. Christopher Wren used the idea for the dome of St Paul's. For a
three-hinged arch (:func:`~engineeringkit.structural_engineering.three_hinged_arch`)
the moment is the simple-beam moment less the thrust times the rise,
:math:`M = M_0 - Hy`, and on a parabolic arch under a load uniform along
the span it vanishes everywhere. The example checks this, then moves the
load off the crown to show the bending that an arch of the wrong shape must
carry.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import parabolic_cable, three_hinged_arch

# %%
# The cable and the arch carry the same load with the same force
# --------------------------------------------------------------

L, f, w = 40.0, 8.0, 20e3
arch = three_hinged_arch(L, f, w=w)
H_cable, T_cable, _ = parabolic_cable(w, L, f)
print(f"arch thrust {arch.thrust / 1e3:.0f} kN, cable tension {H_cable / 1e3:.0f} kN")
print(f"largest arch moment {np.max(np.abs(arch.moment)):.1e} N m")
print(f"springing compression {-arch.axial_force[0] / 1e3:.0f} kN, cable support tension {T_cable / 1e3:.0f} kN")
assert np.isclose(arch.thrust, H_cable)
assert np.isclose(-arch.axial_force[0], T_cable)
assert np.max(np.abs(arch.moment)) < 1e-6 * w * L**2

# %%
# A point load is not funicular for a parabola
# --------------------------------------------

P = 400e3
patterns = {
    "uniform load": {"w": w},
    "point load at quarter span": {"point_loads": [(L / 4, P)]},
    "both": {"w": w, "point_loads": [(L / 4, P)]},
}
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7, 5))
axes[0].plot(arch.x, arch.y, color="k")
axes[0].plot([0, L / 2, L], [0, f, 0], "o", mfc="white", color="k")
axes[0].plot(arch.x, f - 4 * f * arch.x * (L - arch.x) / L**2 - f, ":", color="gray", label="hanging cable, inverted")
axes[0].set(ylabel="y [m]", title="Three-hinged parabolic arch")
axes[0].set_aspect("equal")
axes[0].legend(loc="lower right")
for label, kwargs in patterns.items():
    r = three_hinged_arch(L, f, **kwargs)
    axes[1].plot(r.x, r.moment / 1e3, label=label)
    assert np.allclose(r.moment[[0, -1]], 0) and np.isclose(np.interp(L / 2, r.x, r.moment), 0, atol=1e-6)
axes[1].axhline(0, color="gray", lw=0.5)
axes[1].set(xlabel="x [m]", ylabel="moment [kN m]")
axes[1].legend()
fig.tight_layout()

plt.show()

r"""
Nikuradse's sand-roughened pipes (1933)
=======================================

Johann Nikuradse, in Prandtl's laboratory at Göttingen, glued graded sand
to the inside of brass pipes and measured friction over six roughnesses,
from a radius 15 to 507 times the grain size, published in 1933. The
curves share the smooth-pipe law at low Reynolds number
(:func:`~engineeringkit.fluid_engineering.smooth_pipe_friction_factor`),
leave it once the grains poke through the viscous sublayer, and settle on a
constant that depends only on the roughness,
:math:`1/\sqrt f = 1.74 + 2\log_{10}(R/k_s)`
(:func:`~engineeringkit.fluid_engineering.nikuradse_friction_factor`). This
example draws the two limits and finds where each roughness becomes fully
rough.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.fluid_engineering import nikuradse_friction_factor, smooth_pipe_friction_factor

# %%
# Smooth and fully rough limits
# -----------------------------
#
# Nikuradse's radius-to-grain ratios :math:`R/k_s` were 15, 30.6, 60, 126,
# 252 and 507, so :math:`k_s/D = 1/(2R/k_s)`.

ratios = np.array([15, 30.6, 60, 126, 252, 507])
rough = nikuradse_friction_factor(1 / (2 * ratios))
Re = np.geomspace(4e3, 1e7, 200)
smooth = smooth_pipe_friction_factor(Re)

fig, ax = plt.subplots(figsize=(6.5, 4.2))
ax.loglog(Re, smooth, "k", label="smooth (Prandtl-Kármán)")
for ratio, f_r in zip(ratios, rough):
    # Fully rough once the roughness Reynolds number Re (k/D) sqrt(f/8) exceeds 70.
    Re_rough = 70 / ((1 / (2 * ratio)) * np.sqrt(f_r / 8))
    ax.loglog([Re_rough, 1e7], [f_r, f_r], label=f"R/k = {ratio:g}")
    ax.plot(Re_rough, f_r, "o", color=ax.lines[-1].get_color())
    print(f"R/k = {ratio:5g}: f = {f_r:.4f}, fully rough above Re = {Re_rough:.2e}")
ax.set(xlabel="Reynolds number", ylabel="friction factor", title="Nikuradse's limits")
ax.legend(fontsize="small")
fig.tight_layout()

# %%
# The rough limit does not depend on the Reynolds number at all; a pipe
# twice as rough relative to its size has about a third more friction. The
# smooth law falls without limit as the Reynolds number grows, so a pipe of
# any roughness eventually becomes rough.

crossing = Re[np.argmax(smooth < rough[-1])]
print(f"the smoothest pipe's rough limit meets the smooth law near Re = {crossing:.1e}")
assert np.all(np.diff(rough) < 0)
assert np.all(np.diff(smooth) < 0)
assert np.isclose(1 / np.sqrt(rough[0]), 1.74 + 2 * np.log10(15))

plt.show()

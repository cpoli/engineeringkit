r"""
Nusselt's unit-sphere view factor (1928)
========================================

How much of the radiation leaving one diffuse surface lands on another? The
view factor is a fourfold integral of :math:`\cos\theta_1\cos\theta_2/\pi r^2`.
In 1928 Wilhelm Nusselt found a construction that a draughtsman could do:
project the second surface onto a unit hemisphere centred on a point of the
first, then project that image straight down onto the base circle. The
fraction of the circle covered is the view factor from the point. Picking
points uniformly in the circle is the same as picking cosine-weighted
directions, which makes the construction a Monte Carlo method. This example
uses it to check the closed forms for opposed rectangles
(:func:`~engineeringkit.heat_transfer.view_factor_parallel_plates`) and
coaxial disks (:func:`~engineeringkit.heat_transfer.view_factor_coaxial_disks`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.heat_transfer import view_factor_coaxial_disks, view_factor_parallel_plates

rng = np.random.default_rng(1928)

# %%
# The unit-sphere construction as a ray tracer
# --------------------------------------------
#
# A uniform point :math:`(u, v)` in the unit circle lifted to the
# hemisphere is the direction :math:`(u, v, \sqrt{1 - u^2 - v^2})`.


def nusselt_directions(n):
    r, phi = np.sqrt(rng.random(n)), 2 * np.pi * rng.random(n)
    u, v = r * np.cos(phi), r * np.sin(phi)
    return u, v, np.sqrt(1 - u**2 - v**2)


def monte_carlo_plates(a, b, L, n=400_000):
    x, y = a * rng.random(n), b * rng.random(n)
    u, v, w = nusselt_directions(n)
    X, Y = x + L * u / w, y + L * v / w
    return np.mean((X >= 0) & (X <= a) & (Y >= 0) & (Y <= b))


def monte_carlo_disks(r1, r2, L, n=400_000):
    radius, angle = r1 * np.sqrt(rng.random(n)), 2 * np.pi * rng.random(n)
    u, v, w = nusselt_directions(n)
    X, Y = radius * np.cos(angle) + L * u / w, radius * np.sin(angle) + L * v / w
    return np.mean(X**2 + Y**2 <= r2**2)


# %%
# Opposed rectangles and coaxial disks
# ------------------------------------
#
# The Monte Carlo estimate agrees with the closed forms to within its
# statistical error :math:`\sqrt{F(1-F)/n}`.

n = 400_000
cases = [
    ("1 x 1 m plates, 1 m apart", monte_carlo_plates(1, 1, 1), view_factor_parallel_plates(1, 1, 1)),
    ("2 x 1 m plates, 0.5 m apart", monte_carlo_plates(2, 1, 0.5), view_factor_parallel_plates(2, 1, 0.5)),
    ("disks r = 1 to r = 2, 1 m apart", monte_carlo_disks(1, 2, 1), view_factor_coaxial_disks(1, 2, 1)),
    ("disks r = 2 to r = 1, 1 m apart", monte_carlo_disks(2, 1, 1), view_factor_coaxial_disks(2, 1, 1)),
]
for name, mc, exact in cases:
    sigma = np.sqrt(exact * (1 - exact) / n)
    print(f"{name:32s}: unit sphere {mc:.4f}, formula {exact:.4f} ({abs(mc - exact) / sigma:.1f} sigma)")
    assert abs(mc - exact) < 4 * sigma

# %%
# Reciprocity
# -----------
#
# :math:`A_1F_{12} = A_2F_{21}`: the small disk sees the large one far
# better than the reverse.

assert np.isclose(np.pi * 1**2 * view_factor_coaxial_disks(1, 2, 1), np.pi * 2**2 * view_factor_coaxial_disks(2, 1, 1))
print("reciprocity A1 F12 = A2 F21 holds for the disks")

gap = np.geomspace(0.01, 10, 200)
fig, ax = plt.subplots(figsize=(6, 4))
ax.semilogx(gap, view_factor_parallel_plates(1.0, 1.0, gap), label="1 m square plates")
ax.semilogx(gap, view_factor_coaxial_disks(0.5, 0.5, gap), label="1 m diameter disks")
ax.semilogx(gap, view_factor_coaxial_disks(0.5, 1.0, gap), label="disk to 2 m disk")
ax.set(xlabel="separation [m]", ylabel="view factor", title="Parallel surfaces")
ax.legend()
fig.tight_layout()

plt.show()

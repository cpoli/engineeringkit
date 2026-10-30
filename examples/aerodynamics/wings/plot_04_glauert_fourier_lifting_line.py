r"""
Glauert's Fourier solution of the lifting line for any planform (1926)
======================================================================

Prandtl's lifting-line equation is an integro-differential equation for the
spanwise circulation. In *The Elements of Aerofoil and Airscrew Theory*
(1926) Hermann Glauert solved it for any planform. He wrote the circulation
as a Fourier sine series in :math:`\theta`, with
:math:`y = -\tfrac b2\cos\theta`, and enforced the equation at as many
stations as there are terms (:func:`~engineeringkit.aerodynamics.lifting_line`).
The first coefficient fixes the lift and the others add induced drag, so the
span efficiency is :math:`e = 1/(1 + \delta)` with
:math:`\delta = \sum_{n>1} nA_n^2/A_1^2`. This example recovers the elliptic
optimum exactly, then sweeps the taper ratio of straight wings to reproduce
Glauert's induced-drag factor :math:`\delta`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics import finite_wing_cl, lifting_line

b, alpha = 10.0, np.radians(5.0)

# %%
# The elliptic wing has uniform downwash and e = 1
# ------------------------------------------------

elliptic = lifting_line(b, lambda y: 1.2 * np.sqrt(1 - (2 * y / b) ** 2), alpha)
print(f"elliptic: AR = {elliptic.aspect_ratio:.2f}, CL = {elliptic.CL:.4f}, e = {elliptic.span_efficiency:.6f}")
assert np.isclose(elliptic.span_efficiency, 1.0)
assert np.isclose(elliptic.CL, finite_wing_cl(alpha, elliptic.aspect_ratio))

# %%
# Taper ratio and Glauert's delta
# -------------------------------
#
# For a fixed aspect ratio of 8, a straight tapered wing with tip-to-root chord
# ratio :math:`\lambda` near 0.4 loads almost elliptically. A rectangular wing
# (:math:`\lambda = 1`) wastes about 6% and a pointed one (:math:`\lambda = 0`)
# far more.

AR = 8.0
taper = np.linspace(0.0, 1.0, 41)
delta = []
for lam in taper:
    root = 2 * b / (AR * (1 + lam))
    wing = lifting_line(b, lambda y, r=root, t=lam: r * (1 - (1 - t) * np.abs(2 * y / b)), alpha, n_terms=60)
    delta.append(1 / wing.span_efficiency - 1)
delta = np.array(delta)
best = taper[np.argmin(delta)]
print(f"AR = 8: delta = {delta[-1]:.4f} rectangular, minimum {delta.min():.4f} at taper ratio {best:.3f}")
assert 0.25 < best < 0.45 and delta.min() < 0.02

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(taper, delta)
axes[0].set(xlabel=r"taper ratio $\lambda = c_t/c_r$", ylabel=r"induced-drag factor $\delta$", title="AR = 8 straight wings")
for lam in (1.0, 0.35, 0.0):
    root = 2 * b / (AR * (1 + lam))
    wing = lifting_line(b, lambda y, r=root, t=lam: r * (1 - (1 - t) * np.abs(2 * y / b)), alpha, n_terms=60)
    axes[1].plot(2 * wing.y / b, wing.cl / wing.CL, label=rf"$\lambda$ = {lam:g}")
axes[1].axhline(1, color="gray", ls="--", lw=0.8)
axes[1].set(xlabel="$2y/b$", ylabel="$c_l / C_L$", title="Section lift: where each wing stalls first")
axes[1].legend()
fig.tight_layout()

plt.show()

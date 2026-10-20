r"""
Poisson's ratio
===============

A bar stretched along its axis gets thinner. Siméon Denis Poisson's 1829
molecular theory of elasticity predicted that the lateral strain is always
one quarter of the axial strain. Measurements soon showed that the ratio
:math:`\nu` depends on the material: about 0.3 for steel, nearly 0.5 for
rubber and nearly zero for cork, which is why cork seals a wine bottle
without jamming in the neck. The lateral strain is
:math:`-\nu\varepsilon` (:func:`~engineeringkit.structural_engineering.poisson_lateral_strain`),
so the volume changes by :math:`(1 - 2\nu)\varepsilon` to first order: an
incompressible material has :math:`\nu = 1/2`.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.structural_engineering import poisson_lateral_strain

# %%
# Volume change of a stretched bar
# --------------------------------

strain = 1e-3
materials = {"cork": 0.0, "Poisson's 1/4": 0.25, "steel": 0.30, "concrete": 0.20, "rubber": 0.499}
for name, nu in materials.items():
    lateral = poisson_lateral_strain(strain, nu)
    volume = (1 + strain) * (1 + lateral) ** 2 - 1
    print(f"{name:14s} nu = {nu:5.3f}: lateral strain {lateral:+.2e}, volume strain {volume:+.3e}")
    assert abs(volume - (1 - 2 * nu) * strain) <= strain**2  # first-order formula, second-order error

# %%
# The exact volume ratio, across the admissible range
# ---------------------------------------------------

nu = np.linspace(-0.9, 0.49, 200)
eps = 0.01
exact = (1 + eps) * (1 + poisson_lateral_strain(eps, nu)) ** 2 - 1
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(nu, exact / eps, label="exact")
ax.plot(nu, 1 - 2 * nu, "--", label=r"$1 - 2\nu$")
ax.axvline(0.25, color="gray", ls=":", label="Poisson's 1/4")
ax.set(xlabel=r"Poisson's ratio $\nu$", ylabel="volume strain / axial strain", title="1% stretch")
ax.legend()
fig.tight_layout()

plt.show()

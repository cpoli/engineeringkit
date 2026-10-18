r"""
Lamé's elastic constants (1852)
===============================

Early elasticity theories built from molecular forces predicted a single
elastic constant for an isotropic solid, with Poisson's ratio fixed at 1/4.
Measurements disagreed. Gabriel Lamé's *Leçons sur la théorie mathématique de
l'élasticité* (1852) wrote the isotropic law with two constants,

.. math:: \sigma_{ij} = \lambda\,\varepsilon_{kk}\,\delta_{ij} + 2\mu\,\varepsilon_{ij},

and every other modulus follows from any two
(:func:`~engineeringkit.materials.lame_parameters`,
:func:`~engineeringkit.materials.shear_modulus`,
:func:`~engineeringkit.materials.bulk_modulus`). The example applies Lamé's
law to a tensile strain state and plots how the moduli depend on
Poisson's ratio.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import bulk_modulus, lame_parameters, shear_modulus

# %%
# Uniaxial tension from Lamé's law
# --------------------------------
#
# A bar stretched by :math:`\varepsilon` contracts sideways by
# :math:`\nu\varepsilon`. Putting that strain into Lamé's law must give a
# stress :math:`E\varepsilon` along the bar and none across it.

E, nu, eps = 70e9, 0.33, 1e-3
lam, mu = lame_parameters(E, nu)
strain = np.diag([eps, -nu * eps, -nu * eps])
stress = lam * np.trace(strain) * np.eye(3) + 2 * mu * strain
print(f"aluminium: lambda = {lam / 1e9:.1f} GPa, mu = G = {mu / 1e9:.1f} GPa, K = {bulk_modulus(E, nu) / 1e9:.1f} GPa")
print(f"stress from Lame's law: {np.round(np.diag(stress) / 1e6, 6)} MPa")
assert np.allclose(np.diag(stress), [E * eps, 0, 0], atol=1e-3)
assert np.isclose(lam + 2 * mu / 3, bulk_modulus(E, nu)) and mu == shear_modulus(E, nu)

# %%
# The moduli against Poisson's ratio
# ----------------------------------
#
# As :math:`\nu \to 1/2` the material becomes incompressible:
# :math:`\lambda` and :math:`K` diverge while :math:`G \to E/3`. A negative
# :math:`\nu` (auxetic foams) is allowed down to -1, where :math:`G` diverges
# and :math:`K \to 0`.

nus = np.linspace(-0.9, 0.49, 300)
lams, mus = lame_parameters(1.0, nus)
fig, ax = plt.subplots(figsize=(6, 3.8))
ax.semilogy(nus, mus, label="G / E")
ax.semilogy(nus, bulk_modulus(1.0, nus), label="K / E")
ax.semilogy(nus[lams > 0], lams[lams > 0], label="lambda / E")
ax.axvline(0.25, color="0.6", ls=":", label="Poisson's 1/4")
ax.set(xlabel="Poisson's ratio", ylabel="modulus / E", title="Isotropic elastic constants")
ax.legend()
fig.tight_layout()
assert np.allclose(mus * (3 * lams + 2 * mus) / (lams + mus), 1)

plt.show()

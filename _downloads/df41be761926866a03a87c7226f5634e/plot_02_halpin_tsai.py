r"""
The Halpin-Tsai equations (1969)
================================

The Voigt and Reuss bounds are far apart for fibre composites, and the
exact elasticity solutions for fibre arrays were too complicated for
design. John Halpin and Stephen Tsai, at the US Air Force Materials
Laboratory, fitted those solutions with one interpolating formula,

.. math:: \frac{P}{P_m} = \frac{1 + \xi\eta V_f}{1 - \eta V_f}, \qquad
          \eta = \frac{P_f/P_m - 1}{P_f/P_m + \xi}

(:func:`~engineeringkit.materials.halpin_tsai`), in which a single
parameter :math:`\xi` measures how well the reinforcement carries load:
:math:`\xi = 0` is the Reuss bound, :math:`\xi \to \infty` the Voigt bound,
and :math:`\xi = 2` fits the transverse modulus of round fibres. The
example estimates the transverse and shear moduli of glass-epoxy and the
stiffness of short fibres.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import halpin_tsai, rule_of_mixtures_E_longitudinal, rule_of_mixtures_E_transverse

# %%
# Transverse and shear moduli
# ---------------------------
#
# For E-glass (:math:`E_f = 72` GPa, :math:`G_f = 30` GPa) in epoxy
# (:math:`E_m = 3.5` GPa, :math:`G_m = 1.3` GPa) at :math:`V_f = 0.6`,
# Halpin-Tsai gives about 15 GPa across the fibres, well above the
# Reuss value of 8 GPa that springs in series predict.

Ef, Em, Gf, Gm, Vf = 72e9, 3.5e9, 30e9, 1.3e9, 0.6
E2 = halpin_tsai(Ef, Em, Vf, 2.0)
G12 = halpin_tsai(Gf, Gm, Vf, 1.0)
reuss = rule_of_mixtures_E_transverse(Ef, Em, Vf)
print(f"E2 = {E2 / 1e9:.1f} GPa (Reuss {reuss / 1e9:.1f}), G12 = {G12 / 1e9:.2f} GPa")
assert reuss < E2 < rule_of_mixtures_E_longitudinal(Ef, Em, Vf)
assert np.isclose(halpin_tsai(Ef, Em, Vf, 0.0), reuss)

# %%
# One parameter between the bounds
# --------------------------------
#
# For aligned short fibres of aspect ratio :math:`\ell/d`, the longitudinal
# modulus uses :math:`\xi = 2\ell/d`: short stubby fibres sit near the lower
# bound, and long ones approach the Voigt rule of mixtures.

fractions = np.linspace(0, 0.7, 100)
fig, ax = plt.subplots(figsize=(6, 3.8))
ax.plot(fractions, rule_of_mixtures_E_longitudinal(Ef, Em, fractions) / 1e9, "k--", label="Voigt")
for aspect in (1, 10, 100):
    ax.plot(fractions, halpin_tsai(Ef, Em, fractions, 2 * aspect) / 1e9, label=f"l/d = {aspect}")
ax.plot(fractions, rule_of_mixtures_E_transverse(Ef, Em, fractions) / 1e9, "k:", label="Reuss")
ax.set(xlabel="fibre volume fraction", ylabel="modulus [GPa]", title="Halpin-Tsai for glass-epoxy")
ax.legend()
fig.tight_layout()
short, long_ = halpin_tsai(Ef, Em, 0.5, 2.0), halpin_tsai(Ef, Em, 0.5, 200.0)
assert short < long_ < rule_of_mixtures_E_longitudinal(Ef, Em, 0.5)

plt.show()

r"""
Reissner and Stavsky: classical laminate theory (1961)
======================================================

A ply of aligned fibres is stiff along them and soft across, so composite
parts stack plies at several angles. Eric Reissner and Yehuda Stavsky, at
MIT, extended Kirchhoff's plate theory to such stacks and found that in an
unsymmetric laminate stretching and bending are coupled: pulling it makes
it curl. The laminate relates force and moment resultants to midplane
strains and curvatures through three matrices,

.. math:: \begin{pmatrix} \mathbf N \\ \mathbf M \end{pmatrix} =
          \begin{pmatrix} \mathbf A & \mathbf B \\ \mathbf B & \mathbf D \end{pmatrix}
          \begin{pmatrix} \boldsymbol\varepsilon^0 \\ \boldsymbol\kappa \end{pmatrix}

(:func:`~engineeringkit.materials.classical_laminate`), built from each
ply's rotated stiffness (:func:`~engineeringkit.materials.transformed_stiffness`).
The example checks a quasi-isotropic carbon-epoxy layup, finds its ply
stresses, and shows the coupling of an unsymmetric one.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.materials import classical_laminate, lamina_stiffness, transformed_stiffness

# %%
# One ply turned through 90 degrees
# ---------------------------------
#
# Carbon-epoxy has :math:`E_1 = 140` GPa along the fibres and only
# :math:`E_2 = 10` GPa across. The stiffness in the loading direction falls
# from :math:`Q_{11}` to :math:`Q_{22}` as the ply turns.

E1, E2, nu12, G12, t = 140e9, 10e9, 0.3, 5e9, 0.125e-3
Q = lamina_stiffness(E1, E2, nu12, G12)
angles = np.radians(np.linspace(0, 90, 91))
q11 = np.array([transformed_stiffness(Q, a)[0, 0] for a in angles])
print(f"Qbar_11 at 0, 45 and 90 deg: {q11[0] / 1e9:.1f}, {q11[45] / 1e9:.1f}, {q11[90] / 1e9:.1f} GPa")
assert np.isclose(q11[90], Q[1, 1])

# %%
# A quasi-isotropic laminate
# --------------------------
#
# Eight plies :math:`[0/45/-45/90]_s` have :math:`\mathbf B = 0` because
# the stack is symmetric, and an in-plane stiffness that is the same in
# every direction: :math:`A_{11} = A_{22}` and
# :math:`A_{66} = (A_{11} - A_{12})/2`, as for an isotropic sheet.

quasi = classical_laminate(E1, E2, nu12, G12, np.radians([0, 45, -45, 90, 90, -45, 45, 0]), t)
A = quasi.A
summary = quasi.summary()
print(f"Ex = {summary['Ex_Pa'] / 1e9:.1f} GPa, Ey = {summary['Ey_Pa'] / 1e9:.1f} GPa, Gxy = {summary['Gxy_Pa'] / 1e9:.1f} GPa")
assert np.allclose(quasi.B, 0, atol=1e-6) and np.isclose(A[2, 2], (A[0, 0] - A[0, 1]) / 2)

N = (100e3, 0.0, 0.0)  # 100 kN per metre of width along x
stresses = quasi.ply_stresses(N=N)
for angle, (s1, s2, t12) in zip((0, 45, -45, 90), stresses[:4]):
    print(f"{angle:4d} deg ply: sigma_1 = {s1 / 1e6:6.1f} MPa, sigma_2 = {s2 / 1e6:5.1f} MPa, tau_12 = {t12 / 1e6:5.1f} MPa")
assert stresses[0, 0] == stresses[:, 0].max()

fig, axs = plt.subplots(1, 2, figsize=(8.5, 3.6))
axs[0].plot(np.degrees(angles), q11 / 1e9)
axs[0].set(xlabel="ply angle [deg]", ylabel="Qbar_11 [GPa]", title="One ply")
quasi.plot(ax=axs[1])
fig.tight_layout()

# %%
# Stretching makes an unsymmetric laminate curl
# ---------------------------------------------
#
# Two plies, :math:`[0/90]`, are not symmetric about the midplane:
# :math:`\mathbf B \neq 0`, and a pure tension gives a curvature.

cross = classical_laminate(E1, E2, nu12, G12, np.radians([0, 90]), t)
eps0, kappa = cross.midplane_response(N=N)
print(f"[0/90] under tension: strain {eps0[0]:.2e}, curvature {kappa[0]:.2f} 1/m")
assert not cross.summary()["symmetric"] and abs(kappa[0]) > 0

plt.show()

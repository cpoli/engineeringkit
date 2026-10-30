r"""
Daniel Bernoulli's superposition of normal modes (1753)
=======================================================

Daniel Bernoulli argued in 1753 that a vibrating string sounds many simple
vibrations at once, each with its own shape and frequency, and that any
motion is their sum. For masses joined by springs the simple vibrations
are the solutions of the eigenproblem :math:`K\phi = \omega^2 M\phi`
(:func:`~engineeringkit.mechanical_engineering.normal_modes`). The modes are
orthogonal through :math:`M` and :math:`K`, so the equations
:math:`M\ddot x + Kx = 0` split into independent oscillators, one per mode
(:func:`~engineeringkit.mechanical_engineering.modal_response`), and a
harmonic force excites each mode in proportion to its share of the force
(:func:`~engineeringkit.mechanical_engineering.harmonic_response`). The
example builds a three-storey shear frame and checks both superpositions
against direct solutions.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.integrators import integrate_ode
from engineeringkit.mechanical_engineering import harmonic_response, modal_response, normal_modes

# %%
# Modes of a three-storey frame
# -----------------------------
#
# Floors of mass :math:`m` on storeys of lateral stiffness :math:`k`.

m, k = 2.0e4, 3.0e7
M = m * np.diag([1.0, 1.0, 0.5])
K = k * np.array([[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 1.0]])
modes = normal_modes(M, K)
phi = modes.mode_shapes.T
print("natural frequencies [Hz]:", np.round(modes.frequencies / (2 * np.pi), 3).tolist())
assert np.allclose(phi.T @ M @ phi, np.eye(3)) and np.allclose(phi.T @ K @ phi, np.diag(modes.frequencies**2))

# %%
# Free vibration as a sum of modes
# --------------------------------
#
# Displace the top floor alone and let go. Modal superposition gives the
# motion in closed form; a direct integration of :math:`M\ddot x + Kx = 0`
# agrees with it.

x0 = np.array([0.0, 0.0, 0.02])
t = np.linspace(0, 2, 801)
modal = modal_response(M, K, x0, t=t)
A = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.linalg.solve(M, K), np.zeros((3, 3))]])
direct = integrate_ode(lambda _, y: A @ y, np.concatenate((x0, np.zeros(3))), t, rtol=1e-11, atol=1e-14).history[:, :3]
print(f"largest difference from direct integration: {np.max(np.abs(modal.displacement - direct)):.2e} m")
assert np.allclose(modal.displacement, direct, atol=1e-9)
q0 = phi.T @ M @ x0
print("modal participation of the initial shape:", np.round(q0 / np.linalg.norm(q0), 3).tolist())

# %%
# Forced response as a sum of modes
# ---------------------------------
#
# :math:`X = \sum_r \phi_r\,\phi_r^T F/(\omega_r^2 - \omega^2)`.

F = np.array([0.0, 0.0, 1e4])
w = np.linspace(0.5, 80, 2000)
X = harmonic_response(M, K, F, w)
by_modes = (phi @ ((phi.T @ F)[:, None] / (modes.frequencies[:, None] ** 2 - w**2))).T
assert np.allclose(X.real, by_modes) and np.allclose(X.imag, 0)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
modes.plot(ax=ax1)
ax1.set(title="Mode shapes", xlabel="floor")
ax2.semilogy(w / (2 * np.pi), np.abs(X[:, 2]))
for f in modes.frequencies / (2 * np.pi):
    ax2.axvline(f, color="gray", lw=0.5)
ax2.set(xlabel="forcing frequency [Hz]", ylabel="top-floor amplitude [m]", title="Response to a force on the roof")
fig.tight_layout()

plt.show()

r"""
Wonham's pole assignment
========================

For a single input, controllability lets state feedback :math:`u = -Kx`
place the poles of :math:`A - BK` anywhere, and :math:`K` is unique
(Ackermann's formula). Murray Wonham (1967) proved that the same holds with
several inputs: the closed-loop poles can be assigned arbitrarily *if and
only if* :math:`(A, B)` is controllable, even when no single input could do
it (:func:`~engineeringkit.control_systems.state_feedback_poles`). A
precompensator :math:`N` then scales the reference so the output settles on
it (:func:`~engineeringkit.control_systems.feedforward_gain`).

The first plant is an inverted pendulum on a cart, unstable in open loop,
driven by the force on the cart.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.control_systems import (
    feedforward_gain,
    is_controllable,
    state_feedback_poles,
    state_space,
    step_response,
)

M, m, ell, g = 1.0, 0.2, 0.5, 9.81  # cart mass, bob mass, rod length, gravity
A = np.array([[0, 1, 0, 0], [0, 0, -m * g / M, 0], [0, 0, 0, 1], [0, 0, (M + m) * g / (M * ell), 0]])
B = np.array([[0.0], [1 / M], [0.0], [-1 / (M * ell)]])
C = np.array([[1.0, 0.0, 0.0, 0.0]])  # cart position
print("open-loop poles:", np.sort(np.linalg.eigvals(A).real).round(3).tolist())
assert is_controllable(A, B) and np.linalg.eigvals(A).real.max() > 0

desired = [-2.0, -2.5, -3.0, -3.5]
K = state_feedback_poles(A, B, desired)
N = feedforward_gain(A, B, C, K)
assert np.allclose(np.sort(np.linalg.eigvals(A - B @ K).real), sorted(desired))
t = np.linspace(0, 6, 600)
_, y_cart = step_response(state_space(A - B @ K, B * N, C, [[0.0]]), t)
_, y_angle = step_response(state_space(A - B @ K, B * N, [[0.0, 0.0, 1.0, 0.0]], [[0.0]]), t)
print(f"K = {K.round(3).tolist()}, N = {N:.3f}, final cart position {y_cart[-1]:.4f} m")
assert np.isclose(y_cart[-1], 1, atol=1e-3) and abs(y_angle[-1]) < 1e-3

# %%
# Two inputs that are each insufficient
# -------------------------------------
#
# Two unstable modes, each driven by its own actuator. Neither input alone
# reaches both modes, but together they are controllable, so Wonham's
# theorem promises any pole set, here a complex pair.

A2 = np.diag([1.0, 2.0])
B2 = np.eye(2)
assert not is_controllable(A2, B2[:, :1]) and not is_controllable(A2, B2[:, 1:])
assert is_controllable(A2, B2)
K2 = state_feedback_poles(A2, B2, [-1 + 2j, -1 - 2j])
assert np.allclose(np.sort_complex(np.linalg.eigvals(A2 - B2 @ K2)), [-1 - 2j, -1 + 2j])
print("two-input gain:", K2.round(3).tolist())

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
ax1.plot(t, y_cart)
ax1.axhline(1, color="k", lw=0.8)
ax1.set(ylabel="cart position [m]", title=f"Poles placed at {desired}")
ax2.plot(t, np.degrees(y_angle), color="C1")
ax2.set(xlabel="time [s]", ylabel="pendulum angle [deg]")

plt.show()

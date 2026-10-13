import numpy as np

from engineeringkit.control_systems import *


def test_routh():
    assert routh_stable([1, 3, 3, 1])


def test_state_space_properties():
    A = np.array([[0, 1], [-2, -3.0]])
    B = np.array([[0], [1.0]])
    C = np.array([[1.0, 0]])
    assert is_controllable(A, B) and is_observable(A, C)
    K = state_feedback_poles(A, B, [-4, -5])
    assert np.allclose(np.sort(np.linalg.eigvals(A - B @ K)), [-5, -4])


def test_lqr_kalman_stable():
    A = np.array([[0, 1], [-2, -0.5]])
    B = np.array([[0], [1.0]])
    C = np.array([[1.0, 0]])
    K, _, e = lqr(A, B, np.eye(2), np.array([[1.0]]))
    assert np.all(np.real(e) < 0)
    L, _, e2 = kalman_gain(A, C, np.eye(2), np.array([[1.0]]))
    assert np.all(np.real(e2) < 0)


def test_second_order_metrics():
    assert np.isclose(second_order_metrics(0.5, 2)["settling_time_2pct"], 4)

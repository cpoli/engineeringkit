import numpy as np
import pytest

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


def test_repeated_poles_single_input_use_ackermann():
    A = np.array([[0, 1, 0], [0, 0, 1], [-1, -5, -6.0]])
    B = np.array([[0], [0], [1.0]])
    K = state_feedback_poles(A, B, [-2, -2, -2])
    np.testing.assert_allclose(K, [[7, 7, 0]], atol=1e-12)
    np.testing.assert_allclose(np.poly(A - B @ K), [1, 6, 12, 8], atol=1e-9)
    # With two inputs place_poles itself can repeat a pole twice.
    np.testing.assert_allclose(state_feedback_poles(np.zeros((2, 2)), np.eye(2), [-1, -1]), np.eye(2))
    for bad in ([-1, -1], [-1j, -1j, -1]):
        with pytest.raises(ValueError, match="conjugation"):
            state_feedback_poles(A, B, bad)
    with pytest.raises(ValueError, match="not controllable"):
        state_feedback_poles([[-1, 0], [0, -2]], [[1], [0]], [-3, -3])


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

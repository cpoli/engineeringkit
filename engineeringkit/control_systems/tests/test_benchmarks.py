import numpy as np
import pytest

from engineeringkit import control_systems as c


def test_scalar_lqr_against_analytic_optimum():
    K, P, poles = c.lqr([[0]], [[1]], [[4]], [[1]])
    np.testing.assert_allclose(K, [[2]])
    np.testing.assert_allclose(P, [[2]])
    np.testing.assert_allclose(poles, [-2])


def test_matrix_riccati_residual_and_cost_symmetry():
    A = np.array([[0, 1], [-2, -0.5]])
    B = np.array([[0], [1.0]])
    Q = np.diag([4.0, 1.0])
    R = np.array([[2.0]])
    K, P, _ = c.lqr(A, B, Q, R)
    residual = A.T @ P + P @ A - P @ B @ np.linalg.solve(R, B.T @ P) + Q
    np.testing.assert_allclose(residual, 0, atol=1e-12)
    np.testing.assert_allclose(P, P.T, atol=1e-12)
    np.testing.assert_allclose(R @ K, B.T @ P, atol=1e-12)


def test_scalar_kalman_analytic_covariance():
    L, P, eigen = c.kalman_gain([[0]], [[1]], [[4]], [[1]])
    np.testing.assert_allclose(L, [[2]])
    np.testing.assert_allclose(P, [[2]])
    np.testing.assert_allclose(eigen, [-2])


@pytest.mark.parametrize("Q,R", [([[-1]], [[1]]), ([[1]], [[0]]), ([[1]], [[np.nan]])])
def test_invalid_riccati_weights(Q, R):
    for fn in (c.lqr, c.kalman_gain):
        with pytest.raises(ValueError):
            fn([[0]], [[1]], Q, R)


def test_state_matrices_accept_lists_and_reject_shapes():
    assert c.is_controllable([[0]], [[1]])
    assert c.is_observable([[0]], [[1]])
    with pytest.raises(ValueError):
        c.controllability_matrix([[0, 1]], [[1]])
    with pytest.raises(ValueError):
        c.observability_matrix([[0]], [[1, 0]])

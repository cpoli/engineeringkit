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


def test_tustin_against_scipy_bilinear():
    from scipy import signal

    num, den = [2.0, 3.0, 1.0], [1.0, 4.0, 5.0, 2.0]
    dt = 0.05
    Gd = c.discretize_tustin(c.transfer_function(num, den), dt)
    bz, az = signal.bilinear(num, den, fs=1 / dt)
    np.testing.assert_allclose(Gd.num, bz / az[0], rtol=1e-10, atol=1e-14)
    np.testing.assert_allclose(Gd.den, az / az[0], rtol=1e-10)
    # Prewarping at w0 is the bilinear transform with fs = w0 / (2 tan(w0 dt / 2)).
    w0 = 20.0
    Gp = c.discretize_tustin(c.transfer_function(num, den), dt, prewarp_frequency=w0)
    bz, az = signal.bilinear(num, den, fs=w0 / (2 * np.tan(w0 * dt / 2)))
    np.testing.assert_allclose(Gp.den, az / az[0], rtol=1e-10)


@pytest.mark.parametrize("p", [-1.0, 0.5, 2.0])
def test_bode_sensitivity_integral(p):
    # Relative degree two: the integral of ln|S| over (0, inf) is pi times the sum of the unstable open-loop poles.
    from scipy import integrate

    L = c.transfer_function([30, 30], np.polymul([1, -p], [1, 7, 10]))
    assert np.all(c.poles(c.feedback(L)).real < 0)
    integrand = lambda w: np.log(abs(1 / (1 + c.nyquist(L, w))))  # noqa: E731
    value = sum(integrate.quad(integrand, a, b, limit=400)[0] for a, b in [(0, 1), (1, 100), (100, np.inf)])
    assert value == pytest.approx(np.pi * max(p, 0), abs=2e-4)

import numpy as np
import pytest

from engineeringkit.integrators import integrate_ode


@pytest.mark.parametrize("method", ["euler", "rk4", "RK45", "Radau"])
def test_exponential_decay(method):
    t = np.linspace(0, 2, 401)
    r = integrate_ode(lambda t, y: -y, [1.0], t, method)
    tolerance = 0.002 if method == "euler" else 1e-7
    np.testing.assert_allclose(r.history[:, 0], np.exp(-t), atol=tolerance)
    assert r.history.shape == (401, 1)
    assert r.evaluations > 0


def test_rk4_fourth_order_convergence():
    errors = []
    for n in (11, 21):
        r = integrate_ode(lambda t, y: y, [1.0], np.linspace(0, 1, n), "rk4")
        errors.append(abs(r.history[-1, 0] - np.e))
    assert 14 < errors[0] / errors[1] < 18


@pytest.mark.parametrize("times", [[0], [0, 0], [1, 0], [0, np.nan], [[0, 1]]])
def test_invalid_grid(times):
    with pytest.raises(ValueError):
        integrate_ode(lambda t, y: y, [1.0], times)


def test_invalid_derivative():
    with pytest.raises(ValueError, match="derivative"):
        integrate_ode(lambda t, y: [1, 2], [1.0], [0, 1])


def test_solver_failure_is_visible():
    with pytest.raises(RuntimeError):
        integrate_ode(lambda t, y: y * y, [1.0], [0, 2])

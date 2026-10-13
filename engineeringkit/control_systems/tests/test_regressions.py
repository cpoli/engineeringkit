import numpy as np
import pytest
from scipy import signal

from engineeringkit.control_systems import PID, root_locus, routh_first_column, routh_stable


@pytest.mark.parametrize(
    "coeffs,stable",
    [([1, 3, 2], True), ([-1, -3, -2], True), ([1, 0, 1], False), ([1, 0, 2, 0, 1], False), ([1, 1, -2], False), ([1, 0], False)],
)
def test_routh_strict_stability(coeffs, stable):
    assert routh_stable(coeffs) == stable
    assert np.all(np.isfinite(routh_first_column(coeffs)))


def test_proportional_pid_has_no_spurious_integrator():
    plant = signal.TransferFunction([1], [1, 1])
    closed = PID(2).closed_loop(plant)
    np.testing.assert_allclose(closed.den, [1, 3])
    r = PID(2).simulate(plant, np.linspace(0, 5, 101))
    np.testing.assert_allclose(r.y, 2 / 3 * (1 - np.exp(-3 * r.t)), atol=1e-12)


def test_derivative_pid_handles_polynomial_degrees():
    closed = PID(1, 0, 1).closed_loop(signal.TransferFunction([1], [1, 1]))
    np.testing.assert_allclose(closed.num, [0.5, 0.5])
    np.testing.assert_allclose(closed.den, [1, 1])


def test_root_locus_known_characteristic_roots():
    roots = root_locus([1], [1, 3, 2], [0, 2])
    np.testing.assert_allclose(np.sort(roots[0].real), [-2, -1])
    np.testing.assert_allclose(np.sort(roots[1].imag), [-np.sqrt(7) / 2, np.sqrt(7) / 2])
    with pytest.raises(ValueError):
        root_locus([1, 1], [1, 2], [0])

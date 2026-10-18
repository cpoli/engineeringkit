import numpy as np
import pytest

from engineeringkit.integrators import integrate_ode


def test_summary_reports_evaluations_and_final_state():
    r = integrate_ode(lambda t, y: -y, [1], [0, 0.1], "euler")
    assert r.summary() == {"method": "euler", "samples": 2, "evaluations": 1, "final_state": [0.9]}


@pytest.mark.parametrize("state", [[], [np.nan], [[1]]])
def test_invalid_initial_state(state):
    with pytest.raises(ValueError):
        integrate_ode(lambda t, y: y, state, [0, 1])


@pytest.mark.parametrize("kwargs", [{"rtol": 0}, {"atol": np.inf}])
def test_invalid_tolerance(kwargs):
    with pytest.raises(ValueError):
        integrate_ode(lambda t, y: y, [1], [0, 1], **kwargs)


def test_explicit_overflow_is_detected():
    # Deliberately exceed float64 range; the solver must reject its next state.
    with np.errstate(over="ignore"):
        with pytest.raises(RuntimeError, match="non-finite"):
            integrate_ode(lambda t, y: np.array([1e308]), [1e308], [0, 10], "euler")

import numpy as np
import pytest

from engineeringkit import control_systems as c


def test_classical_first_order_responses_and_frequency():
    sys = c.transfer_function([2, 2], [1, 3, 2])
    np.testing.assert_allclose(np.sort(c.poles(sys)), [-2, -1])
    np.testing.assert_allclose(c.zeros(sys), [-1])
    assert c.dc_gain(sys) == 1
    times = np.linspace(0, 5, 101)
    ts, y = c.step_response(sys, times)
    np.testing.assert_allclose(y, 1 - np.exp(-2 * times), atol=1e-12)
    _, impulse = c.impulse_response(sys, times)
    np.testing.assert_allclose(impulse, 2 * np.exp(-2 * times), atol=1e-12)
    w, mag, phase = c.bode(sys, [2])
    assert np.isclose(mag[0], 20 * np.log10(1 / np.sqrt(2)))
    assert np.isclose(phase[0], -45)
    np.testing.assert_allclose(c.nyquist(sys, [2]), [0.5 - 0.5j])
    ss = c.state_space([[-2]], [[2]], [[1]], [[0]])
    assert (-ss.C @ np.linalg.solve(ss.A, ss.B) + ss.D).item() == 1
    assert c.feedforward_gain([[-1]], [[1]], [[1]], [[2]]) == 3


def test_pid_integral_response_summary_and_plot():
    import matplotlib.pyplot as plt

    plant = c.transfer_function([1], [1, 1])
    r = c.PID(1, 1).simulate(plant, np.linspace(0, 5, 101))
    np.testing.assert_allclose(r.y, 1 - np.exp(-r.t), atol=1e-12)
    assert r.summary()["final"] == r.y[-1]
    fig, ax = r.plot()
    assert ax.get_xlabel() == "time [s]"
    plt.close(fig)
    with pytest.raises(ValueError):
        c.PID(np.nan).closed_loop(plant)


def test_frequency_margins_at_known_crossovers():
    # 10/(s*(s+1)*(s+2)): phase crossover sqrt(2), gain margin 0.6.
    sys = c.transfer_function([10], [1, 3, 2, 0])
    margins = c.stability_margins(sys)
    assert np.isclose(margins["phase_cross_frequency"], np.sqrt(2), rtol=0.002)
    assert np.isclose(margins["gain_margin"], 0.6, rtol=0.003)
    assert np.isfinite(margins["phase_margin_deg"])
    no_cross = c.stability_margins(c.transfer_function([0.1], [1, 1]), np.logspace(-2, 2, 101))
    assert no_cross["gain_margin"] == np.inf and np.isnan(no_cross["phase_margin_deg"])


def test_matrix_and_routh_special_failures():
    with pytest.raises(ValueError):
        c.controllability_matrix([[0]], [[1], [2]])
    with pytest.raises(ValueError):
        c.lqr(np.eye(2), np.eye(2), [[1, 1], [0, 1]], np.eye(2))
    with pytest.raises(ValueError):
        c.routh_first_column([0, 1])
    # An isolated zero pivot invokes epsilon; it does not imply stability.
    col = c.routh_first_column([1, 0, 1, 1])
    assert np.all(np.isfinite(col)) and np.any(col < 0)
    with pytest.raises(ValueError):
        c.second_order_metrics(1, 1)


def test_nyquist_encirclements_count_and_contour_checks():
    # 12/(s+1)^3 is unstable in closed loop (K_u = 8): two clockwise encirclements.
    assert c.nyquist_encirclements(c.transfer_function([12], [1, 3, 3, 1])) == 2
    assert c.nyquist_encirclements(c.transfer_function([4], [1, 3, 3, 1]), np.logspace(-3, 3, 4001)) == 0
    with pytest.raises(ValueError, match="proper"):
        c.nyquist_encirclements(c.transfer_function([1, 0, 0], [1, 1]))
    with pytest.raises(ValueError, match="imaginary axis"):
        c.nyquist_encirclements(c.transfer_function([1], [1, 1, 0]))
    with pytest.raises(ValueError, match="marginally stable"):
        c.nyquist_encirclements(c.transfer_function([8], [1, 3, 3, 1]), np.sqrt(3) * np.array([0.5, 1.0, 2.0]))
    with pytest.raises(ValueError):
        c.nyquist(c.transfer_function([1], [1, 1]), [np.nan])


def test_hurwitz_determinants_and_lyapunov_stable_errors():
    np.testing.assert_allclose(c.hurwitz_determinants([2, 6, 6, 2]), [3, 8, 8])  # normalised (s+1)^3
    assert c.hurwitz_determinants([1, 2, 1, 2])[1] == 0  # roots at +-j and -2
    with pytest.raises(ValueError):
        c.lyapunov_stable([[1.0, 2.0]])  # a non-square matrix is not a 'singular equation'


def test_discrete_lqr_and_zoh_against_scalar_closed_forms():
    Ad, Bd = c.discretize_zoh([[-1.0]], [[1.0]], 0.5)
    assert np.isclose(Ad[0, 0], np.exp(-0.5)) and np.isclose(Bd[0, 0], 1 - np.exp(-0.5))
    K, P, poles = c.dlqr([[1.0]], [[1.0]], [[1.0]], [[1.0]])
    golden = (1 + np.sqrt(5)) / 2
    assert np.isclose(P[0, 0], golden) and np.isclose(K[0, 0], golden / (1 + golden))
    assert abs(poles[0]) < 1


def test_state_space_contracts():
    with pytest.raises(ValueError):
        c.state_space([[0]], [[1]], [[1]], [[0, 0]])
    L = c.observer_gain([[0, 1], [-2, -3]], [[1, 0]], [-5, -6])
    np.testing.assert_allclose(np.sort(np.linalg.eigvals(np.array([[0, 1], [-2, -3]]) - L @ [[1, 0]]).real), [-6, -5])
    with pytest.raises(ValueError, match="SISO"):
        c.feedforward_gain(np.eye(2), np.eye(2), [[1, 0]], np.zeros((2, 2)))
    with pytest.raises(ValueError, match="zero DC gain"):
        c.feedforward_gain([[-1, 0], [0, -2]], [[1], [1]], [[1, -2]], [[0, 0]])  # zero at s = 0


def test_feedback_dc_gain_pid_and_model_conversion_errors():
    with pytest.raises(ValueError, match="sign"):
        c.feedback(c.transfer_function([1], [1, 1]), sign=0)
    with pytest.raises(ValueError, match="unbounded"):
        c.dc_gain(c.transfer_function([1], [1, 0]))
    with pytest.raises(ValueError, match="all gains zero"):
        c.PID(0).transfer_function()
    with pytest.raises(ValueError):
        c.ziegler_nichols(8, 4, controller="PD")
    num, den = c.num_den(c.state_space([[-2]], [[2]], [[1]], [[0]]))
    np.testing.assert_allclose(num, [2]) and np.testing.assert_allclose(den, [1, 2])
    with pytest.raises(ValueError, match="single-input"):
        c.num_den(c.state_space(-np.eye(2), np.eye(2), np.eye(2), np.zeros((2, 2))))

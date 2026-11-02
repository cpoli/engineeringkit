"""Worked examples from K. Ogata, *Modern Control Engineering*, 5th ed. (Prentice Hall, 2010).

Each test reproduces the answer the book prints, to the precision it prints
it. Where Ogata reads a value off a plot or rounds an intermediate result,
the test says so and the control_systems guide records the difference.
"""

import numpy as np

from engineeringkit import control_systems as c


def test_ogata_example_5_1_second_order_specifications():
    # zeta = 0.6, omega_n = 5 rad/s: tp = 0.785 s, Mp = 0.095, ts = 1.33 s (2 %).
    m = c.second_order_metrics(0.6, 5.0)
    assert round(float(m["peak_time"]), 3) == 0.785
    assert round(float(m["percent_overshoot"]) / 100, 3) == 0.095
    assert round(float(m["settling_time_2pct"]), 2) == 1.33


def test_ogata_example_5_12_routh_array():
    # s^4 + 2s^3 + 3s^2 + 4s + 5: first column 1, 2, 1, -6, 5, two sign changes.
    column = c.routh_first_column([1, 2, 3, 4, 5])
    np.testing.assert_allclose(column, [1, 2, 1, -6, 5])
    assert np.count_nonzero(np.diff(np.sign(column))) == 2
    assert np.count_nonzero(np.roots([1, 2, 3, 4, 5]).real > 0) == 2
    assert not c.routh_stable([1, 2, 3, 4, 5])


def test_ogata_example_6_1_root_locus():
    # K/(s(s+1)(s+2)): breakaway at s = -0.4226 with K = 0.3849, and the
    # locus crosses the imaginary axis at s = +-j sqrt(2) with K = 6.
    den = [1, 3, 2, 0]
    crossing = c.root_locus([1], den, [6.0])[0]
    imaginary = np.sort(crossing[np.abs(crossing.real) < 1e-9].imag)
    np.testing.assert_allclose(imaginary, [-np.sqrt(2), np.sqrt(2)])
    assert round(np.sqrt(2), 4) == round(float(imaginary[1]), 4) == 1.4142
    assert c.routh_stable([1, 3, 2, 5.99]) and not c.routh_stable([1, 3, 2, 6])

    breakaway = np.sort(c.root_locus([1], den, [0.3849])[0].real)[1:]
    assert abs(breakaway[1] - breakaway[0]) < 1e-3
    assert round(float(breakaway.mean()), 4) == -0.4226
    assert np.count_nonzero(c.root_locus([1], den, [0.39])[0].imag) == 2


def test_ogata_example_7_26_lead_compensator():
    # G = 4/(s(s+2)) with K = 10 for Kv = 20 /s, phase margin >= 50 deg.
    # Ogata reads a 17 deg margin off the Bode plot (exactly 17.96 deg) and
    # adds 5 deg, so needs 38 deg of lead: alpha = 0.24, crossover 9 rad/s,
    # Gc = 41.7 (s + 4.41)/(s + 18.4).
    loop = c.transfer_function([40], [1, 2, 0])
    assert round(c.error_constants(loop)["Kv"]) == 20
    pm0 = c.stability_margins(loop)["phase_margin_deg"]
    assert abs(pm0 - 17) < 1

    # Ogata's 38 deg of lead is design_lead's 50 - pm0 + extra_deg.
    r = c.design_lead(loop, 50.0, extra_deg=38 - (50 - pm0))
    alpha = 1 / r.summary()["high_frequency_gain"]
    zero, pole = -r.compensator.zeros[0], -r.compensator.poles[0]
    assert round(alpha, 2) == 0.24
    assert round(r.gain_cross_frequency) == 9
    # Ogata places the zero and pole at 9 sqrt(0.24) = 4.41 and
    # 9/sqrt(0.24) = 18.4 from the rounded crossover; the exact crossover,
    # 8.95 rad/s, and alpha = 0.238 put them at 4.36 and 18.34.
    assert abs(zero - 4.41) < 0.06
    assert abs(pole - 18.4) < 0.06
    assert r.phase_margin_deg >= 50

    # Ogata's printed compensator gives Kv = 20 /s and about 50 deg of margin.
    compensated = c.transfer_function(np.polymul([4 * 41.7], [1, 4.41]), np.polymul([1, 18.4], [1, 2, 0]))
    assert round(c.error_constants(compensated)["Kv"], 1) == 20.0
    m = c.stability_margins(compensated)
    assert abs(m["phase_margin_deg"] - 50) < 1
    assert m["gain_margin"] == np.inf


def test_ogata_example_8_1_ziegler_nichols_second_method():
    # 1/(s(s+1)(s+5)): Kcr = 30, Pcr = 2.8099 s, so Kp = 18, Ti = 1.405 s,
    # Td = 0.35124 s and Gc = 6.3223 (s + 1.4235)^2 / s.
    m = c.stability_margins(c.transfer_function([1], [1, 6, 5, 0]))
    kcr, pcr = m["gain_margin"], 2 * np.pi / m["phase_cross_frequency"]
    assert round(kcr, 4) == 30
    assert round(pcr, 4) == 2.8099

    pid = c.ziegler_nichols(kcr, pcr)
    assert round(float(pid.kp), 4) == 18
    assert round(float(pid.kp / pid.ki), 3) == 1.405
    assert round(float(pid.kd / pid.kp), 5) == 0.35124
    assert round(float(pid.kd), 4) == 6.3223
    double_zero = -float(pid.kp) / (2 * float(pid.kd))
    assert round(double_zero, 4) == -1.4235
    assert abs(float(pid.kp) ** 2 - 4 * float(pid.kd) * float(pid.ki)) < 1e-4


def test_ogata_example_10_1_pole_placement():
    # Closed-loop poles -2 +- j4 and -10 give K = [199 55 8].
    A = [[0, 1, 0], [0, 0, 1], [-1, -5, -6]]
    B = [[0], [0], [1]]
    K = c.state_feedback_poles(A, B, [-2 + 4j, -2 - 4j, -10])
    np.testing.assert_allclose(K, [[199, 55, 8]], atol=1e-9)
    closed = np.sort_complex(np.linalg.eigvals(np.asarray(A) - np.asarray(B) @ K))
    np.testing.assert_allclose(closed, [-10, -2 - 4j, -2 + 4j], atol=1e-9)


def test_ogata_example_10_6_full_order_observer_with_a_double_pole():
    # A = [[0, 20.6], [1, 0]], C = [0 1], observer poles -10, -10: Ke = [120.6, 20].
    # SciPy's place_poles refuses a pole repeated more often than rank(C), so
    # this needs the Ackermann fallback for one output.
    A = [[0, 20.6], [1, 0]]
    L = c.observer_gain(A, [[0, 1]], [-10, -10])
    np.testing.assert_allclose(L, [[120.6], [20]], atol=1e-9)
    np.testing.assert_allclose(np.poly(np.asarray(A) - L @ [[0, 1]]), [1, 20, 100], atol=1e-9)

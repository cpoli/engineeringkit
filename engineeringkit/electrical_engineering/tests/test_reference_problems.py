"""Worked examples from Glover, Overbye and Sarma, checked to the precision printed.

J. D. Glover, T. J. Overbye and M. S. Sarma, *Power System Analysis and
Design*, 6th ed. (Boston: Cengage, 2017), and J. D. Glover, T. J. Overbye,
A. B. Birchfield and M. S. Sarma, 7th ed. (Boston: Cengage, 2022), which
moves transient stability from chapter 11 to chapter 12. Each test uses the
book's per-unit inputs and accepts half a unit in the last printed digit,
except Example 12.4, whose printed root the guide explains.

The printed answers were not read from a copy of either book: Examples 6.9,
12.4 and 12.5 from Overbye's ECEN 460 slides (Texas A&M, 2019 and 2025),
Example 9.3 from K. Y. Lee's ELC 4340 slides (Baylor, 2022, which teach
from the 6th edition) and Example 11.1 from two reproductions of it.
Hayt, Kemmerly and Durbin, the other book named in ROADMAP D4, could not
be checked against a copy and is left for later.
"""

import numpy as np
import pytest

from engineeringkit import electrical_engineering as e


def _five_bus_case():
    # 7th ed. Example 6.9 (the five-bus case of 6th ed. ch. 6), 100 MVA base; buses 1..5 are 0..4 here.
    # Lines 2-4, 2-5 and 4-5 with total charging B'; transformers 1-5 and 3-4.
    lines = [
        (1, 3, 0.0090 + 0.100j, 1.72j),
        (1, 4, 0.0045 + 0.050j, 0.88j),
        (3, 4, 0.00225 + 0.025j, 0.44j),
        (0, 4, 0.00150 + 0.02j),
        (2, 3, 0.00075 + 0.01j),
    ]
    return e.admittance_matrix(5, lines)


def test_example_6_9_five_bus_ybus_entry():
    Y = _five_bus_case()
    assert Y[1, 1].real == pytest.approx(2.67828, abs=5e-5)
    assert Y[1, 1].imag == pytest.approx(-28.4590, abs=5e-4)


def test_example_6_9_five_bus_newton_raphson():
    # Bus 2 load 8.0 + j2.8 pu; bus 3 generates 5.2 pu at |V| = 1.05 and loads 0.8 + j0.4 pu.
    r = e.newton_raphson_load_flow(_five_bus_case(), [8.0, 0.8 - 5.2, 0.0, 0.0], [2.8, 0.4, 0.0, 0.0], pv_buses={2: 1.05})
    assert r.converged
    np.testing.assert_allclose(np.abs(r.voltages), [1.000, 0.834, 1.050, 1.019, 0.974], atol=5e-4)
    np.testing.assert_allclose(np.degrees(np.angle(r.voltages)), [0.0, -22.406, -0.597, -2.834, -4.548], atol=5e-4)
    # Slack generation 395 MW, 114 Mvar; bus 3 generator reactive output 337 Mvar.
    assert 100 * r.power[0].real == pytest.approx(395, abs=0.5)
    assert 100 * r.power[0].imag == pytest.approx(114, abs=0.5)
    assert 100 * (r.power[2].imag + 0.4) == pytest.approx(337, abs=0.5)


def test_example_9_3_bolted_single_line_to_ground_fault():
    # 6th ed. Example 9.3, fault at bus 2 of Example 9.1: Vf = 1.05 pu, Z1 = j0.13893, Z2 = j0.14562, Z0 = j0.25.
    ia, ib, ic = e.fault_currents(1.05, 0.13893j, 0.14562j, 0.25j, kind="slg")
    i0, i1, i2 = e.symmetrical_components(ia, ib, ic)
    assert i1.real == pytest.approx(0, abs=1e-12)
    assert i1.imag == pytest.approx(-1.964, abs=5e-4)
    assert i0 == pytest.approx(i1, abs=1e-12) and i2 == pytest.approx(i1, abs=1e-12)
    assert ia == pytest.approx(3 * i1, abs=1e-12)
    assert abs(ib) < 1e-12 and abs(ic) < 1e-12


def test_example_11_1_swing_angle_three_cycles_after_terminal_fault():
    # 6th ed. Example 11.1: H = 2.0 pu-s, 60 Hz, pm = pe = 1.0 at delta = 10 degrees; a bolted
    # terminal fault drops pe to zero. The angle three cycles (0.05 s) later is 16.75 degrees.
    p_pre = 1 / np.sin(np.radians(10))
    r = e.swing_response(2.0, 1.0, p_pre, 0.0, p_pre, clearing_time=0.05, t_end=0.1)
    assert np.degrees(r.delta[0]) == pytest.approx(10, abs=1e-9)
    assert np.degrees(r.delta[np.searchsorted(r.t, 0.05)]) == pytest.approx(16.75, abs=0.005)


def test_example_12_4_maximum_swing_after_cleared_fault():
    # 7th ed. Example 12.4: pm = 1.0, Pmax = 2.4638 pu before and after a bolted fault at bus 1,
    # delta0 = 0.4179 rad, cleared at delta1 = 0.4964 rad; the rotor swings to delta2 = 0.7003 rad.
    # The slide gives only delta1; H = 3 s cleared at 0.05 s (three cycles) reaches it, and delta2
    # depends on delta1 alone.
    r = e.swing_response(3.0, 1.0, 2.4638, 0.0, 2.4638, clearing_time=0.05, t_end=1.0, n_points=20001)
    assert r.delta[0] == pytest.approx(0.4179, abs=5e-5)
    assert r.delta[np.searchsorted(r.t, 0.05)] == pytest.approx(0.4964, abs=1e-4)
    assert r.stable
    # 0.70025 rad; the book's iterative root of its own equal-area equation is 0.70020 (see the guide).
    assert r.delta.max() == pytest.approx(0.7003, abs=1e-4)


def test_example_12_5_critical_clearing_time():
    # 7th ed. Example 12.5, the system of Example 12.4: the critical clearing time is slightly less than 0.19 s.
    delta0, _ = e.equal_area_clearing_angle(1.0, 2.4638, 0.0, 2.4638)
    assert delta0 == pytest.approx(0.4179, abs=5e-5)
    t_cr = e.critical_clearing_time(3.0, 1.0, 2.4638)
    assert 0.185 < t_cr < 0.19
    assert e.swing_response(3.0, 1.0, 2.4638, 0.0, 2.4638, t_cr - 1e-3, 2.0).stable
    assert not e.swing_response(3.0, 1.0, 2.4638, 0.0, 2.4638, t_cr + 1e-3, 2.0).stable

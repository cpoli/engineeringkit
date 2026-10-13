import numpy as np
import pytest
from scipy import signal

from engineeringkit.control_systems import *
from engineeringkit.electrical_engineering import *
from engineeringkit.fluid_engineering import *
from engineeringkit.structural_engineering import *
from engineeringkit.thermodynamics import *


def test_hooke():
    assert hooke_stress(200e9, 1e-3) == pytest.approx(200e6)


def test_beam_midpoint_point_load():
    b = SimplySupportedBeam(4, 200e9, 8e-6)
    r = b.solve_point_load(1000)
    assert r.reactions == (500, 500)
    assert np.max(r.moment) == pytest.approx(1000)
    assert np.min(r.deflection) == pytest.approx(-1000 * 4**3 / (48 * 200e9 * 8e-6), rel=2e-4)


def test_beam_udl_moment():
    r = SimplySupportedBeam(6, 200e9, 8e-6).solve_uniform_load(1000)
    assert np.max(r.moment) == pytest.approx(1000 * 6**2 / 8, rel=1e-4)


def test_buckling():
    assert euler_buckling_load(200e9, 1e-6, 2) == pytest.approx(np.pi**2 * 200e9 * 1e-6 / 4)


def test_mohr():
    assert mohr_circle(100, 0, 0)[2:] == (100, 0)


def test_control_state_space():
    A = np.array([[0, 1], [-2, -3.0]])
    B = np.array([[0], [1.0]])
    C = np.array([[1.0, 0]])
    assert is_controllable(A, B)
    assert is_observable(A, C)


def test_lqr_stabilizes():
    A = np.array([[0, 1], [1, 0.0]])
    B = np.array([[0], [1.0]])
    K, _, eig = lqr(A, B, np.eye(2), np.eye(1))
    assert np.all(eig.real < 0)


def test_routh():
    assert routh_stable([1, 3, 3, 1])


def test_pid_runs():
    assert PID(2, 1, 0.1).simulate(signal.TransferFunction([1], [1, 1])).y[-1] > 0.8


def test_ohm():
    assert ohms_law(voltage=12, resistance=6) == 2


def test_parallel():
    assert parallel_resistance(10, 10) == 5


def test_resonance():
    assert resonance_frequency(1, 1) == pytest.approx(1 / (2 * np.pi))


def test_three_phase():
    assert three_phase_power(400, 10, 0.8) == pytest.approx(np.sqrt(3) * 3200)


def test_otto():
    assert otto_efficiency(10) == pytest.approx(1 - 10**-0.4)


def test_carnot():
    assert carnot_efficiency(600, 300) == 0.5


def test_first_law():
    assert first_law_heat(100, 20) == 120


def test_reynolds():
    assert reynolds_number(1000, 1, 0.01, 0.001) == pytest.approx(10000)


def test_poiseuille():
    assert poiseuille_flow(8 / np.pi, 1, 1, 1) == pytest.approx(1)


def test_colebrook_laminar():
    assert colebrook_friction_factor(1000, 0) == pytest.approx(0.064)


def test_continuity():
    assert continuity_velocity(2, 3, 1) == 6

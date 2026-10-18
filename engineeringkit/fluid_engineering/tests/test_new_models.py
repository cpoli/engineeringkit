"""Analytical checks for the models added in the fluid_engineering rebuild."""

import numpy as np
import pytest
from scipy import integrate, optimize

from engineeringkit import fluid_engineering as f
from engineeringkit.fluid_engineering.utils import integer_null_space, rational_rref


def test_plane_force_and_centre_of_pressure_match_pressure_integrals():
    rho, top, height, width = 1000.0, 2.0, 3.0, 1.5
    force = integrate.quad(lambda y: f.hydrostatic_pressure(rho, y) * width, top, top + height)[0]
    moment = integrate.quad(lambda y: y * f.hydrostatic_pressure(rho, y) * width, top, top + height)[0]
    yc, area = top + height / 2, width * height
    assert f.hydrostatic_force_plane(rho, yc, area) == pytest.approx(force)
    assert f.center_of_pressure_depth(yc, area, width * height**3 / 12) == pytest.approx(moment / force)
    # A surface-piercing rectangle: the centre of pressure is two thirds of the way down.
    assert f.center_of_pressure_depth(1.5, 3.0, 3.0**3 / 12) == pytest.approx(2.0)
    with pytest.raises(ValueError):
        f.center_of_pressure_depth(0.0, 1.0, 1.0)


def test_poiseuille_profile_integrates_to_flow_and_peaks_at_twice_the_mean():
    dp, R, L, mu = 50.0, 0.003, 2.0, 1e-3
    Q = integrate.quad(lambda r: 2 * np.pi * r * f.poiseuille_velocity(r, dp, R, L, mu), 0, R)[0]
    assert Q == pytest.approx(f.poiseuille_flow(dp, R, L, mu))
    assert f.poiseuille_velocity(0.0, dp, R, L, mu) == pytest.approx(2 * Q / (np.pi * R**2))
    assert f.poiseuille_velocity(np.array([R]), dp, R, L, mu)[0] == 0
    with pytest.raises(ValueError):
        f.poiseuille_velocity(2 * R, dp, R, L, mu)


def test_smooth_and_fully_rough_limits_of_colebrook():
    Re = np.array([1e4, 1e5, 1e7])
    smooth = f.smooth_pipe_friction_factor(Re)
    assert np.allclose(1 / np.sqrt(smooth), 2 * np.log10(Re * np.sqrt(smooth)) - 0.8)
    assert np.allclose(smooth, f.colebrook_friction_factor(Re, 0.0), rtol=2e-3)
    eps = np.array([1e-3, 0.01, 0.033])
    rough = f.nikuradse_friction_factor(eps)
    assert np.allclose(rough, f.colebrook_friction_factor(1e8, eps), rtol=5e-3)
    assert isinstance(f.smooth_pipe_friction_factor(1e5), float)
    with pytest.raises(ValueError):
        f.smooth_pipe_friction_factor(3000)
    with pytest.raises(ValueError):
        f.nikuradse_friction_factor(0.0)


def test_friction_factors_vectorize_and_reject_transition_anywhere():
    Re = np.array([[1000.0, 1e5], [1e6, 1e7]])
    out = f.colebrook_friction_factor(Re, 1e-4)
    assert out.shape == (2, 2) and out[0, 0] == pytest.approx(0.064)
    assert out[1, 1] == pytest.approx(f.colebrook_friction_factor(1e7, 1e-4))
    with pytest.raises(ValueError, match="transition"):
        f.colebrook_friction_factor([1e5, 3000.0], 1e-4)
    assert np.allclose(f.laminar_friction_factor(np.array([100.0, 1000.0])), [0.64, 0.064])
    with pytest.raises(ValueError):
        f.laminar_friction_factor([1000.0, 2500.0])
    with pytest.raises(ValueError):
        f.euler_number(1.0, 1.0, [1.0, 0.0])
    with pytest.raises(ValueError):
        f.cavitation_number(1.0, 0.0, 1.0, [0.0])


def test_pipe_resistance_reproduces_darcy_weisbach():
    fr, L, D, Q = 0.021, 250.0, 0.15, 0.03
    v = Q / (np.pi * D**2 / 4)
    assert f.pipe_resistance(fr, L, D) * Q**2 == pytest.approx(f.darcy_head_loss(fr, L, D, v))


def _two_loop_network():
    # Junctions A (inflow 0.1), B, C, D (outflows 0.03, 0.04, 0.03).
    # Pipes: 0 A-B, 1 B-D, 2 A-C, 3 C-D, 4 B-C.
    r = np.array([500.0, 800.0, 1200.0, 400.0, 2000.0])
    loops = np.array([[1, 0, -1, 0, 1], [0, 1, 0, -1, -1]])
    Q0 = np.array([0.06, 0.02, 0.04, 0.0, 0.01])
    return r, loops, Q0


def test_hardy_cross_matches_a_newton_solve_and_keeps_continuity():
    r, loops, Q0 = _two_loop_network()
    result = f.hardy_cross(r, loops, Q0)
    incidence = np.array([[-1, 0, -1, 0, 0], [1, -1, 0, 0, -1], [0, 0, 1, -1, 1], [0, 1, 0, 1, 0]])
    assert np.allclose(incidence @ result.flows, incidence @ Q0)
    assert np.max(np.abs(result.loop_residuals)) < 1e-8

    def equations(q):
        return np.concatenate([(incidence @ (q - Q0))[:3], loops @ (r * q * np.abs(q))])

    assert np.allclose(result.flows, optimize.fsolve(equations, Q0, xtol=1e-13), atol=1e-9)
    s = result.summary()
    assert s["iterations"] == result.corrections.size and s["max_loop_residual_m"] < 1e-8
    # Hazen-Williams exponent and a linear (n = 1) network also converge.
    assert f.hardy_cross(r, loops, Q0, n=1.852).corrections[-1] < 1e-10
    lin = f.hardy_cross([1.0, 3.0], [1, -1], [0.2, 0.2], n=1)
    assert np.allclose(lin.flows, [0.3, 0.1])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"loops": [[1, 0, -1, 0]]},
        {"loops": [[2, 0, -1, 0, 1]]},
        {"loops": [[0, 0, 0, 0, 0]]},
        {"initial_flows": [0.1, 0.0]},
        {"initial_flows": [0.0] * 5},
        {"n": 0.5},
        {"tol": 0.0},
        {"max_iter": 0},
        {"resistances": [1.0, -1.0, 1.0, 1.0, 1.0]},
    ],
)
def test_hardy_cross_rejects_bad_input(kwargs):
    r, loops, Q0 = _two_loop_network()
    args = {"resistances": r, "loops": loops, "initial_flows": Q0} | kwargs
    with pytest.raises(ValueError):
        f.hardy_cross(**args)


def test_hardy_cross_reports_nonconvergence():
    r, loops, Q0 = _two_loop_network()
    with pytest.raises(RuntimeError):
        f.hardy_cross(r, loops, Q0, tol=1e-30, max_iter=3)


PIPE = {
    "dp": {"M": 1, "L": -1, "T": -2},
    "rho": {"M": 1, "L": -3},
    "V": {"L": 1, "T": -1},
    "D": {"L": 1},
    "mu": {"M": 1, "L": -1, "T": -1},
    "eps": {"L": 1},
    "L": {"L": 1},
}


def test_buckingham_pi_finds_the_pipe_flow_groups():
    pi = f.buckingham_pi(PIPE, repeating=["rho", "V", "D"])
    assert pi.rank == 3 and pi.exponents.shape == (4, 7)
    assert not np.any(pi.dimension_matrix @ pi.exponents.T)
    assert pi.groups() == ["rho^-1 V^-2 dp^1", "rho^-1 V^-1 D^-1 mu^1", "D^-1 eps^1", "D^-1 L^1"]
    assert pi.summary() == {"variables": 7, "rank": 3, "groups": 4, "pi": pi.groups()}
    assert pi.bases == ("M", "L", "T") and pi.repeating == ("rho", "V", "D")
    # The default picks the first independent variables; V is skipped because dp ~ rho V^2.
    default = f.buckingham_pi(PIPE)
    assert default.repeating == ("dp", "rho", "D")
    assert not np.any(default.dimension_matrix @ default.exponents.T)


def test_buckingham_pi_edge_cases_and_errors():
    alone = f.buckingham_pi({"a": {}, "b": {"X": 0}})
    assert alone.rank == 0 and alone.groups() == ["a^1", "b^1"]
    custom = f.buckingham_pi({"T": {"Θ": 1}, "x": {"L": 1, "Q": 1}, "y": {"L": 1, "Q": 1}})
    assert custom.bases == ("L", "Θ", "Q") and custom.groups() == ["x^-1 y^1"]
    for variables, repeating in [
        ([1, 2], None),
        ({"a": {"L": 1}}, None),
        ({"a": {"L": 0.5}, "b": {}}, None),
        ({"a": {"L": True}, "b": {}}, None),
        ({"a": [1], "b": {}}, None),
        (PIPE, ["rho", "V"]),
        (PIPE, ["D", "L", "rho"]),
        (PIPE, ["rho", "V", "nope"]),
        (PIPE, ["rho", "rho", "V"]),
    ]:
        with pytest.raises(ValueError):
            f.buckingham_pi(variables, repeating)


def test_rational_null_space_is_exact_and_coprime():
    rref, pivots = rational_rref([[2, 4, 6], [1, 2, 4]])
    assert pivots == [0, 2] and rref[0][0] == 1
    basis = integer_null_space([[2, 4, 6], [1, 2, 4]])
    assert basis.tolist() == [[-2, 1, 0]]
    assert integer_null_space([[1, 0], [0, 1]]).shape == (0, 2)


def test_manning_is_chezy_with_a_depth_dependent_coefficient():
    n, S = 0.015, 4e-4
    for depth in [0.5, 1.0, 2.0]:
        area, perim = 3.0 * depth, 3.0 + 2 * depth
        R = f.hydraulic_radius(area, perim)
        assert f.manning_velocity(n, R, S) == pytest.approx(R ** (1 / 6) / n * np.sqrt(R * S))
        assert f.manning_flow(n, area, perim, S) == pytest.approx(area * f.manning_velocity(n, R, S))
    # A full circular pipe: R = D/4, a quarter of the hydraulic diameter.
    D = 0.6
    assert 4 * f.hydraulic_radius(np.pi * D**2 / 4, np.pi * D) == pytest.approx(f.hydraulic_diameter(np.pi * D**2 / 4, np.pi * D))
    assert f.manning_velocity(0.013, 1.0, 0.0) == 0


def test_wave_speed_limits_and_closure_time():
    K, rho = 2.2e9, 998.0
    assert f.pipe_wave_speed(K, rho, 0.5, 0.01, 200e9, c1=0) == pytest.approx(np.sqrt(K / rho))
    steel, pvc = f.pipe_wave_speed(K, rho, 0.3, 0.01, [200e9, 3e9])
    assert 1000 < steel < np.sqrt(K / rho) and pvc < 500
    assert f.critical_closure_time(600.0, 1200.0) == 1.0
    with pytest.raises(ValueError):
        f.pipe_wave_speed(K, rho, 0.3, 0.01, 200e9, c1=-1)


def test_npsh_thoma_and_cavitation_number_agree():
    p, pv, rho, v = 90e3, 2.3e3, 998.0, 3.0
    npsh = f.npsh_available(p, 0.0, pv, rho)
    assert npsh == pytest.approx(f.cavitation_number(p, pv, rho, v) * f.velocity_head(v))
    assert f.npsh_available(p, v, pv, rho) == pytest.approx(npsh + f.velocity_head(v))
    assert f.thoma_sigma(npsh, 40.0) == pytest.approx(npsh / 40)
    with pytest.raises(ValueError):
        f.npsh_available(0.0, 0.0, pv, rho)
    with pytest.raises(ValueError):
        f.thoma_sigma(1.0, 0.0)


def test_venturi_discharge_coefficient():
    ideal = f.venturi_flow(2e4, 998.0, 0.0177, 0.0044)
    assert f.venturi_flow(2e4, 998.0, 0.0177, 0.0044, 0.98) == pytest.approx(0.98 * ideal)
    for cd in [0.0, 1.1]:
        with pytest.raises(ValueError):
            f.venturi_flow(2e4, 998.0, 0.0177, 0.0044, cd)
    with pytest.raises(ValueError):
        f.venturi_flow(2e4, 998.0, [0.02, 0.01], 0.015)

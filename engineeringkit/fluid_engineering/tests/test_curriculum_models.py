"""Analytical checks for the curriculum models of ROADMAP item D5 (White; Munson et al.)."""

import numpy as np
import pytest
from scipy import integrate, optimize

from engineeringkit import fluid_engineering as f

# Momentum equation


def _momentum_balance_on_bend(rho, Q, A1, A2, p1, p2, theta):
    # Independent vector control-volume balance: force of fluid on bend = -(force of bend on fluid).
    e1, e2 = np.array([1.0, 0.0]), np.array([np.cos(theta), np.sin(theta)])
    v1, v2 = Q / A1 * e1, Q / A2 * e2
    on_fluid_from_bend = rho * Q * (v2 - v1) - p1 * A1 * e1 + p2 * A2 * e2
    return -on_fluid_from_bend


@pytest.mark.parametrize("theta", [0.0, np.pi / 4, np.pi / 2, 2.0, np.pi])
def test_bend_force_matches_a_vector_momentum_balance(theta):
    rho, Q, A1, A2, p1 = 998.0, 0.08, 0.0314, 0.0154, 2.0e5
    v1, v2 = Q / A1, Q / A2
    p2 = f.bernoulli_pressure(p1, rho, v1, v2)
    Fx, Fy = f.bend_force(rho, Q, A1, A2, p1, p2, theta)
    assert np.allclose([Fx, Fy], _momentum_balance_on_bend(rho, Q, A1, A2, p1, p2, theta))


def test_bend_force_limits():
    # A straight pipe of constant section with no pressure drop feels no force.
    assert np.allclose(f.bend_force(1000.0, 0.05, 0.02, 0.02, 1e5, 1e5, 0.0), 0.0)
    # The reversal pushes back with both ends' pressure and momentum fluxes.
    Fx, Fy = f.bend_force(1000.0, 0.05, 0.02, 0.02, 1e5, 1e5, np.pi)
    assert Fx == pytest.approx(2 * (1e5 * 0.02 + 1000.0 * 0.05**2 / 0.02)) and Fy == pytest.approx(0.0, abs=1e-9)
    with pytest.raises(ValueError):
        f.bend_force(1000.0, 0.05, 0.02, 0.02, 1e5, 1e5, 4.0)


def test_jet_force_on_plates_and_moving_vanes():
    rho, A, V = 1000.0, 0.002, 30.0
    plate = f.jet_force(rho, A, V, np.pi / 2)
    assert plate[0] == pytest.approx(rho * A * V**2) and plate[1] == pytest.approx(rho * A * V**2)
    reversed_ = f.jet_force(rho, A, V, np.pi)
    assert reversed_[0] == pytest.approx(2 * rho * A * V**2)
    # A vane running at u sees the jet at V - u: momentum flux (rho A (V-u)) (V-u).
    u = 12.0
    Fx, _ = f.jet_force(rho, A, V, np.pi / 2, vane_speed=u)
    assert Fx == pytest.approx(rho * A * (V - u) ** 2)
    assert np.allclose(f.jet_force(rho, A, V, np.array([0.0, np.pi])), [[0.0, 2 * rho * A * V**2], [0.0, 0.0]], atol=1e-9)
    with pytest.raises(ValueError):
        f.jet_force(rho, A, V, np.pi, vane_speed=31.0)


def test_pelton_wheel_power_peaks_at_half_the_jet_speed():
    rho, Q, V, theta = 1000.0, 0.05, 60.0, np.radians(165)
    best = optimize.minimize_scalar(lambda u: -f.pelton_wheel_power(rho, Q, V, u, theta), bounds=(0, V), method="bounded")
    assert best.x == pytest.approx(V / 2, rel=1e-5)
    jet_power = 0.5 * rho * Q * V**2
    assert f.pelton_wheel_power(rho, Q, V, V / 2, theta) / jet_power == pytest.approx((1 - np.cos(theta)) / 2)
    assert f.pelton_wheel_power(rho, Q, V, V / 2, np.pi) == pytest.approx(jet_power)
    # The wheel turns the whole flow: its force is one vane's force times V/(V - u).
    u, A = 20.0, Q / V
    one_vane = f.jet_force(rho, A, V, theta, u)[0]
    assert f.pelton_wheel_power(rho, Q, V, u, theta) == pytest.approx(one_vane * V / (V - u) * u)
    assert f.pelton_wheel_power(rho, Q, V, 0.0) == 0
    with pytest.raises(ValueError):
        f.pelton_wheel_power(rho, Q, V, 61.0)


# Pumps


def test_euler_head_equals_shaft_power_over_weight_flow():
    rho, Q, omega, r1, r2 = 998.0, 0.03, 150.0, 0.05, 0.15
    vt1, vt2 = 2.0, 16.0
    torque = rho * Q * (r2 * vt2 - r1 * vt1)  # angular momentum balance
    H = f.euler_turbomachine_head(omega * r2, vt2, omega * r1, vt1)
    assert rho * f.G * Q * H == pytest.approx(torque * omega)
    assert f.euler_turbomachine_head(20.0, 10.0) == pytest.approx(200.0 / f.G)
    with pytest.raises(ValueError):
        f.euler_turbomachine_head(20.0, 10.0, g=0.0)


@pytest.mark.parametrize("series,parallel", [(1, 1), (2, 1), (1, 2), (2, 3)])
def test_operating_point_lies_on_both_curves(series, parallel):
    H0, a, b, Hs, k = 45.0, 2500.0, 30.0, 15.0, 4000.0
    op = f.pump_operating_point(H0, a, Hs, k, b=b, series=series, parallel=parallel)
    each = f.pump_curve_head(op.flow / parallel, H0, a, b)
    assert series * each == pytest.approx(op.head)
    assert f.system_curve_head(op.flow, Hs, k) == pytest.approx(op.head)
    root = optimize.brentq(lambda Q: series * f.pump_curve_head(Q / parallel, H0, a, b) - f.system_curve_head(Q, Hs, k), 0, 1)
    assert op.flow == pytest.approx(root)
    assert op.pump_head(op.flow) == pytest.approx(op.system_head(op.flow))
    assert op.summary() == {"flow_m3_s": op.flow, "head_m": op.head, "series": series, "parallel": parallel}


def test_series_wins_on_a_steep_system_and_parallel_on_a_flat_one():
    H0, a = 30.0, 3000.0
    steep = [f.pump_operating_point(H0, a, 5.0, 20000.0, series=2).flow, f.pump_operating_point(H0, a, 5.0, 20000.0, parallel=2).flow]
    flat = [f.pump_operating_point(H0, a, 5.0, 200.0, series=2).flow, f.pump_operating_point(H0, a, 5.0, 200.0, parallel=2).flow]
    assert steep[0] > steep[1] and flat[1] > flat[0]


@pytest.mark.parametrize(
    "kwargs",
    [{"static_head": 40.0}, {"a": 0.0, "k": 0.0}, {"series": 0}, {"parallel": 1.5}, {"shutoff_head": 0.0}, {"k": -1.0}],
)
def test_operating_point_rejects_bad_input(kwargs):
    args = {"shutoff_head": 40.0, "a": 1000.0, "static_head": 10.0, "k": 2000.0} | kwargs
    with pytest.raises(ValueError):
        f.pump_operating_point(**args)


def test_affinity_laws_keep_the_dimensionless_coefficients():
    Q, H, P, N, D = 0.05, 32.0, 21e3, 150.0, 0.3
    n, d = np.array([0.5, 1.2, 1.0]), np.array([1.0, 0.9, 1.1])
    Q2, H2, P2 = f.affinity_laws(Q, H, P, n, d)
    N2, D2 = N * n, D * d
    assert np.allclose(Q2 / (N2 * D2**3), Q / (N * D**3))
    assert np.allclose(f.G * H2 / (N2**2 * D2**2), f.G * H / (N**2 * D**2))
    assert np.allclose(P2 / (N2**3 * D2**5), P / (N**3 * D**5))
    # Hydraulic power over shaft power, the efficiency, is unchanged.
    rho = 1000.0
    assert np.allclose(f.pump_hydraulic_power(rho, Q2, H2) / P2, f.pump_hydraulic_power(rho, Q, H) / P)
    with pytest.raises(ValueError):
        f.affinity_laws(Q, H, P, 0.0)


def test_pump_and_system_curves():
    assert f.pump_curve_head(0.0, 40.0, 1000.0) == 40.0
    assert f.pump_curve_head(np.array([0.1, 0.2]), 40.0, 1000.0, 50.0) == pytest.approx([35.0, 10.0])
    assert f.system_curve_head(0.0, -3.0, 500.0) == -3.0
    with pytest.raises(ValueError):
        f.pump_curve_head(-0.1, 40.0, 1000.0)
    with pytest.raises(ValueError):
        f.system_curve_head(0.1, 5.0, -1.0)


# Open channels


def test_specific_energy_is_least_at_the_critical_depth():
    q = 3.0
    yc = f.critical_depth(q)
    least = optimize.minimize_scalar(lambda y: f.specific_energy(y, q), bounds=(0.1, 5.0), method="bounded", options={"xatol": 1e-10})
    assert least.x == pytest.approx(yc, rel=1e-6)
    assert f.specific_energy(yc, q) == pytest.approx(1.5 * yc)
    assert f.froude_number(q / yc, yc) == pytest.approx(1.0)
    assert f.specific_energy(2.0, 0.0) == 2.0


def test_alternate_depths_share_the_energy_on_either_side_of_critical():
    q, E = 2.5, 2.0
    deep, shallow = f.alternate_depths(E, q)
    yc = f.critical_depth(q)
    assert shallow < yc < deep
    assert f.specific_energy(deep, q) == pytest.approx(E) and f.specific_energy(shallow, q) == pytest.approx(E)
    assert f.froude_number(q / deep, deep) < 1 < f.froude_number(q / shallow, shallow)
    # At the minimum energy both depths are critical; below it there is no flow.
    assert f.alternate_depths(1.5 * yc, q) == (pytest.approx(yc), pytest.approx(yc))
    with pytest.raises(ValueError):
        f.alternate_depths(1.4 * yc, q)


def test_hydraulic_jump_conserves_momentum_and_loses_energy():
    y1, v1 = 0.25, 8.0
    q = y1 * v1
    y2 = f.sequent_depth(y1, v1)

    def momentum_function(y):  # pressure force plus momentum flux per unit width, over rho g
        return y * y / 2 + q * q / (f.G * y)

    assert momentum_function(y2) == pytest.approx(momentum_function(y1))
    loss = f.hydraulic_jump_head_loss(y1, y2)
    assert loss == pytest.approx(f.specific_energy(y1, q) - f.specific_energy(y2, q))
    assert loss > 0
    # Bélanger's relation is symmetric: from y2 it returns y1.
    assert f.sequent_depth(y2, q / y2) == pytest.approx(y1)
    # Fr1 = 1 is no jump at all.
    assert f.sequent_depth(1.0, np.sqrt(f.G)) == pytest.approx(1.0)
    assert np.allclose(f.sequent_depth(np.array([0.2, 0.4]), np.array([6.0, 5.0])), [f.sequent_depth(0.2, 6.0), f.sequent_depth(0.4, 5.0)])
    with pytest.raises(ValueError):
        f.hydraulic_jump_head_loss(0.0, 1.0)


# Orifices and weirs


def test_tank_drain_time_matches_an_ode_integration():
    At, Ao, h0, cd = 2.0, 0.003, 3.0, 0.61

    def emptied(t, h):
        return h[0] - 0.75

    emptied.terminal = True
    sol = integrate.solve_ivp(
        lambda t, h: [-f.orifice_flow(Ao, max(h[0], 0.0), cd) / At], (0, 1e4), [h0], events=emptied, rtol=1e-10, atol=1e-12
    )
    assert f.tank_drain_time(At, Ao, h0, 0.75, cd) == pytest.approx(sol.t_events[0][0], rel=1e-6)
    # Half the depth goes in under a third of the time.
    assert f.tank_drain_time(At, Ao, h0, h0 / 2, cd) < f.tank_drain_time(At, Ao, h0, 0.0, cd) / 3
    assert f.orifice_flow(0.01, 2.0, 1.0) == pytest.approx(0.01 * np.sqrt(2 * f.G * 2.0))
    with pytest.raises(ValueError):
        f.tank_drain_time(At, Ao, 1.0, 2.0)
    with pytest.raises(ValueError):
        f.orifice_flow(0.01, 2.0, 0.0)


def test_weir_formulas_integrate_torricelli_over_the_nappe():
    H, b, cd, g = 0.35, 1.8, 0.62, f.G
    rect = integrate.quad(lambda z: b * np.sqrt(2 * g * z), 0, H)[0]
    assert f.rectangular_weir_flow(H, b, cd) == pytest.approx(cd * rect)
    theta = np.radians(60)
    vee = integrate.quad(lambda z: 2 * np.tan(theta / 2) * (H - z) * np.sqrt(2 * g * z), 0, H)[0]
    assert f.v_notch_weir_flow(H, theta, cd) == pytest.approx(cd * vee)
    # Doubling the head: H^(3/2) for the rectangle, H^(5/2) for the notch.
    assert f.rectangular_weir_flow(2 * H, b) / f.rectangular_weir_flow(H, b) == pytest.approx(2**1.5)
    assert f.v_notch_weir_flow(2 * H) / f.v_notch_weir_flow(H) == pytest.approx(2**2.5)
    assert f.rehbock_discharge_coefficient(0.0, 1.0) == 0.611
    assert f.rehbock_discharge_coefficient(0.4, 0.8) == pytest.approx(0.6485)
    with pytest.raises(ValueError):
        f.v_notch_weir_flow(0.1, np.pi)
    with pytest.raises(ValueError):
        f.rectangular_weir_flow(0.1, 1.0, 1.2)
    with pytest.raises(ValueError):
        f.rehbock_discharge_coefficient(0.1, 0.0)


# Drag


def test_stokes_settling_balances_drag_with_weight_less_buoyancy():
    d, rho_p, rho, mu = 50e-6, 2650.0, 998.0, 1.0e-3
    V = f.stokes_settling_velocity(d, rho_p, rho, mu)
    weight_less_buoyancy = (rho_p - rho) * f.G * np.pi * d**3 / 6
    assert f.stokes_drag(mu, V, d) == pytest.approx(weight_less_buoyancy)
    assert f.reynolds_number(rho, V, d, mu) < 1
    assert f.stokes_settling_velocity(1e-3, 50.0, 998.0, 1e-3) < 0  # a light sphere rises
    # The drag coefficient on frontal area is 24/Re.
    Re = 0.05
    v = Re * mu / (rho * d)
    assert f.stokes_drag(mu, v, d) / (0.5 * rho * v**2 * np.pi * d**2 / 4) == pytest.approx(24 / Re)


def test_sphere_drag_curve_limits():
    Re = np.array([1e-4, 1e-3])
    assert np.allclose(f.sphere_drag_coefficient(Re), 24 / Re, rtol=1e-3)
    newton = f.sphere_drag_coefficient(np.array([1e3, 1e4, 1e5]))
    assert np.all((newton > 0.38) & (newton < 0.65))
    assert np.all(np.diff(f.sphere_drag_coefficient(np.geomspace(0.01, 2e5, 50))) < 0)
    with pytest.raises(ValueError):
        f.sphere_drag_coefficient(3e5)
    with pytest.raises(ValueError):
        f.sphere_drag_coefficient(0.0)


def test_terminal_velocity_balances_the_drag_curve_and_reduces_to_stokes():
    rho, mu, rho_p = 998.0, 1.0e-3, 2650.0
    for d in [2e-4, 1e-3, 5e-3]:
        V = f.terminal_velocity(d, rho_p, rho, mu)
        Re = f.reynolds_number(rho, V, d, mu)
        drag = f.sphere_drag_coefficient(Re) * 0.5 * rho * V**2 * np.pi * d**2 / 4
        assert drag == pytest.approx((rho_p - rho) * f.G * np.pi * d**3 / 6, rel=1e-9)
    tiny = 5e-6
    assert f.terminal_velocity(tiny, rho_p, rho, mu) == pytest.approx(f.stokes_settling_velocity(tiny, rho_p, rho, mu), rel=1e-3)
    with pytest.raises(ValueError):
        f.terminal_velocity(1e-3, 900.0, rho, mu)
    with pytest.raises(ValueError, match="drag crisis"):
        f.terminal_velocity(0.5, 7800.0, 1.2, 1.8e-5)


# Couette flow


@pytest.mark.parametrize("dp_dx", [0.0, -2e4, 3e4])
def test_couette_profile_integrates_to_the_flow_and_satisfies_the_momentum_equation(dp_dx):
    h, U, mu = 2e-3, 1.5, 0.08
    q = integrate.quad(lambda y: f.couette_velocity(y, h, U, mu, dp_dx), 0, h)[0]
    assert f.couette_flow_rate(h, U, mu, dp_dx) == pytest.approx(q)
    # mu u'' = dp/dx everywhere, with u(0) = 0 and u(h) = U.
    y, dy = 0.7e-3, 1e-6
    u = [f.couette_velocity(y + s * dy, h, U, mu, dp_dx) for s in (-1, 0, 1)]
    assert mu * (u[0] - 2 * u[1] + u[2]) / dy**2 == pytest.approx(dp_dx, abs=1e-3 * abs(dp_dx) + 1e-3)
    assert f.couette_velocity(0.0, h, U, mu, dp_dx) == 0 and f.couette_velocity(h, h, U, mu, dp_dx) == pytest.approx(U)


def test_couette_shear_stress_and_backflow():
    h, U, mu = 1e-3, 2.0, 0.1
    y = np.linspace(0, h, 5)
    assert np.allclose(np.gradient(f.couette_velocity(y, h, U, mu), y) * mu, f.couette_shear_stress(mu, U, h))
    critical = 2 * mu * U / h**2
    assert f.couette_velocity(1e-5, h, U, mu, 1.5 * critical) < 0 < f.couette_velocity(1e-5, h, U, mu, 0.5 * critical)
    with pytest.raises(ValueError):
        f.couette_velocity(2 * h, h, U, mu)


def test_couette_viscometer_torque_from_the_circular_solution_and_narrow_gap():
    mu, omega, R1, R2, L = 0.9, 3.0, 0.02, 0.025, 0.1
    # Outer cylinder turning: v = A r + B/r with v(R1) = 0, v(R2) = omega R2; tau = -2 mu B / r^2.
    A = omega * R2**2 / (R2**2 - R1**2)
    B = -A * R1**2
    for r in (R1, 0.5 * (R1 + R2), R2):
        torque = abs(-2 * mu * B / r**2) * 2 * np.pi * r * L * r
        assert f.couette_viscometer_torque(mu, omega, R1, R2, L) == pytest.approx(torque)
    narrow = f.couette_viscometer_torque(mu, omega, 0.1, 0.1001, L)
    tau = f.couette_shear_stress(mu, omega * 0.1, 1e-4)
    assert narrow == pytest.approx(tau * 2 * np.pi * 0.1 * L * 0.1, rel=2e-3)
    with pytest.raises(ValueError):
        f.couette_viscometer_torque(mu, omega, R2, R1, L)

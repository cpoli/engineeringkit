"""Analytical checks for the models added against Shigley's and Rao's tables of contents (ROADMAP D8)."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import eigh

from engineeringkit import mechanical_engineering as m

CHAIN_M = np.diag([2.0, 1.0, 1.5])
CHAIN_K = np.array([[3e3, -1e3, 0.0], [-1e3, 3e3, -2e3], [0.0, -2e3, 2e3]])


def test_magnification_and_phase_match_the_simulated_steady_state():
    mass, k, c, F0 = 2.0, 800.0, 8.0, 10.0
    wn, zeta = np.sqrt(k / mass), c / (2 * np.sqrt(k * mass))
    w = 1.3 * wn
    t = np.linspace(0, 20, 20001)
    x = m.simulate_sdof(mass, c, k, t=t, force=lambda s: F0 * np.cos(w * s)).displacement
    tail = t > 15
    X, phi = F0 / k * m.magnification_factor(1.3, zeta), m.phase_angle(1.3, zeta)
    np.testing.assert_allclose(x[tail], X * np.cos(w * t[tail] - phi), atol=1e-6)
    r = np.linspace(0, 3, 30001)
    peak = m.magnification_factor(r, 0.1)
    assert np.isclose(peak.max(), 1 / (2 * 0.1 * np.sqrt(1 - 0.01)), rtol=1e-6)
    assert np.isclose(r[peak.argmax()], np.sqrt(1 - 2 * 0.01), atol=1e-4)
    assert m.magnification_factor(0.0, 0.3) == 1 and np.isinf(m.magnification_factor(1.0, 0.0))
    np.testing.assert_allclose(m.phase_angle(np.array([0.0, 1.0]), [0.05, 0.5]), [0.0, np.pi / 2])
    assert m.phase_angle(100.0, 0.1) > 0.99 * np.pi
    with pytest.raises(ValueError):
        m.magnification_factor(-1.0, 0.1)
    with pytest.raises(ValueError):
        m.phase_angle(1.0, -0.1)


def test_normal_modes_are_mass_and_stiffness_orthogonal():
    modes = m.normal_modes(CHAIN_M, CHAIN_K)
    phi = modes.mode_shapes.T
    np.testing.assert_allclose(phi.T @ CHAIN_M @ phi, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(phi.T @ CHAIN_K @ phi, np.diag(modes.frequencies**2), atol=1e-8)
    np.testing.assert_allclose(modes.frequencies**2, eigh(CHAIN_K, CHAIN_M, eigvals_only=True))
    assert np.all(modes.mode_shapes[np.arange(3), np.argmax(np.abs(modes.mode_shapes), axis=1)] > 0)
    free = m.normal_modes(np.eye(2), [[1.0, -1.0], [-1.0, 1.0]])
    np.testing.assert_allclose(free.frequencies, [0.0, np.sqrt(2)], atol=1e-7)
    assert modes.summary()["modes"] == 3
    for M, K in (
        (np.eye(2), np.eye(3)),
        ([1.0, 2.0], np.eye(2)),
        ([[1.0, 0.5], [0.0, 1.0]], np.eye(2)),
        (-np.eye(2), np.eye(2)),
        (np.eye(2), -np.eye(2)),
        (np.eye(2) * np.nan, np.eye(2)),
    ):
        with pytest.raises(ValueError):
            m.normal_modes(M, K)


def test_modal_superposition_matches_direct_integration():
    x0, v0 = np.array([0.01, -0.02, 0.005]), np.array([0.0, 0.3, -0.1])
    t = np.linspace(0, 1, 201)
    zeta = [0.02, 0.05, 0.1]
    modal = m.modal_response(CHAIN_M, CHAIN_K, x0, v0, t, modal_damping=zeta)
    modes = m.normal_modes(CHAIN_M, CHAIN_K)
    phi, w = modes.mode_shapes.T, modes.frequencies
    C = CHAIN_M @ phi @ np.diag(2 * np.array(zeta) * w) @ phi.T @ CHAIN_M
    from engineeringkit.integrators import integrate_ode

    A = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.linalg.solve(CHAIN_M, CHAIN_K), -np.linalg.solve(CHAIN_M, C)]])
    direct = integrate_ode(lambda _, y: A @ y, np.concatenate((x0, v0)), t, rtol=1e-11, atol=1e-13).history
    np.testing.assert_allclose(modal.displacement, direct[:, :3], atol=1e-9)
    np.testing.assert_allclose(modal.velocity, direct[:, 3:], atol=1e-8)
    rigid = m.modal_response(np.eye(2), [[1.0, -1.0], [-1.0, 1.0]], [0.0, 0.0], [1.0, 1.0], t=[0.0, 2.0], modal_damping=0.1)
    np.testing.assert_allclose(rigid.displacement[-1], [2.0, 2.0], atol=1e-9)
    default = m.modal_response(CHAIN_M, CHAIN_K, x0)
    assert default.displacement.shape == (1001, 3) and default.summary()["coordinates"] == 3
    for kwargs in ({"x0": [1.0]}, {"x0": x0, "v0": [1.0]}, {"x0": x0, "modal_damping": [0.1]}, {"x0": x0, "modal_damping": 1.0}):
        with pytest.raises(ValueError):
            m.modal_response(CHAIN_M, CHAIN_K, **kwargs)


def test_harmonic_response_equals_the_modal_sum():
    F = np.array([1.0, 0.0, -2.0])
    omegas = np.array([5.0, 30.0, 80.0])
    X = m.harmonic_response(CHAIN_M, CHAIN_K, F, omegas)
    modes = m.normal_modes(CHAIN_M, CHAIN_K)
    phi, w = modes.mode_shapes.T, modes.frequencies
    for omega, x in zip(omegas, X):
        np.testing.assert_allclose(x, phi @ ((phi.T @ F) / (w**2 - omega**2)))
    np.testing.assert_allclose(m.harmonic_response(CHAIN_M, CHAIN_K, F, 0.0).real, np.linalg.solve(CHAIN_K, F))
    damped = m.harmonic_response(CHAIN_M, CHAIN_K, F, w[0], C=0.01 * CHAIN_K)
    assert np.all(np.isfinite(damped)) and np.any(np.abs(damped.imag) > 0)
    with pytest.raises(ValueError):
        m.harmonic_response(CHAIN_M, CHAIN_K, F, w[1])
    with pytest.raises(ValueError):
        m.harmonic_response(CHAIN_M, CHAIN_K, F[:2], 1.0)
    with pytest.raises(ValueError):
        m.harmonic_response(CHAIN_M, CHAIN_K, F, -1.0)


def test_rayleigh_and_dunkerley_bracket_the_fundamental():
    flexibility = np.linalg.inv(CHAIN_K)
    masses = np.diag(CHAIN_M)
    exact = m.normal_modes(CHAIN_M, CHAIN_K).frequencies[0]
    static = flexibility @ (masses * 9.81)
    upper = m.rayleigh_frequency(CHAIN_M, CHAIN_K, static)
    lower = m.dunkerley_frequency(masses, flexibility)
    assert lower < exact < upper and upper / exact - 1 < 0.02
    assert np.isclose(1 / lower**2, np.trace(flexibility @ CHAIN_M))
    assert np.isclose(m.rayleigh_frequency(CHAIN_M, CHAIN_K, m.normal_modes(CHAIN_M, CHAIN_K).mode_shapes[0]), exact)
    for shape in ([0.0, 0.0, 0.0], [1.0, 1.0]):
        with pytest.raises(ValueError):
            m.rayleigh_frequency(CHAIN_M, CHAIN_K, shape)
    with pytest.raises(ValueError):
        m.dunkerley_frequency(masses, -flexibility)


def _double_pendulum(m1, m2, l1, l2, g=9.81):
    def L(q, v):
        T = 0.5 * (m1 + m2) * l1**2 * v[0] ** 2 + 0.5 * m2 * l2**2 * v[1] ** 2 + m2 * l1 * l2 * v[0] * v[1] * np.cos(q[0] - q[1])
        return T + (m1 + m2) * g * l1 * np.cos(q[0]) + m2 * g * l2 * np.cos(q[1])

    return L


def test_lagrange_reproduces_the_hand_derived_double_pendulum():
    m1, m2, l1, l2, g = 1.0, 0.5, 1.0, 0.7, 9.81
    L = _double_pendulum(m1, m2, l1, l2, g)
    q, v = np.array([0.6, -0.4]), np.array([1.1, -0.3])
    d = q[0] - q[1]
    H = np.array([[(m1 + m2) * l1**2, m2 * l1 * l2 * np.cos(d)], [m2 * l1 * l2 * np.cos(d), m2 * l2**2]])
    rhs = np.array(
        [
            -m2 * l1 * l2 * v[1] ** 2 * np.sin(d) - (m1 + m2) * g * l1 * np.sin(q[0]),
            m2 * l1 * l2 * v[0] ** 2 * np.sin(d) - m2 * g * l2 * np.sin(q[1]),
        ]
    )
    np.testing.assert_allclose(m.lagrange_accelerations(L, q, v), np.linalg.solve(H, rhs), rtol=1e-6)
    t = np.linspace(0, 3, 301)
    run = m.simulate_lagrangian(L, q, v, t)
    # T + V = L(q, qdot) - 2 L(q, 0), since L(q, 0) = -V.
    energy = [L(qq, vv) - 2 * L(qq, 0 * vv) for qq, vv in zip(run.displacement, run.velocity)]
    np.testing.assert_allclose(energy, energy[0], rtol=1e-6)
    # Linearised about the hanging position, Lagrange's equations give the M-K eigenproblem.
    M_lin = np.array([[(m1 + m2) * l1**2, m2 * l1 * l2], [m2 * l1 * l2, m2 * l2**2]])
    K_lin = np.diag([(m1 + m2) * g * l1, m2 * g * l2])
    w = m.normal_modes(M_lin, K_lin).frequencies
    small = m.simulate_lagrangian(L, 1e-4 * m.normal_modes(M_lin, K_lin).mode_shapes[1], [0.0, 0.0], np.linspace(0, 2 * np.pi / w[1], 3))
    np.testing.assert_allclose(small.displacement[-1], small.displacement[0], atol=1e-8)


def test_lagrange_forces_validation_and_damping():
    L = lambda q, v: 0.5 * 2.0 * v[0] ** 2 - 0.5 * 50.0 * q[0] ** 2  # noqa: E731
    damper = lambda t, q, v: [-4.0 * v[0]]  # noqa: E731
    t = np.linspace(0, 2, 51)
    np.testing.assert_allclose(
        m.simulate_lagrangian(L, [0.1], [0.0], t, forces=damper).displacement[:, 0],
        m.simulate_sdof(2.0, 4.0, 50.0, x0=0.1, t=t).displacement,
        atol=1e-7,
    )
    assert np.isclose(m.lagrange_accelerations(L, [0.1], [0.0], forces=damper)[0], -2.5)
    assert m.simulate_lagrangian(L, [0.1], [0.0]).t.size == 1001
    for bad in ({"forces": lambda t, q, v: [1.0, 2.0]}, {"forces": lambda t, q, v: [np.nan]}, {"step": 0.0}, {"t": np.nan}):
        with pytest.raises(ValueError):
            m.lagrange_accelerations(L, [0.1], [0.0], **bad)
    with pytest.raises(ValueError):
        m.lagrange_accelerations(L, [0.1], [0.0, 1.0])
    with pytest.raises(ValueError):
        m.simulate_lagrangian(L, [0.1], [0.0, 1.0])


def test_lewis_stress_is_the_inscribed_parabola_cantilever():
    # A cantilever of depth l and root thickness t under tip load W has 6Wl/(F t^2);
    # Lewis's x = t^2/(4l) and Y = 2x/(3m) give the same stress.
    W, F, module, t, l = 2000.0, 0.025, 0.004, 0.0065, 0.0072
    Y = 2 * (t * t / (4 * l)) / (3 * module)
    assert np.isclose(m.lewis_bending_stress(W, F, module, Y), 6 * W * l / (F * t * t))
    assert np.isclose(m.lewis_bending_stress(W, F, module, Y, velocity_factor=1.5), 1.5 * 6 * W * l / (F * t * t))
    assert m.lewis_form_factor(12) == 0.245 and m.lewis_form_factor(400) == 0.480
    assert np.isclose(m.lewis_form_factor(23), (0.331 + 0.337) / 2)
    assert np.all(np.diff(m.lewis_form_factor(np.arange(12, 401))) >= 0)
    np.testing.assert_allclose(
        [m.barth_velocity_factor(9.0, p) for p in ("cast", "cut", "hobbed", "ground")],
        [(3.05 + 9) / 3.05, (6.1 + 9) / 6.1, (3.56 + 3) / 3.56, np.sqrt((5.56 + 3) / 5.56)],
    )
    assert m.barth_velocity_factor(0.0, "hobbed") == 1
    for call in (
        lambda: m.lewis_form_factor(11),
        lambda: m.lewis_bending_stress(W, F, module, Y, velocity_factor=0.9),
        lambda: m.lewis_bending_stress(-1.0, F, module, Y),
        lambda: m.barth_velocity_factor(5.0, "planed"),
        lambda: m.barth_velocity_factor(-5.0),
    ):
        with pytest.raises(ValueError):
            call()


def test_shaft_stresses_and_diameters_are_consistent():
    d, M, T = 0.035, 420.0, 610.0
    sigma, tau = 32 * M / (np.pi * d**3), 16 * T / (np.pi * d**3)
    assert np.isclose(m.shaft_von_mises_stress(d, M, T), np.sqrt(sigma**2 + 3 * tau**2))
    assert np.isclose(m.shaft_max_shear_stress(d, M, T), np.hypot(sigma / 2, tau))
    Sy, n = 400e6, 2.0
    d_vm = m.shaft_diameter(M, T, Sy, n)
    d_tr = m.shaft_diameter(M, T, Sy, n, criterion="tresca")
    assert np.isclose(m.shaft_von_mises_stress(d_vm, M, T), Sy / n)
    assert np.isclose(m.shaft_max_shear_stress(d_tr, M, T), Sy / (2 * n))
    assert d_tr > d_vm
    Se, Sut, Kf, Kfs = 180e6, 620e6, 1.6, 1.3
    loads = {"Ma": 300.0, "Mm": 50.0, "Ta": 40.0, "Tm": 500.0}
    dg = m.shaft_diameter_goodman(Se, Sut, safety_factor=n, Kf=Kf, Kfs=Kfs, **loads)
    a = m.shaft_von_mises_stress(dg, Kf * loads["Ma"], Kfs * loads["Ta"])
    mean = m.shaft_von_mises_stress(dg, Kf * loads["Mm"], Kfs * loads["Tm"])
    assert np.isclose(a / Se + mean / Sut, 1 / n)
    for call in (
        lambda: m.shaft_diameter(M, T, Sy, criterion="rankine"),
        lambda: m.shaft_diameter_goodman(Se, Sut),
        lambda: m.shaft_diameter_goodman(Se, Sut, Ma=1.0, Kf=0.5),
        lambda: m.shaft_von_mises_stress(0.0, M, T),
    ):
        with pytest.raises(ValueError):
            call()


def test_spring_rate_from_castigliano_and_wahl_limits():
    d, D, N, G, F = 0.004, 0.032, 8.0, 79.3e9, 150.0
    # Energy of the wire in torsion, U = T^2 L / 2GJ, differentiated by F.
    J, length = np.pi * d**4 / 32, np.pi * D * N
    deflection = F * (D / 2) ** 2 * length / (G * J)
    assert np.isclose(F / m.spring_rate(d, D, N, G), deflection)
    C = m.spring_index(D, d)
    assert C == 8
    nominal = 8 * F * D / (np.pi * d**3)
    assert np.isclose(m.spring_shear_stress(F, D, d, "none"), nominal)
    assert np.isclose(m.spring_shear_stress(F, D, d, "direct"), (1 + 0.5 / C) * nominal)
    assert np.isclose(m.spring_shear_stress(F, D, d), m.wahl_factor(C) * nominal)
    assert np.isclose(m.spring_shear_stress(F, D, d, "bergstrasser"), m.bergstrasser_factor(C) * nominal)
    c = np.linspace(3, 12, 91)
    np.testing.assert_allclose(m.bergstrasser_factor(c), m.wahl_factor(c), rtol=0.02)
    assert np.all(m.bergstrasser_factor(c) < m.wahl_factor(c))
    assert np.all(m.wahl_factor(c) > 1 + 0.5 / c)  # curvature adds to direct shear
    assert np.isclose(m.wahl_factor(1e6), 1, atol=1e-6)
    for call in (
        lambda: m.spring_index(0.004, 0.005),
        lambda: m.wahl_factor(1.0),
        lambda: m.spring_shear_stress(F, D, d, "goehner"),
        lambda: m.spring_rate(d, D, 0, G),
    ):
        with pytest.raises(ValueError):
            call()


def test_member_stiffness_integrates_the_pressure_cone():
    E, d, grip, alpha = 207e9, 0.010, 0.030, np.radians(30)
    D = 1.5 * d

    def area(x):
        outer = D / 2 + x * np.tan(alpha)
        return np.pi * (outer**2 - (d / 2) ** 2)

    compliance = 2 * quad(lambda x: 1 / (E * area(x)), 0, grip / 2)[0]
    assert np.isclose(m.member_stiffness(E, d, grip), 1 / compliance, rtol=1e-10)
    shigley = 0.5774 * np.pi * E * d / (2 * np.log(5 * (0.5774 * grip + 0.5 * d) / (0.5774 * grip + 2.5 * d)))
    assert np.isclose(m.member_stiffness(E, d, grip), shigley, rtol=1e-4)
    assert m.member_stiffness(E, d, grip, washer_diameter=2 * d) > m.member_stiffness(E, d, grip)
    assert np.isclose(m.bolt_stiffness(E, 1e-4, 0.02, 8e-5, 0.01), 1 / (0.02 / (E * 1e-4) + 0.01 / (E * 8e-5)))
    for call in (
        lambda: m.member_stiffness(E, d, grip, washer_diameter=d),
        lambda: m.member_stiffness(E, d, grip, cone_angle=np.pi / 2),
        lambda: m.bolt_stiffness(E, 1e-4, 0.0, 1e-4, 0.0),
    ):
        with pytest.raises(ValueError):
            call()


def test_bolted_joint_shares_load_by_stiffness_until_separation():
    kb, km, Fi, P = 0.8e9, 3.2e9, 25e3, 10e3
    joint = m.bolted_joint(kb, km, Fi, P)
    C = kb / (kb + km)
    # Compatibility: the bolt stretches and the members relax by the same delta.
    delta = P / (kb + km)
    assert np.isclose(joint.bolt_load, Fi + kb * delta) and np.isclose(joint.member_load, Fi - km * delta)
    assert np.isclose(joint.bolt_load - joint.member_load, P)
    assert np.isclose(joint.separation_load, Fi / (1 - C))
    proof = 60e3
    assert np.isclose(joint.load_factor(proof), (proof - Fi) / (C * P))
    near = 28e3  # a proof load reached before the members separate
    assert np.isclose(m.bolted_joint(kb, km, Fi, joint.load_factor(near) * P).bolt_load, near)
    apart = m.bolted_joint(kb, km, Fi, 2 * joint.separation_load)
    assert apart.member_load == 0 and apart.bolt_load == apart.load
    summary = joint.summary()
    assert np.isclose(summary["separation_factor"], joint.separation_load / P) and summary["joint_constant"] == C
    assert np.isclose(m.preload_torque(Fi, 0.012, 0.18), 0.18 * Fi * 0.012)
    with pytest.raises(ValueError):
        joint.load_factor(Fi)
    with pytest.raises(ValueError):
        m.bolted_joint(kb, km, -1.0, P)


def test_two_plane_balance_cancels_force_and_couple():
    masses, radii = [1.2, 0.8, 1.5], [0.10, 0.12, 0.08]
    angles, z = np.radians([0.0, 120.0, 250.0]), [0.10, 0.25, 0.40]
    planes = [0.0, 0.5]
    result = m.balance_rotating_masses(masses, radii, angles, z, planes)
    U = result.mass_radius * np.exp(1j * result.angle)
    u = np.array(masses) * radii * np.exp(1j * angles)
    assert np.isclose(abs(u.sum() + U.sum()), 0, atol=1e-14)
    assert np.isclose(abs((u * z).sum() + (U * planes).sum()), 0, atol=1e-14)
    force, couple = m.unbalance(
        list(masses) + [1.0, 1.0], list(radii) + list(result.mass_radius), list(angles) + list(result.angle), z + planes
    )
    assert force < 1e-14 and couple < 1e-14
    static = m.balance_rotating_masses(masses, radii, angles)
    assert np.isclose(static.mass_radius[0], abs(u.sum())) and static.plane_positions.size == 0
    assert np.isclose(m.unbalance(masses, radii, angles), abs(u.sum()))
    # Static balance alone leaves the couple.
    assert m.unbalance(masses + [1.0], radii + [static.mass_radius[0]], list(angles) + [static.angle[0]], z + [0.0])[1] > 0.01
    assert len(result.summary()["angle_deg"]) == 2
    for call in (
        lambda: m.balance_rotating_masses(masses, radii, angles, z),
        lambda: m.balance_rotating_masses(masses, radii, angles, z, [0.1, 0.1]),
        lambda: m.balance_rotating_masses(masses, radii, angles, z[:2], planes),
        lambda: m.unbalance(masses, radii[:2], angles),
        lambda: m.unbalance(masses, radii, angles, z[:2]),
        lambda: m.unbalance([-1.0], [0.1], [0.0]),
    ):
        with pytest.raises(ValueError):
            call()


def test_fourbar_velocities_and_accelerations_match_finite_differences():
    a, b, c, d = 1.0, 3.0, 2.0, 3.0
    omega2, alpha2, h = 3.0, 2.0, 1e-5
    theta2 = np.linspace(0, 2 * np.pi, 25)
    for crossed in (False, True):
        motion = m.fourbar_kinematics(a, b, c, d, theta2, omega2, alpha2, crossed=crossed)

        def angles(t):
            return m.fourbar_kinematics(a, b, c, d, t, omega2, crossed=crossed)

        plus, minus = angles(theta2 + h), angles(theta2 - h)
        dtheta4 = np.angle(np.exp(1j * (plus.theta4 - minus.theta4))) / (2 * h)
        dtheta3 = np.angle(np.exp(1j * (plus.theta3 - minus.theta3))) / (2 * h)
        np.testing.assert_allclose(motion.omega4, omega2 * dtheta4, atol=1e-6)
        np.testing.assert_allclose(motion.omega3, omega2 * dtheta3, atol=1e-6)
        # alpha = omega2^2 d2theta/dtheta2^2 + alpha2 dtheta/dtheta2
        np.testing.assert_allclose(motion.alpha4, omega2 * (plus.omega4 - minus.omega4) / (2 * h) + alpha2 * dtheta4, atol=1e-4)
        np.testing.assert_allclose(motion.alpha3, omega2 * (plus.omega3 - minus.omega3) / (2 * h) + alpha2 * dtheta3, atol=1e-4)
        closure = a * np.exp(1j * theta2) + b * np.exp(1j * motion.theta3) - d - c * np.exp(1j * motion.theta4)
        np.testing.assert_allclose(np.abs(closure), 0, atol=1e-12)
    single = m.fourbar_kinematics(a, b, c, d, 0.4, omega2)
    assert np.isclose(single.omega4 / omega2, a * np.sin(single.theta3 - 0.4) / (c * np.sin(single.theta3 - single.theta4)))
    assert motion.summary()["peak_omega4_rad_s"] > 0
    with pytest.raises(ValueError):
        m.fourbar_kinematics(a, b, c, d, theta2, np.inf)


def test_slider_crank_kinematics_match_differentiation_and_the_harmonic_series():
    r, l, omega, alpha, h = 0.04, 0.15, 250.0, 30.0, 1e-5
    theta = np.linspace(0, 2 * np.pi, 37)
    motion = m.slider_crank_kinematics(r, l, theta, omega, alpha)
    np.testing.assert_allclose(motion.position, m.slider_crank_displacement(r, l, theta))
    x = lambda t: m.slider_crank_displacement(r, l, t)  # noqa: E731
    dx, d2x = (x(theta + h) - x(theta - h)) / (2 * h), (x(theta + h) - 2 * x(theta) + x(theta - h)) / h**2
    np.testing.assert_allclose(motion.velocity, omega * dx, atol=1e-6)
    np.testing.assert_allclose(motion.acceleration, omega**2 * d2x + alpha * dx, rtol=1e-4, atol=1e-2)
    beta = lambda t: m.slider_crank_kinematics(r, l, t, omega).rod_angle  # noqa: E731
    np.testing.assert_allclose(motion.rod_angular_velocity, omega * (beta(theta + h) - beta(theta - h)) / (2 * h), atol=1e-6)
    np.testing.assert_allclose(np.sin(motion.rod_angle) * l, r * np.sin(theta))
    long_rod = m.slider_crank_kinematics(r, 100 * r, theta, omega)
    n = r / (100 * r)
    np.testing.assert_allclose(long_rod.acceleration, -r * omega**2 * (np.cos(theta) + n * np.cos(2 * theta)), rtol=1e-3, atol=1e-2)
    assert motion.summary()["peak_acceleration_m_s2"] > r * omega**2
    with pytest.raises(ValueError):
        m.slider_crank_kinematics(r, r, theta, omega)


@pytest.mark.parametrize("law", ["uniform", "parabolic", "harmonic", "cycloidal", "polynomial_345", "polynomial_4567"])
def test_cam_laws_rise_the_full_lift_with_consistent_derivatives(law):
    lift, beta = 0.012, np.radians(100)
    theta = np.linspace(0, beta, 2001)
    motion = m.cam_follower_motion(theta, lift, beta, law)
    assert np.isclose(motion.displacement[0], 0, atol=1e-15) and np.isclose(motion.displacement[-1], lift)
    for f, df in ((motion.displacement, motion.velocity), (motion.velocity, motion.acceleration), (motion.acceleration, motion.jerk)):
        interior = np.zeros(theta.size, bool)
        interior[5:-5] = True
        if law == "parabolic":
            interior &= np.abs(theta / beta - 0.5) > 0.01  # the acceleration jumps at mid-rise
        np.testing.assert_allclose(np.gradient(f, theta)[interior], df[interior], rtol=1e-4, atol=1e-4 * np.max(np.abs(df)) + 1e-12)
    if law != "uniform":
        np.testing.assert_allclose(motion.velocity[[0, -1]], 0, atol=1e-15)
    if law in ("cycloidal", "polynomial_345", "polynomial_4567"):
        np.testing.assert_allclose(motion.acceleration[[0, -1]], 0, atol=1e-12)
    assert motion.summary()["law"] == law


def test_cam_peak_accelerations_and_validation():
    lift, beta = 1.0, 1.0
    theta = np.linspace(0, beta, 20001)
    peak = {
        law: m.cam_follower_motion(theta, lift, beta, law).summary()["peak_acceleration_m_rad2"]
        for law in ("parabolic", "harmonic", "cycloidal", "polynomial_345")
    }
    np.testing.assert_allclose([peak["parabolic"], peak["harmonic"], peak["cycloidal"]], [4, np.pi**2 / 2, 2 * np.pi], rtol=1e-6)
    assert np.isclose(peak["polynomial_345"], 10 / np.sqrt(3), rtol=1e-6)
    assert m.cam_follower_motion(0.5, lift, beta, "uniform").velocity == 1
    for call in (lambda: m.cam_follower_motion(1.5, lift, beta), lambda: m.cam_follower_motion(0.5, lift, beta, "dwell")):
        with pytest.raises(ValueError):
            call()

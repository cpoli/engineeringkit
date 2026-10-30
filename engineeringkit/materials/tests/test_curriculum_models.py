"""Analytical and textbook checks for the curriculum models of ROADMAP item D7 (Callister; Dowling)."""

import numpy as np
import pytest

from engineeringkit import materials as m

# Textbook values: Callister and Rethwisch, *Materials Science and Engineering*, chapters 3, 5,
# 6, 9 and 10, and Budynas and Nisbett, *Shigley's Mechanical Engineering Design*, chapter 6.
# Exact example numbers and published answers are the job of ROADMAP item B4.


def test_callister_copper_density_and_iron_diffraction():
    a = m.cubic_lattice_parameter(0.128e-9, "fcc")
    assert np.isclose(m.theoretical_density(4, 63.5e-3, a**3), 8.89e3, rtol=1e-3)
    d = m.interplanar_spacing(0.2866e-9, 2, 2, 0)
    assert np.isclose(d, 0.1013e-9, rtol=1e-3)
    two_theta = np.degrees(2 * m.bragg_angle(d, 0.1790e-9))
    assert np.isclose(two_theta, 124.26, atol=0.2)


def test_callister_fick_flux_diffusivity_and_carburizing_time():
    assert np.isclose(m.fick_flux(3e-11, (0.8 - 1.2) / (10e-3 - 5e-3)), 2.4e-9)
    assert np.isclose(m.diffusion_coefficient(1.2e-4, 130e3, 550 + 273.0), 6.7e-13, rtol=1e-2)
    t = m.diffusion_time(0.5e-3, 0.80, 1.6e-11, 0.25, 1.20)
    assert np.isclose(t, 25400, rtol=5e-3)


def test_callister_lever_rule_and_avrami_half_time():
    W_L, W_alpha = m.lever_rule(35.0, 31.5, 42.5)
    assert np.isclose(W_L, 0.68, atol=5e-3) and np.isclose(W_alpha, 0.32, atol=5e-3)
    n = 3.1
    k = -np.log(1 - 0.30) / (20.0 * 60) ** n
    assert np.isclose(m.avrami_fraction(20.0 * 60, k, n), 0.30)
    assert np.isclose(m.avrami_time(0.5, k, n) / 60, 24.8, atol=0.05)


def test_brinell_number_and_strength_conversion():
    HB = m.brinell_hardness(500 * 9.80665, 10e-3, 1.62e-3)
    assert np.isclose(HB, 241, atol=0.5)
    assert np.isclose(m.brinell_tensile_strength(HB), 3.45e6 * HB)


def test_shigley_marin_factors():
    assert np.isclose(m.marin_surface_factor(520e6, "machined"), 0.860, atol=1e-3)
    assert np.isclose(m.marin_size_factor(32e-3), 0.856, atol=1e-3)
    assert np.isclose(m.marin_reliability_factor(0.99), 0.814, atol=1e-3)
    assert np.isclose(m.marin_reliability_factor(0.5), 1.0)


def test_packing_factors_and_cubic_cells():
    np.testing.assert_allclose(
        [m.atomic_packing_factor(s) for s in ("sc", "bcc", "fcc")], [np.pi / 6, np.pi * np.sqrt(3) / 8, np.pi / np.sqrt(18)]
    )
    assert np.isclose(m.atomic_packing_factor("hcp"), m.atomic_packing_factor("fcc"))
    R = 0.124e-9
    a = m.cubic_lattice_parameter(R, "bcc")
    assert np.isclose(np.sqrt(3) * a, 4 * R)
    with pytest.raises(ValueError):
        m.atomic_packing_factor("bct")


def test_bragg_peaks_follow_the_selection_rules():
    a, lam = 0.3615e-9, 0.1542e-9
    planes, two_theta = m.diffraction_peaks(a, lam, "fcc")
    assert planes[:4].tolist() == [[1, 1, 1], [2, 0, 0], [2, 2, 0], [3, 1, 1]]
    assert np.all(np.diff(two_theta) > 0)
    d = a / np.sqrt(np.sum(planes**2, axis=1))
    np.testing.assert_allclose(2 * d * np.sin(two_theta / 2), lam)
    bcc, _ = m.diffraction_peaks(0.2866e-9, 0.1790e-9, "bcc")
    assert bcc.tolist() == [[1, 1, 0], [2, 0, 0], [2, 1, 1], [2, 2, 0], [3, 1, 0]]
    sc, _ = m.diffraction_peaks(a, lam, "sc", max_index=3)
    assert [1, 0, 0] in sc.tolist() and len({tuple(p) for p in sc}) == len(sc)
    assert np.isclose(m.bragg_angle(1.0, 1.0, order=2), np.pi / 2)
    with pytest.raises(ValueError):
        m.bragg_angle(1.0, 1.0, order=3)
    with pytest.raises(ValueError):
        m.interplanar_spacing(a, 0, 0, 0)
    with pytest.raises(ValueError):
        m.interplanar_spacing(a, 1.5, 0, 0)


def test_error_function_profile_solves_ficks_second_law():
    D, C0, Cs, t = 1e-11, 0.2, 1.0, 3600.0
    x, h = np.linspace(1e-5, 5e-4, 50), 1e-7
    C = m.semi_infinite_concentration(x, t, D, C0, Cs)
    dCdt = (m.semi_infinite_concentration(x, t + 1.0, D, C0, Cs) - m.semi_infinite_concentration(x, t - 1.0, D, C0, Cs)) / 2
    d2C = (m.semi_infinite_concentration(x + h, t, D, C0, Cs) - 2 * C + m.semi_infinite_concentration(x - h, t, D, C0, Cs)) / h**2
    np.testing.assert_allclose(dCdt, D * d2C, rtol=1e-3, atol=1e-12)
    tx = m.diffusion_time(x[10], C[10], D, C0, Cs)
    assert np.isclose(tx, t)
    for bad in ((1e-4, 1.2, D, C0, Cs), (1e-4, 0.5, D, 0.5, 0.5)):
        with pytest.raises(ValueError):
            m.diffusion_time(*bad)


def test_lever_rule_conserves_mass_and_avrami_round_trips():
    W_a, W_b = m.lever_rule(np.array([20.0, 30.0]), 10.0, 40.0)
    np.testing.assert_allclose(W_a * 10 + W_b * 40, [20, 30])
    with pytest.raises(ValueError):
        m.lever_rule(50.0, 10.0, 40.0)
    with pytest.raises(ValueError):
        m.lever_rule(10.0, 10.0, 10.0)
    t = np.geomspace(1.0, 100.0, 9)
    y = m.avrami_fraction(t, 2e-4, 2.5)
    k, n = m.avrami_fit(t, y)
    assert np.isclose(k, 2e-4) and np.isclose(n, 2.5)
    np.testing.assert_allclose(m.avrami_time(y, k, n), t)


def test_vickers_number_and_hardness_inputs():
    assert np.isclose(m.vickers_hardness(9.80665, 1e-3), 1.8544, atol=1e-4)
    with pytest.raises(ValueError):
        m.brinell_hardness(1e3, 1e-3, 2e-3)


def test_mean_stress_lines_are_ordered_and_meet_on_the_axes():
    sa, sm, Se, Sut, Sy = 80e6, np.array([0.0, 100e6, 300e6]), 200e6, 600e6, 450e6
    gerber = m.gerber_safety_factor(sa, sm, Se, Sut)
    goodman = m.goodman_safety_factor(sa, sm, Se, Sut)
    soderberg = m.soderberg_safety_factor(sa, sm, Se, Sy)
    assert np.all(soderberg <= goodman) and np.all(goodman <= gerber)
    assert np.allclose([gerber[0], goodman[0], soderberg[0]], Se / sa)
    np.testing.assert_allclose(gerber * sa / Se + (gerber * sm / Sut) ** 2, 1)
    assert np.isclose(m.gerber_safety_factor(0.0, 100e6, Se, Sut), Sut / 100e6)
    with pytest.raises(ValueError):
        m.soderberg_safety_factor(0.0, 0.0, Se, Sy)


def test_marin_chain_and_endurance_estimate():
    np.testing.assert_allclose(m.endurance_limit_estimate([600e6, 2000e6]), [300e6, 700e6])
    assert m.marin_surface_factor(600e6, "ground") > m.marin_surface_factor(600e6, "forged")
    assert m.marin_surface_factor(600e6, "hot-rolled") < 1
    kb = m.marin_size_factor(np.array([7.62e-3, 100e-3]))
    assert np.isclose(kb[0], 1.0, atol=3e-3) and np.isclose(kb[1], 1.51 * 100**-0.157)
    assert (m.marin_load_factor("bending"), m.marin_load_factor("axial"), m.marin_load_factor("torsion")) == (1.0, 0.85, 0.59)
    assert np.isclose(m.marin_endurance_limit(300e6, ka=0.8, kb=0.9, kc=0.85, ke=0.814), 300e6 * 0.8 * 0.9 * 0.85 * 0.814)
    for bad in (lambda: m.marin_surface_factor(600e6, "polished"), lambda: m.marin_load_factor("shear"), lambda: m.marin_size_factor(1e-3)):
        with pytest.raises(ValueError):
            bad()


def test_notch_factors_lie_between_one_and_kt():
    assert m.elliptical_hole_concentration(1.0, 1.0) == 3
    assert np.isclose(m.elliptical_hole_concentration(4.0, 1.0), 1 + 2 * np.sqrt(4.0 / (1.0 / 4.0)))
    r = np.array([0.1e-3, 1e-3, 10e-3])
    for q in (m.peterson_notch_sensitivity(r, 0.25e-3), m.neuber_notch_sensitivity(r, 0.1e-3)):
        assert np.all((q > 0) & (q < 1)) and np.all(np.diff(q) > 0)
        Kf = m.fatigue_notch_factor(2.5, q)
        assert np.all((Kf >= 1) & (Kf <= 2.5))
    with pytest.raises(ValueError):
        m.fatigue_notch_factor(0.5, 0.5)


def test_neuber_rule_conserves_the_elastic_energy_product():
    E, sy, n, Kt = 200e9, 400e6, 8.0, 3.0
    S = np.array([-300e6, 0.0, 50e6, 300e6])
    stress, strain = m.neuber_notch_stress(S, Kt, E, sy, n)
    np.testing.assert_allclose(stress * strain, (Kt * S) ** 2 / E, rtol=1e-10, atol=1e-30)
    np.testing.assert_allclose(strain, m.ramberg_osgood_strain(stress, E, sy, n))
    assert np.all(np.abs(stress) <= Kt * np.abs(S)) and stress[0] == -stress[-1]
    assert abs(stress[-1]) < sy * 1.2 < Kt * S[-1]


def test_principal_stresses_and_three_dimensional_von_mises():
    p = m.principal_stresses(50.0, -20.0, 10.0, txy=30.0, tyz=-15.0, tzx=5.0)
    assert np.isclose(p.sum(), 40.0) and np.all(np.diff(p) <= 0)
    tensor = np.array([[50.0, 30.0, 5.0], [30.0, -20.0, -15.0], [5.0, -15.0, 10.0]])
    np.testing.assert_allclose(p, np.sort(np.linalg.eigvalsh(tensor))[::-1])
    vm = m.von_mises_stress(50.0, -20.0, 10.0, 30.0, -15.0, 5.0)
    assert np.isclose(vm, np.sqrt(0.5 * ((p[0] - p[1]) ** 2 + (p[1] - p[2]) ** 2 + (p[2] - p[0]) ** 2)))
    assert np.isclose(m.von_mises_stress(80.0, 20.0, 0.0, 40.0), m.von_mises_plane_stress(80.0, 20.0, 40.0))
    assert np.isclose(m.von_mises_stress(150.0, 50.0, 50.0), 100.0)
    batch = m.principal_stresses([1.0, 2.0], 0.0, 0.0)
    assert batch.shape == (2, 3)
    with pytest.raises(ValueError):
        m.von_mises_stress(np.nan, 0.0, 0.0)


def test_halpin_tsai_spans_the_reuss_and_voigt_bounds():
    Ef, Em, Vf = 72e9, 3.5e9, np.linspace(0, 1, 6)
    np.testing.assert_allclose(m.halpin_tsai(Ef, Em, Vf, 0.0), m.rule_of_mixtures_E_transverse(Ef, Em, Vf))
    np.testing.assert_allclose(m.halpin_tsai(Ef, Em, Vf, 1e9), m.rule_of_mixtures_E_longitudinal(Ef, Em, Vf), rtol=1e-6)
    mid = m.halpin_tsai(Ef, Em, Vf, 2.0)
    assert np.all(mid >= m.rule_of_mixtures_E_transverse(Ef, Em, Vf) - 1)


def test_laminate_abd_matches_closed_forms():
    E1, E2, nu12, G12, t = 140e9, 10e9, 0.3, 5e9, 0.125e-3
    Q = m.lamina_stiffness(E1, E2, nu12, G12)
    np.testing.assert_allclose(m.transformed_stiffness(Q, 0.0), Q, atol=1e-3)
    Q90 = m.transformed_stiffness(Q, np.pi / 2)
    assert np.isclose(Q90[0, 0], Q[1, 1]) and np.isclose(Q90[1, 1], Q[0, 0])
    single = m.classical_laminate(E1, E2, nu12, G12, [0.0], 4 * t)
    np.testing.assert_allclose(single.A, Q * 4 * t)
    np.testing.assert_allclose(single.D, Q * (4 * t) ** 3 / 12)
    assert np.isclose(single.effective_moduli()["Ex"], E1) and np.isclose(single.effective_moduli()["nuxy"], nu12)
    quasi = m.classical_laminate(E1, E2, nu12, G12, np.radians([0, 45, -45, 90, 90, -45, 45, 0]), t)
    A = quasi.A
    assert np.isclose(A[0, 0], A[1, 1]) and np.isclose(A[2, 2], (A[0, 0] - A[0, 1]) / 2)
    assert quasi.summary()["symmetric"] and np.isclose(quasi.thickness, 8 * t)
    stress = single.ply_stresses(N=(1e5, 0.0, 0.0))
    np.testing.assert_allclose(stress[0], [1e5 / (4 * t), 0, 0], atol=1.0)
    unsym = m.classical_laminate(E1, E2, nu12, G12, [0.0, np.pi / 2], [t, 2 * t])
    assert not unsym.summary()["symmetric"]
    eps0, kappa = unsym.midplane_response(N=(1e4, 0, 0))
    assert np.any(kappa != 0)
    with pytest.raises(ValueError):
        unsym.midplane_response(N=(1.0, 0.0))
    with pytest.raises(ValueError):
        m.lamina_stiffness(10e9, 140e9, 4.0, 5e9)
    with pytest.raises(ValueError):
        m.transformed_stiffness(np.eye(2), 0.0)


def test_weibull_strength_is_the_reliability_distribution_in_stress():
    from engineeringkit import reliability

    s, s0, mod = np.linspace(0, 500e6, 6), 300e6, 10.0
    np.testing.assert_allclose(m.weibull_failure_probability(s, s0, mod), reliability.weibull_cdf(s, s0, mod))
    V1, V2 = 1e-6, 1e-3
    s2 = m.weibull_size_effect(400e6, V1, V2, mod)
    p1 = m.weibull_failure_probability(400e6, s0, mod, V1, 1e-6)
    p2 = m.weibull_failure_probability(s2, s0, mod, V2, 1e-6)
    assert np.isclose(p1, p2) and s2 < 400e6
    from scipy.special import gamma

    assert np.isclose(m.weibull_mean_strength(s0, mod), s0 * gamma(1 + 1 / mod))

"""Analytical checks for the models added in the heat_transfer rebuild."""

import numpy as np
import pytest
from scipy import integrate, special

from engineeringkit import heat_transfer as h
from engineeringkit.heat_transfer.utils import plane_wall_eigenvalues


def test_plane_wall_eigenvalues_solve_the_characteristic_equation():
    for Bi in [0.01, 1.0, 100.0]:
        z = plane_wall_eigenvalues(Bi, 20)
        assert np.allclose(z * np.tan(z), Bi, rtol=1e-9)
        assert np.all((z > np.arange(20) * np.pi) & (z < np.arange(20) * np.pi + np.pi / 2))
    # Small Bi: zeta_1^2 ~ Bi, which is why the lumped model is exp(-Bi Fo).
    assert np.isclose(plane_wall_eigenvalues(1e-6, 1)[0] ** 2, 1e-6, rtol=1e-6)


def test_plane_wall_transient_limits_and_energy():
    assert h.plane_wall_transient(0.0, 5.0, 0.3, terms=2000) == pytest.approx(1.0, abs=1e-3)
    # Large Bi tends to the fixed-surface-temperature series.
    Fo = 0.3
    fixed = sum(4 / np.pi * (-1) ** n / (2 * n + 1) * np.exp(-(((2 * n + 1) * np.pi / 2) ** 2) * Fo) for n in range(50))
    assert h.plane_wall_transient(Fo, 1e6) == pytest.approx(fixed, rel=1e-5)
    # Heat lost through the surface equals the drop in mean energy.
    Bi, x = 2.0, np.linspace(0, 1, 2001)
    mean = integrate.trapezoid(h.plane_wall_transient(0.4, Bi, x), x)
    lost = integrate.quad(lambda F: Bi * h.plane_wall_transient(F, Bi, 1.0), 0, 0.4, limit=200)[0]
    assert 1 - mean == pytest.approx(lost, rel=1e-4)
    out = h.plane_wall_transient(np.array([0.1, 1.0]), 1.0, np.array([0.0, 1.0]))
    assert out.shape == (2,)


def test_conduction_1d_matches_series_and_conserves_energy_when_insulated():
    L, alpha, n = 0.1, 1e-5, 81
    x = np.linspace(0, L, n)
    T0 = 300 + 100 * np.cos(np.pi * x / L)
    r = h.conduction_1d(T0, L, alpha, 0.5, 400, method="implicit")
    # Adiabatic faces: the mean is conserved (trapezoid weights) and the cosine mode decays.
    assert integrate.trapezoid(r.T[-1], x) == pytest.approx(integrate.trapezoid(T0, x), rel=1e-12)
    expected = 300 + 100 * np.cos(np.pi * x / L) * np.exp(-alpha * (np.pi / L) ** 2 * 200)
    assert np.max(np.abs(r.T[-1] - expected)) < 0.2
    cn = h.conduction_1d([300.0] * 5, 1.0, 1.0, 0.01, 3, T_left=310.0)
    assert cn.T[0, 0] == 310.0 and cn.method == "crank_nicolson"
    summary = cn.summary()
    assert summary["nodes"] == 5 and summary["steps"] == 3 and summary["mesh_fourier"] == pytest.approx(0.16)
    assert summary["final_max_K"] == 310.0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"T_initial": [300.0, 300.0]},
        {"method": "euler"},
        {"steps": 0},
        {"dt": 0.0},
        {"T_left": -1.0},
        {"T_initial": [300.0, np.nan, 300.0]},
    ],
)
def test_conduction_1d_rejects_bad_input(kwargs):
    args = {"T_initial": [300.0] * 5, "length": 1.0, "alpha": 1.0, "dt": 0.1, "steps": 2} | kwargs
    with pytest.raises(ValueError):
        h.conduction_1d(**args)


def test_correlations():
    assert h.dittus_boelter_nusselt(1e4, 1.0) == pytest.approx(0.023 * 1e4**0.8)
    assert h.dittus_boelter_nusselt(1e4, 2.0, heating=False) == pytest.approx(0.023 * 1e4**0.8 * 2**0.3)
    assert h.churchill_chu_nusselt(0.0, 0.7) == pytest.approx(0.825**2)
    assert h.churchill_chu_nusselt(0.0, 0.7, laminar=True) == pytest.approx(0.68)
    assert h.rayleigh_number(1.0, 1.0, 1.0, 1.0, 1.0, g=2.0) == 2.0
    with pytest.raises(ValueError):
        h.dittus_boelter_nusselt(-1, 0.7)
    with pytest.raises(ValueError):
        h.churchill_chu_nusselt(1e6, 0)


def test_planck_integrates_to_stefan_boltzmann_and_has_no_overflow():
    T = 800.0
    total = integrate.quad(h.planck_spectral_emissive_power, 1e-8, 1e-2, args=(T,), points=[h.WIEN_DISPLACEMENT / T], limit=200)[0]
    assert total == pytest.approx(h.stefan_boltzmann_emissive_power(T), rel=1e-6)
    assert h.planck_spectral_emissive_power(1e-9, 10.0) == 0.0
    # Rayleigh-Jeans limit at long wavelength: 2 pi c k T / lambda^4.
    lam = 1.0
    assert h.planck_spectral_emissive_power(lam, 300.0) == pytest.approx(2 * np.pi * 299792458.0 * 1.380649e-23 * 300 / lam**4, rel=1e-4)


def test_view_factors():
    # Very close plates see only each other; far plates see little.
    assert h.view_factor_parallel_plates(1, 1, 1e-4) == pytest.approx(1, abs=1e-3)
    assert h.view_factor_coaxial_disks(1, 1, 1e-6) == pytest.approx(1, abs=1e-5)
    far = h.view_factor_coaxial_disks(0.1, 0.1, 100.0)
    assert far == pytest.approx(0.1**2 / 100**2, rel=1e-3)  # r2^2 / L^2 for small disks
    assert h.view_factor_parallel_plates(2.0, 1.0, 0.5) == pytest.approx(h.view_factor_parallel_plates(1.0, 2.0, 0.5))


def test_neumann_front_and_condensation_and_chf():
    Ste = 0.5
    s = h.neumann_stefan_front(1.0, 1.0, Ste)
    lam = s / 2
    assert lam * np.exp(lam**2) * special.erf(lam) == pytest.approx(Ste / np.sqrt(np.pi))
    assert h.neumann_stefan_front(0.0, 1.0, Ste) == 0.0
    assert h.stefan_number(1.0, 2.0, 4.0) == 0.5
    # h ~ L^(-1/4) and dT^(-1/4)
    base = h.nusselt_film_condensation(960.0, 0.6, 3e-4, 0.68, 2.26e6, 373.0, 363.0, 0.2)
    assert h.nusselt_film_condensation(960.0, 0.6, 3e-4, 0.68, 2.26e6, 373.0, 363.0, 3.2) == pytest.approx(base / 2)
    with pytest.raises(ValueError):
        h.nusselt_film_condensation(960.0, 0.6, 3e-4, 0.68, 2.26e6, 373.0, 380.0, 0.2)
    with pytest.raises(ValueError):
        h.nusselt_film_condensation(1.0, 2.0, 3e-4, 0.68, 2.26e6, 373.0, 363.0, 0.2)
    q = h.zuber_critical_heat_flux(957.9, 0.5955, 2.257e6, 0.0589)
    assert h.zuber_critical_heat_flux(957.9, 0.5955, 2.257e6, 0.0589, C=0.149) == pytest.approx(q * 0.149 / (np.pi / 24))
    with pytest.raises(ValueError):
        h.zuber_critical_heat_flux(1.0, 2.0, 2.257e6, 0.0589)


def test_vectorized_legacy_helpers():
    assert np.allclose(h.fin_efficiency_adiabatic(np.array([0.0, 1.0]), 1.0), [1.0, np.tanh(1.0)])
    r = h.cylinder_resistance(0.01, np.array([0.02, 0.04]), 1.0)
    assert r.shape == (2,)
    with pytest.raises(ValueError):
        h.sphere_resistance(0.01, np.array([0.02, 0.005]), 1.0)


def test_flat_plate_average_is_the_integral_of_the_local_value():
    for Re_L in [2e5, 3e6]:
        Re_c = 5e5
        local = lambda Re: h.flat_plate_nusselt(Re, 0.7, local=True) / Re  # noqa: E731
        points = [Re_c] if Re_L > Re_c else None
        # 0.037 = 0.0296/0.8: the average of the turbulent law, laminar part kept exact.
        integral = integrate.quad(local, 0, Re_L, points=points, limit=200)[0]
        assert h.flat_plate_nusselt(Re_L, 0.7) == pytest.approx(integral, rel=2e-3)
    assert h.flat_plate_nusselt(1e6, 1.0, Re_c=0.0) == pytest.approx(0.037 * 1e6**0.8)
    assert h.flat_plate_nusselt(1e6, 1.0, Re_c=np.inf) == pytest.approx(0.664 * 1e3)
    assert h.flat_plate_nusselt(np.array([1e4, 1e6]), 0.7).shape == (2,)


def test_cylinder_and_tube_bank_correlations():
    assert h.churchill_bernstein_nusselt(0.0, 0.7) == 0.3
    assert h.churchill_bernstein_nusselt(6071.0, 0.7) == pytest.approx(40.6, abs=0.05)
    # Incropera Example 7.7 (staggered bank, 7 rows); the book rounds C to 0.34 and gets 87.9.
    Nu = h.zukauskas_nusselt(13943.0, 0.71, 0.0313, 0.0343, "staggered", rows=7, Pr_s=0.70)
    assert Nu == pytest.approx(87.9, rel=0.015)
    deep = h.zukauskas_nusselt(5e4, 0.7, 0.03, 0.03, rows=25)
    assert deep == pytest.approx(0.27 * 5e4**0.63 * 0.7**0.36)
    assert h.zukauskas_nusselt(5e4, 0.7, 0.03, 0.03, rows=6) == pytest.approx(deep * 0.935)
    assert h.zukauskas_nusselt(5e4, 0.7, 0.07, 0.03, "staggered") == pytest.approx(0.40 * 5e4**0.6 * 0.7**0.36)
    assert h.zukauskas_nusselt(np.array([50.0, 500.0, 5e5]), 0.7, 0.03, 0.03).shape == (3,)
    # Staggered with a tight diagonal: the diagonal gap sets the maximum velocity.
    V = h.tube_bank_max_velocity(1.0, 0.01, 0.03, 0.005, "staggered")
    assert V == pytest.approx(0.03 / (2 * (np.hypot(0.005, 0.015) - 0.01)))
    assert h.tube_bank_max_velocity(1.0, 0.01, 0.03, 0.03, "staggered") == pytest.approx(1.5)


@pytest.mark.parametrize(
    "call",
    [
        lambda: h.flat_plate_nusselt(1e5, 0.7, Re_c=-1.0),
        lambda: h.flat_plate_nusselt(1e5, 0.7, Re_c=np.nan),
        lambda: h.zukauskas_nusselt(5.0, 0.7, 0.03, 0.03),
        lambda: h.zukauskas_nusselt(1e4, 0.7, 0.03, 0.03, "inline"),
        lambda: h.zukauskas_nusselt(1e4, 0.7, 0.03, 0.03, rows=0),
        lambda: h.tube_bank_max_velocity(1.0, 0.01, 0.01, 0.03),
        lambda: h.tube_bank_max_velocity(1.0, 0.01, 0.03, 0.03, "inline"),
        lambda: h.tube_bank_max_velocity(1.0, 0.02, 0.025, 0.001, "staggered"),
        lambda: h.laminar_tube_nusselt("wall"),
        lambda: h.exchanger_effectiveness_crossflow(1.0, 0.5, "both"),
        lambda: h.lmtd_correction_factor(0.8, 1.5),
        lambda: h.lmtd_correction_factor(0.6, 1.5, 2),
        lambda: h.lmtd_correction_factor(1.0, 0.5),
        lambda: h.rohsenow_boiling_heat_flux(2.8e-4, 2.257e6, 1.0, 2.0, 0.0589, 4217.0, 1.76, 10.0),
        lambda: h.shape_factor_buried_sphere(1.0, 0.4),
        lambda: h.shape_factor_buried_cylinder(1.0, 0.5),
        lambda: h.shape_factor_eccentric_cylinders(0.1, 0.2, 0.06),
        lambda: h.shape_factor_parallel_cylinders(0.2, 0.2, 0.1),
        lambda: h.view_factor_reciprocity(0.5, 4.0, 1.0),
        lambda: h.two_surface_radiation(500.0, 300.0, 0.0, 0.5, 1.0, 1.0),
        lambda: h.radiation_shield_heat_flux(500.0, 300.0, 0.5, 0.5, [(0.5, 0.0)]),
        lambda: h.radiation_enclosure([500.0, 300.0], [0.5, 0.5], [1.0, 1.0], [[0.0, 0.9], [0.9, 0.0]]),
        lambda: h.radiation_enclosure([500.0, 300.0], [0.5, 0.5], [1.0, 2.0], [[0.0, 1.0], [1.0, 0.0]]),
        lambda: h.radiation_enclosure([500.0, 300.0], [0.5, 0.5], [1.0], [[0.0, 1.0], [1.0, 0.0]]),
        lambda: h.radiation_enclosure([500.0, 300.0], [0.5, 0.5], [1.0, 1.0], [[0.0, 1.0], [1.0, 0.0]], reradiating=[True]),
        lambda: h.conduction_2d(1.0, 1.0, 2, 5, 300.0, 300.0, 300.0, 400.0),
        lambda: h.conduction_2d(1.0, 1.0, 5, 5, 300.0, 300.0, -1.0, 400.0),
        lambda: h.cylinder_transient(0.1, 1.0, 1.5),
        lambda: h.sphere_transient(-0.1, 1.0),
    ],
)
def test_curriculum_models_reject_bad_input(call):
    with pytest.raises(ValueError):
        call()


def test_internal_flow_and_the_heat_mass_analogy():
    assert h.laminar_tube_nusselt() == pytest.approx(4.364, abs=1e-3)
    assert h.hausen_nusselt(1000.0, 5.0, 0.01, 1e9) == pytest.approx(3.66, abs=1e-6)
    assert h.hausen_nusselt(1000.0, 5.0, 0.01, 0.1) > h.hausen_nusselt(1000.0, 5.0, 0.01, 1.0)
    assert h.entry_length(1000.0, 0.01) == pytest.approx(0.5)
    assert h.entry_length(1e4, 0.01, 7.0) == pytest.approx(0.1)
    assert h.lewis_number(2.0, 4.0) == 0.5
    # Laminar flat plate: Sh from the Sc^(1/3) correlation equals the analogy exactly.
    nu, alpha, D_AB, k, rho, L, Re = 1.5e-5, 2.2e-5, 2.6e-5, 0.026, 1.16, 0.5, 1e5
    cp = k / (rho * alpha)
    h_heat = h.flat_plate_nusselt(Re, nu / alpha) * k / L
    h_mass = h.flat_plate_nusselt(Re, nu / D_AB) * D_AB / L
    assert h.chilton_colburn_mass_transfer(h_heat, rho, cp, h.lewis_number(alpha, D_AB)) == pytest.approx(h_mass, rel=1e-12)


def test_new_exchanger_arrangements():
    for mixed in ["none", "cmax", "cmin"]:
        assert h.exchanger_effectiveness_crossflow(2.0, 0.0, mixed) == pytest.approx(1 - np.exp(-2.0))
        cross = h.exchanger_effectiveness_crossflow(2.0, 0.5, mixed)
        assert h.exchanger_effectiveness_parallel(2.0, 0.5) < cross < h.exchanger_effectiveness_counterflow(2.0, 0.5)
    # With equal capacity rates it does not matter which stream is mixed.
    assert h.exchanger_effectiveness_crossflow(1.3, 1.0, "cmax") == pytest.approx(h.exchanger_effectiveness_crossflow(1.3, 1.0, "cmin"))
    assert h.exchanger_effectiveness_shell_tube(0.0, 0.5) == 0.0
    assert h.exchanger_effectiveness_shell_tube(np.array([0.0, 1.0]), 0.5, 2).shape == (2,)
    assert h.exchanger_effectiveness_shell_tube(3.0, 0.0) == pytest.approx(1 - np.exp(-3.0))
    many = h.exchanger_effectiveness_shell_tube(2.0, 0.7, 60)
    assert many == pytest.approx(h.exchanger_effectiveness_counterflow(2.0, 0.7), rel=1e-3)
    assert h.lmtd_correction_factor(0.0, 2.0) == 1.0
    assert h.lmtd_correction_factor(0.4, 1.5) == pytest.approx(h.lmtd_correction_factor(0.6, 1 / 1.5))
    assert h.lmtd_correction_factor(0.3, 1.0, 2) > h.lmtd_correction_factor(0.3, 1.0)


def test_fin_with_convecting_tip_and_rohsenow():
    args = (25.0, 0.04, 200.0, 5e-5)
    assert h.fin_heat_rate_convective_tip(*args, 0.1, 400.0, 300.0, h_tip=0.0) == pytest.approx(
        h.fin_heat_rate_adiabatic(*args, 0.1, 400.0, 300.0)
    )
    assert h.fin_heat_rate_convective_tip(*args, 0.0, 400.0, 300.0) == pytest.approx(25.0 * 5e-5 * 100.0)
    corrected = h.fin_heat_rate_adiabatic(*args, 0.1 + 5e-5 / 0.04, 400.0, 300.0)
    assert h.fin_heat_rate_convective_tip(*args, 0.1, 400.0, 300.0) == pytest.approx(corrected, rel=1e-5)
    props = (2.79e-4, 2.257e6, 957.9, 0.5955, 0.0589, 4217.0, 1.76)
    q = h.rohsenow_boiling_heat_flux(*props, np.array([10.0, 20.0]))
    assert q[1] / q[0] == pytest.approx(8.0)
    assert h.rohsenow_boiling_heat_flux(*props, 10.0, C_sf=0.0065) == pytest.approx(q[0] * 8)
    assert h.rohsenow_boiling_heat_flux(*props, 0.0) == 0.0


def test_transient_eigenvalues_and_view_factor_reciprocity():
    from scipy import special

    from engineeringkit.heat_transfer.utils import cylinder_eigenvalues, sphere_eigenvalues

    for Bi in [1e-6, 0.01, 1.0, 100.0]:
        z = cylinder_eigenvalues(Bi, 15)
        assert np.allclose(z * special.j1(z), Bi * special.j0(z), atol=1e-10 * max(Bi, 1))
        z = sphere_eigenvalues(Bi, 15)
        assert np.all((z > np.arange(15) * np.pi) & (z < np.arange(1, 16) * np.pi))
        assert np.allclose(z * np.cos(z) + (Bi - 1) * np.sin(z), 0, atol=1e-10 * max(Bi, 1))
    assert sphere_eigenvalues(1e-6, 1)[0] ** 2 == pytest.approx(3e-6, rel=1e-5)
    assert h.cylinder_transient(np.array([0.1, 1.0]), 1.0, np.array([0.0, 1.0])).shape == (2,)
    F21 = h.view_factor_reciprocity(h.view_factor_coaxial_disks(0.5, 1.0, 0.4), np.pi * 0.25, np.pi)
    assert F21 == pytest.approx(h.view_factor_coaxial_disks(1.0, 0.5, 0.4))

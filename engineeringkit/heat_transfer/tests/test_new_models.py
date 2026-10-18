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

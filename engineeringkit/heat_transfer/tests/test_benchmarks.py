import numpy as np
import pytest

from engineeringkit import heat_transfer as h


def test_wall_energy_and_interface_continuity():
    result = h.composite_wall(400, 300, [(0.1, 1), (0.2, 2)], A=2, h_hot=10, h_cold=20)
    resistances = np.array([0.05, 0.05, 0.05, 0.025])
    assert np.isclose(result.heat_rate, 100 / resistances.sum())
    np.testing.assert_allclose(-np.diff(result.interface_temperatures), result.heat_rate * resistances)
    assert np.isclose(result.interface_temperatures[-1], 300)


def test_cylinder_thin_wall_limit():
    radius, thickness, k, L = 1.0, 1e-6, 20.0, 2.0
    planar = h.plane_wall_resistance(thickness, k, 2 * np.pi * radius * L)
    assert np.isclose(h.cylinder_resistance(radius, radius + thickness, k, L), planar, rtol=1e-6)


def test_lumped_cooling_initial_and_asymptotic():
    assert h.lumped_temperature(400, 300, 10, 2, 1000, 1000, 0.01, 0) == 400
    assert np.isclose(h.lumped_temperature(400, 300, 10, 2, 1000, 1000, 0.01, 1e6), 300)
    assert h.radiation_heat_rate(1, 1, 300, 300) == 0
    assert h.radiation_heat_rate(0, 1, 400, 300) == 0


def test_heat_exchanger_singular_and_small_ntu_limits():
    assert h.log_mean_temperature_difference(20, 20) == 20
    assert np.isclose(h.log_mean_temperature_difference(20 + 1e-8, 20), 20 + 5e-9, rtol=1e-12)
    for Cr in [0, 0.5, 1 - 1e-12, 1]:
        assert h.exchanger_effectiveness_counterflow(0, Cr) == 0
        assert np.isclose(h.exchanger_effectiveness_counterflow(1e-12, Cr), 1e-12, rtol=1e-10)
    assert np.isclose(h.exchanger_effectiveness_counterflow(2, 1), 2 / 3)
    assert h.exchanger_effectiveness_counterflow(2, 0.5) > h.exchanger_effectiveness_parallel(2, 0.5)


@pytest.mark.parametrize("kwargs", [{"h_hot": 0}, {"h_cold": -1}, {"A": 0}])
def test_wall_invalid_boundary_data(kwargs):
    with pytest.raises(ValueError):
        h.composite_wall(400, 300, [(0.1, 1)], **kwargs)


def _mason_crossflow(NTU, Cr):
    """Exact effectiveness of crossflow with both streams unmixed, as a series."""
    from math import factorial

    total = 0.0
    for n in range(80):
        a = 1 - np.exp(-NTU) * sum(NTU**m / factorial(m) for m in range(n + 1))
        b = 1 - np.exp(-Cr * NTU) * sum((Cr * NTU) ** m / factorial(m) for m in range(n + 1))
        total += a * b
    return total / (Cr * NTU)


def test_pohlhausen_coefficient_from_the_blasius_profile():
    from scipy.integrate import cumulative_trapezoid, solve_ivp, trapezoid

    blasius = solve_ivp(
        lambda eta, y: [y[1], y[2], -0.5 * y[0] * y[2]], [0, 15], [0, 0, 0.332057336], dense_output=True, rtol=1e-11, atol=1e-12
    )
    eta = np.linspace(0, 15, 30001)
    F = cumulative_trapezoid(blasius.sol(eta)[0], eta, initial=0)
    for Pr in [0.6, 0.7, 1.0, 3.0, 10.0, 50.0]:
        # The energy equation theta'' + Pr f theta'/2 = 0 integrates in closed form.
        exact = 1 / trapezoid(np.exp(-Pr / 2 * F), eta)
        assert h.flat_plate_nusselt(1.0, Pr, local=True) == pytest.approx(exact, rel=0.02)
    assert h.flat_plate_nusselt(1.0, 1.0, local=True) == pytest.approx(1 / trapezoid(np.exp(-F / 2), eta), rel=1e-3)


def test_graetz_eigenvalue_gives_the_uniform_temperature_nusselt_number():
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    def wall_value(lam):
        sol = solve_ivp(lambda r, y: [y[1], -y[1] / r - lam**2 * (1 - r**2) * y[0]], [1e-8, 1], [1, 0], rtol=1e-11, atol=1e-12)
        return sol.y[0, -1]

    lam = brentq(wall_value, 2, 3.5)
    assert h.laminar_tube_nusselt("temperature") == pytest.approx(lam**2 / 2, rel=1e-4)
    # Uniform flux: theta = r^2 - r^4/4 solves the fully developed energy equation
    # with a unit wall gradient; Nu_D = 2 / (theta_wall - theta_bulk).
    from scipy.integrate import trapezoid

    r = np.linspace(0, 1, 200001)
    u, theta = 2 * (1 - r**2), r**2 - r**4 / 4
    mean = trapezoid(u * theta * r, r) / trapezoid(u * r, r)
    assert h.laminar_tube_nusselt("flux") == pytest.approx(2 / (0.75 - mean), rel=1e-8)


def test_crossflow_unmixed_fit_against_the_exact_series():
    for NTU in [0.25, 0.5, 1.0, 2.0, 3.0, 5.0]:
        for Cr in [0.25, 0.5, 0.75, 1.0]:
            tol = 0.04 if NTU < 2 else 0.012
            assert h.exchanger_effectiveness_crossflow(NTU, Cr) == pytest.approx(_mason_crossflow(NTU, Cr), rel=tol)


@pytest.mark.parametrize("Cr", [0.3, 0.6, 1.0])
@pytest.mark.parametrize("shells", [1, 2, 3])
def test_lmtd_correction_factor_agrees_with_effectiveness_ntu(Cr, shells):
    NTU = 1.5
    # Tube stream is C_min: P is the effectiveness, R the capacity ratio; inlets at 1 and 0.
    P = float(h.exchanger_effectiveness_shell_tube(NTU, Cr, shells))
    F = h.lmtd_correction_factor(P, Cr, shells)
    dt_lm = h.log_mean_temperature_difference(1 - P, 1 - Cr * P)
    assert NTU * F * dt_lm == pytest.approx(P, rel=1e-12)


def test_heisler_cylinder_and_sphere_limits_and_energy():
    from scipy import integrate, special

    for transient, shape in ((h.cylinder_transient, 1), (h.sphere_transient, 2)):
        # Lumped limit: V/A = r/(shape + 1), so theta = exp(-(shape + 1) Bi Fo).
        assert transient(1000.0, 1e-4) == pytest.approx(np.exp(-(shape + 1) * 1e-4 * 1000.0), rel=1e-4)
        assert transient(0.0, 5.0, 0.3, terms=2000) == pytest.approx(1.0, abs=1e-3)
        # Heat lost through the surface equals the drop in mean energy (volume weight r^shape).
        Bi, r = 2.0, np.linspace(0, 1, 4001)
        weight = (shape + 1) * r**shape
        mean = integrate.trapezoid(weight * transient(0.3, Bi, r), r)
        lost = integrate.quad(lambda Fo: (shape + 1) * Bi * transient(Fo, Bi, 1.0), 0, 0.3, limit=200)[0]
        assert 1 - mean == pytest.approx(lost, rel=1e-4)
    # Fixed surface temperature (Bi -> infinity) series.
    Fo, j0 = 0.2, special.jn_zeros(0, 60)
    fixed_cyl = np.sum(2 / (j0 * special.j1(j0)) * np.exp(-(j0**2) * Fo))
    fixed_sph = sum(2 * (-1) ** (n + 1) * np.exp(-((n * np.pi) ** 2) * Fo) for n in range(1, 60))
    assert h.cylinder_transient(Fo, 1e7) == pytest.approx(fixed_cyl, rel=1e-5)
    assert h.sphere_transient(Fo, 1e7) == pytest.approx(fixed_sph, rel=1e-5)


def test_conduction_2d_against_the_fourier_series_and_second_order():
    L, W = 2.0, 1.0
    x0, y0 = 1.0, 0.5

    def series(x, y):
        n = np.arange(1, 400, 2)
        return 2 / np.pi * np.sum(2 / n * np.sin(n * np.pi * x / L) * np.sinh(n * np.pi * y / L) / np.sinh(n * np.pi * W / L))

    exact = 300 + 100 * series(x0, y0)
    errors = []
    for nodes in (11, 21, 41):
        r = h.conduction_2d(L, W, 2 * nodes - 1, nodes, 300.0, 300.0, 300.0, 400.0)
        errors.append(abs(r.T[(nodes - 1) // 2, nodes - 1] - exact))
        assert sum(r.heat_rates.values()) == pytest.approx(0.0, abs=1e-9)
    assert errors[-1] < 0.05
    assert errors[0] / errors[1] == pytest.approx(4, rel=0.1)
    assert errors[1] / errors[2] == pytest.approx(4, rel=0.1)


def test_shape_factor_identities():
    # Equal parallel cylinders: the midplane is isothermal, so S is half that of one buried at w/2.
    assert h.shape_factor_parallel_cylinders(0.3, 0.3, 1.2, 2.0) == pytest.approx(h.shape_factor_buried_cylinder(0.3, 0.6, 2.0) / 2)
    # Concentric limit and the isolated sphere.
    S = h.shape_factor_eccentric_cylinders(0.1, 0.3, 0.0, 2.0)
    assert 1 / S == pytest.approx(h.cylinder_resistance(0.05, 0.15, 1.0, 2.0))
    assert h.shape_factor_buried_sphere(1.0, 1e9) == pytest.approx(2 * np.pi, rel=1e-8)
    assert h.shape_factor_buried_cylinder(0.01, 100.0) == pytest.approx(2 * np.pi / np.log(4 * 100 / 0.01), rel=1e-8)


def test_enclosure_network_against_closed_forms():
    sigma = h.STEFAN_BOLTZMANN
    two = h.radiation_enclosure([500.0, 300.0], [0.8, 0.6], [1.0, 4.0], [[0.0, 1.0], [0.25, 0.75]])
    assert two.heat_rate[0] == pytest.approx(h.two_surface_radiation(500.0, 300.0, 0.8, 0.6, 1.0, 4.0))
    # Three surfaces, the third reradiating (Incropera eq. 13.30).
    A, F = np.array([1.0, 1.0, 2.0]), np.array([[0.0, 0.2, 0.8], [0.2, 0.0, 0.8], [0.4, 0.4, 0.2]])
    r = h.radiation_enclosure([1000.0, 500.0, 0.0], [0.7, 0.4, 0.5], A, F, reradiating=[False, False, True])
    network = 0.3 / 0.7 + 1 / (0.2 + 1 / (1 / 0.8 + 1 / 0.8)) + 0.6 / 0.4
    assert r.heat_rate[0] == pytest.approx(sigma * (1000.0**4 - 500.0**4) / network, rel=1e-12)
    assert r.heat_rate[2] == pytest.approx(0.0, abs=1e-9)
    assert r.radiosity[2] == pytest.approx(sigma * r.temperature[2] ** 4)
    # Black surfaces exchange A_i F_ij sigma (T_i^4 - T_j^4) pairwise.
    T = np.array([800.0, 400.0, 600.0])
    black = h.radiation_enclosure(T, [1.0, 1.0, 1.0], A, F)
    pairwise = sigma * np.sum(A[:, None] * F * (T[:, None] ** 4 - T[None, :] ** 4), axis=1)
    assert np.allclose(black.heat_rate, pairwise)
    # N identical shields divide the flux by N + 1.
    q0 = h.radiation_shield_heat_flux(500.0, 300.0, 0.3, 0.3)
    assert h.radiation_shield_heat_flux(500.0, 300.0, 0.3, 0.3, [0.3, (0.3, 0.3)]) == pytest.approx(q0 / 3)
    assert q0 == pytest.approx(h.two_surface_radiation(500.0, 300.0, 0.3, 0.3, 1.0, 1.0))

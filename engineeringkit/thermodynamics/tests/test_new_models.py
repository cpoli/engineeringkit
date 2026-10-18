"""Real fluids, steam, moist air and finite-time engines: reference and analytical checks."""

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

from engineeringkit import thermodynamics as t
from engineeringkit.thermodynamics.systems import steam
from engineeringkit.thermodynamics.utils import if97

CO2_A, CO2_B = 0.3640, 4.267e-5  # Pa m^6/mol^2, m^3/mol
N2_A, N2_B, N2_CP = 0.137, 3.87e-5, 29.1


# IAPWS-IF97 verification tables 5 and 15: T [K], p [MPa], v [m^3/kg], h [kJ/kg], u [kJ/kg], s [kJ/(kg K)], cp [kJ/(kg K)].
@pytest.mark.parametrize(
    "region, T, p, v, h, u, s, cp",
    [
        (1, 300, 3, 0.100215168e-2, 0.115331273e3, 0.112324818e3, 0.392294792, 0.417301218e1),
        (1, 300, 80, 0.971180894e-3, 0.184142828e3, 0.106448356e3, 0.368563852, 0.401008987e1),
        (1, 500, 3, 0.120241800e-2, 0.975542239e3, 0.971934985e3, 0.258041912e1, 0.465580682e1),
        (2, 300, 0.0035, 0.394913866e2, 0.254991145e4, 0.241169160e4, 0.852238967e1, 0.191300162e1),
        (2, 700, 0.0035, 0.923015898e2, 0.333568375e4, 0.301262819e4, 0.101749996e2, 0.208141274e1),
        (2, 700, 30, 0.542946619e-2, 0.263149474e4, 0.246861076e4, 0.517540298e1, 0.103505092e2),
    ],
)
def test_if97_regions_1_and_2_match_the_release_verification_tables(region, T, p, v, h, u, s, cp):
    st = t.steam_properties(p * 1e6, T)
    assert st.region == region and math.isnan(st.x)
    got = (st.v, st.h / 1e3, st.u / 1e3, st.s / 1e3, st.cp / 1e3)
    np.testing.assert_allclose(got, (v, h, u, s, cp), rtol=1e-8)


def test_if97_region_4_and_b23_match_the_release_verification_tables():
    np.testing.assert_allclose(t.saturation_pressure([300, 500, 600]) / 1e6, [0.353658941e-2, 0.263889776e1, 0.123443146e2], rtol=1e-8)
    np.testing.assert_allclose(t.saturation_temperature([0.1e6, 1e6, 10e6]), [0.372755919e3, 0.453035632e3, 0.584149488e3], rtol=1e-8)
    assert np.isclose(if97.b23_pressure_mpa(0.62315e3), 0.165291643e2)
    assert np.isclose(if97.b23_temperature_k(0.165291643e2), 0.62315e3)


def test_steam_gibbs_relations_hold_numerically():
    # (dh/dT)_p = cp, (dh/ds)_p = T and (dh/dp)_s = v follow from the Gibbs function.
    p, T, dT = 5e6, 700.0, 1e-3
    mid, hi, lo = (t.steam_properties(p, T + d) for d in (0, dT, -dT))
    assert np.isclose((hi.h - lo.h) / (2 * dT), mid.cp, rtol=1e-6)
    assert np.isclose((hi.h - lo.h) / (hi.s - lo.s), T, rtol=1e-6)
    up, down = (t.steam_properties(p + d, s=mid.s) for d in (100.0, -100.0))
    assert np.isclose((up.h - down.h) / 200.0, mid.v, rtol=1e-5)
    assert np.isclose(mid.u, mid.h - p * mid.v)


def test_saturation_states_satisfy_clapeyron_and_equal_gibbs_energy():
    T = 450.0
    p = float(t.saturation_pressure(T))
    f, g = t.steam_properties(p, x=0.0), t.steam_properties(p, x=1.0)
    assert f.region == g.region == 4 and math.isnan(f.cp)
    assert np.isclose(f.T, T, rtol=1e-9)
    assert np.isclose(f.h - T * f.s, g.h - T * g.s, rtol=1e-6, atol=10)
    dpdT = (t.saturation_pressure(T + 0.01) - t.saturation_pressure(T - 0.01)) / 0.02
    assert np.isclose(dpdT, (g.h - f.h) / (T * (g.v - f.v)), rtol=1e-4)
    wet = t.steam_properties(p, x=0.25)
    assert np.isclose(wet.h, 0.75 * f.h + 0.25 * g.h) and wet.x == 0.25


def test_steam_inverse_from_entropy_and_enthalpy_round_trips_every_branch():
    cases = [(10e3, 320.0), (10e3, 700.0), (10e6, 500.0), (10e6, 800.0), (20e6, 400.0), (20e6, 900.0), (300.0, 280.0)]
    for p, T in cases:
        st = t.steam_properties(p, T)
        for key in ("s", "h"):
            back = t.steam_properties(p, **{key: getattr(st, key)})
            assert back.region == st.region and np.isclose(back.T, T, rtol=1e-9)
    wet = t.steam_properties(10e3, x=0.8)
    assert np.isclose(t.steam_properties(10e3, s=wet.s).x, 0.8)
    assert np.isclose(t.steam_properties(10e3, h=wet.h).x, 0.8)
    assert t.steam_properties(10e3, T=320.0).summary()["region"] == 2


@pytest.mark.parametrize(
    "kwargs",
    [
        {"p": 10e6, "s": -500.0},  # below 0 °C liquid
        {"p": 10e3, "h": 5e6},  # above 1073.15 K
        {"p": 1e5},
        {"p": 1e5, "T": 400.0, "x": 1.0},
        {"p": 1e5, "T": 200.0},
        {"p": 25e6, "x": 0.5},  # above the critical pressure
        {"p": 1e5, "x": 1.5},
        {"p": [1e5, 2e5], "T": 400.0},
        {"p": 200e6, "T": 400.0},
        {"p": 1e5, "T": [400.0, 500.0]},
        {"p": 1e5, "s": math.nan},
        {"p": -1.0, "T": 400.0},
    ],
)
def test_steam_properties_rejects_states_outside_regions_1_to_4(kwargs):
    p = kwargs.pop("p")
    with pytest.raises(ValueError):
        t.steam_properties(p, **kwargs)


def test_region3_reproduces_the_if97_verification_table():
    # IF97 Table 33: (rho [kg/m^3], T [K]) -> p [MPa], h, u [kJ/kg], s, cp [kJ/(kg K)].
    table = {
        (500.0, 650.0): (25.5837018, 1863.43019, 1812.26279, 4.05427273, 13.8935717),
        (200.0, 650.0): (22.2930643, 2375.12401, 2263.65868, 4.85438792, 44.6579342),
        (500.0, 750.0): (78.3095639, 2258.68845, 2102.06932, 4.46971906, 6.34165359),
    }
    for (rho, T), (p, h, u, s, cp) in table.items():
        props = if97.region3(rho, T)
        got = [props["p"] / 1e6, props["h"] / 1e3, props["u"] / 1e3, props["s"] / 1e3, props["cp"] / 1e3]
        assert np.allclose(got, [p, h, u, s, cp], rtol=1e-8)


@pytest.mark.parametrize(
    ("p", "T", "v"),
    [
        # IAPWS backward equations v(p, T), 2014 revision, Table 5: one point per subregion 3a-3h.
        (50e6, 630.0, 1.470853100e-3),
        (50e6, 710.0, 2.204728587e-3),
        (20e6, 630.0, 1.761696406e-3),
        (26e6, 656.0, 2.245587720e-3),
        (26e6, 661.0, 2.970225962e-3),
        (26e6, 671.0, 5.019029401e-3),
        (23.6e6, 649.0, 2.163198378e-3),
        (30e6, 690.0, 4.656470142e-3),
    ],
)
def test_region3_density_from_pressure_and_temperature_matches_the_backward_equations(p, T, v):
    st = t.steam_properties(p, T)
    assert st.region == 3 and np.isclose(st.v, v, rtol=1e-5)


def test_region3_saturation_joins_regions_1_and_2_and_closes_at_the_critical_point():
    p13 = float(t.saturation_pressure(623.15))
    T, liquid, vapour = steam._saturated(p13 * (1 + 1e-12))
    # IF97 lets neighbouring regions differ by up to 0.05% in v at their boundary.
    assert np.isclose(liquid["v"], if97.region1(T, p13 / 1e6)["v"], rtol=5e-4)
    assert np.isclose(vapour["v"], if97.region2(T, p13 / 1e6)["v"], rtol=5e-4)
    f, g = t.steam_properties(if97.P_CRIT, x=0.0), t.steam_properties(if97.P_CRIT, x=1.0)
    assert np.isclose(1 / f.v, if97.RHO_CRIT, rtol=1e-3) and np.isclose(f.h, g.h, rtol=1e-3)
    # Equal Gibbs energies and Clapeyron's slope in the near-critical dome.
    T = 640.0
    p = float(t.saturation_pressure(T))
    f, g = t.steam_properties(p, x=0.0), t.steam_properties(p, x=1.0)
    assert f.region == 4 and 1 / f.v > 450 > 200 > 1 / g.v
    assert np.isclose(f.h - T * f.s, g.h - T * g.s, rtol=1e-5)
    dpdT = (t.saturation_pressure(T + 0.01) - t.saturation_pressure(T - 0.01)) / 0.02
    assert np.isclose(dpdT, (g.h - f.h) / (T * (g.v - f.v)), rtol=1e-3)


def test_region3_states_round_trip_from_entropy_and_enthalpy():
    # Compressed liquid, subcritical liquid and vapour branches, supercritical fluid.
    for p, T in [(50e6, 630.0), (20e6, 630.0), (20e6, 642.0), (25e6, 650.0), (40e6, 700.0), (100e6, 860.0)]:
        st = t.steam_properties(p, T)
        assert st.region == 3
        for key in ("s", "h"):
            back = t.steam_properties(p, **{key: getattr(st, key)})
            assert back.region == 3 and np.isclose(back.T, T, rtol=1e-9)
    wet = t.steam_properties(20e6, x=0.4)
    assert np.isclose(t.steam_properties(20e6, s=wet.s).x, 0.4)
    assert np.isclose(t.steam_properties(20e6, h=wet.h).x, 0.4)
    # A value on the region 1-3 boundary returns the boundary temperature.
    edge = t.steam_properties(50e6, h=float(steam._region3("liquid")(623.15, 50.0)["h"]))
    assert np.isclose(edge.T, 623.15)
    with pytest.raises(ValueError):
        if97.region3_density(650.0, 1000.0, "liquid")


def test_saturation_functions_reject_points_off_the_curve():
    with pytest.raises(ValueError):
        t.saturation_pressure(700.0)
    with pytest.raises(ValueError):
        t.saturation_temperature(100.0)


def test_van_der_waals_reduces_to_ideal_gas_and_has_a_critical_inflection():
    assert t.van_der_waals_pressure(1, 300, 0.024, 0, 0) == t.ideal_gas_pressure(1, 300, 0.024)
    assert np.isclose(t.van_der_waals_pressure(1, 300, 10.0, CO2_A, CO2_B), t.ideal_gas_pressure(1, 300, 10.0), rtol=1e-4)
    Tc, pc, vc = t.van_der_waals_critical_point(CO2_A, CO2_B)
    p = lambda v: t.van_der_waals_pressure(1.0, Tc, v, CO2_A, CO2_B)  # noqa: E731
    h = 1e-3 * vc
    assert np.isclose(p(vc), pc)
    assert abs(p(vc + h) - p(vc - h)) / (2 * h) < 1e-5 * pc / vc
    assert abs(p(vc + h) - 2 * p(vc) + p(vc - h)) / h**2 < 1e-3 * pc / vc**2
    assert np.isclose(pc * vc / (t.R_UNIVERSAL * Tc), 3 / 8)


def test_van_der_waals_rejects_volumes_inside_the_molecules():
    with pytest.raises(ValueError):
        t.van_der_waals_pressure(1.0, 300.0, 4e-5, CO2_A, CO2_B)
    with pytest.raises(ValueError):
        t.van_der_waals_critical_point(0.0, CO2_B)


def test_joule_thomson_coefficient_matches_an_isenthalpic_vdw_throttle():
    # Molar enthalpy of a vdW gas: H = cv T - 2a/v + RT v/(v - b), cv = cp - R at low pressure.
    T1, p1, dp = 300.0, 1e5, 1e3
    cv = N2_CP - t.R_UNIVERSAL

    def volume(T, p):
        return brentq(lambda v: t.van_der_waals_pressure(1.0, T, v, N2_A, N2_B) - p, 1e-3, 1.0)

    def enthalpy(T, p):
        v = volume(T, p)
        return cv * T - 2 * N2_A / v + t.R_UNIVERSAL * T * v / (v - N2_B)

    H1 = enthalpy(T1, p1)
    T2 = brentq(lambda T: enthalpy(T, p1 - dp) - H1, T1 - 5, T1 + 5)
    mu = (T1 - T2) / dp
    assert np.isclose(mu, t.joule_thomson_coefficient_vdw(T1, N2_A, N2_B, N2_CP), rtol=0.01)
    T_inv = 2 * N2_A / (t.R_UNIVERSAL * N2_B)
    assert np.isclose(t.joule_thomson_coefficient_vdw(T_inv, N2_A, N2_B, N2_CP), 0, atol=1e-20)
    assert t.joule_thomson_coefficient_vdw(2 * T_inv, N2_A, N2_B, N2_CP) < 0


def test_clausius_clapeyron_slope_and_reference_point():
    T = np.array([350.0, 373.15, 400.0])
    p = t.clausius_clapeyron_pressure(T, 373.15, 101325.0, 2.257e6, 461.5)
    assert np.isclose(p[1], 101325.0)
    slope = np.diff(np.log(p)) / np.diff(1 / T)
    np.testing.assert_allclose(slope, -2.257e6 / 461.5)
    # Within 2% of IF97 within 5 K of the reference point.
    assert np.isclose(t.clausius_clapeyron_pressure(368.15, 373.15, 101418.0, 2.257e6, 461.5), t.saturation_pressure(368.15), rtol=0.02)


def test_antoine_water_constants_track_if97_from_1_to_100_celsius():
    A, B, C = 8.07131 + np.log10(133.322), 1730.63, 233.426 - 273.15
    T = np.linspace(274.15, 373.15, 50)
    np.testing.assert_allclose(t.antoine_vapor_pressure(T, A, B, C), t.saturation_pressure(T), rtol=0.01)
    with pytest.raises(ValueError):
        t.antoine_vapor_pressure(30.0, A, B, C)


def test_psychrometric_round_trips_and_dew_point_saturation():
    T, phi = 303.15, np.array([0.2, 0.5, 0.9])
    w = t.humidity_ratio(T, phi)
    np.testing.assert_allclose(t.relative_humidity(T, w), phi)
    Td = t.dew_point(w)
    np.testing.assert_allclose(t.relative_humidity(Td, w), 1.0)
    assert np.all(Td < T)
    # Saturated air at 20 °C and sea level: about 14.7 g/kg (ASHRAE tables 14.75 g/kg).
    assert np.isclose(t.humidity_ratio(293.15, 1.0), 14.75e-3, rtol=0.005)
    assert t.humidity_ratio(293.15, 0.0) == 0
    for bad in (lambda: t.humidity_ratio(T, 1.2), lambda: t.humidity_ratio(373.15, 1.0, 50e3), lambda: t.dew_point(1e-6)):
        with pytest.raises(ValueError):
            bad()


def test_curzon_ahlborn_is_the_maximum_power_efficiency_of_an_endoreversible_engine():
    T_hot, T_cold, K_hot, K_cold = 800.0, 300.0, 2.0, 5.0

    def power(T_hi):
        # Reversible core between T_hi and T_lo; entropy balance K_h(T_h-T_hi)/T_hi = K_c(T_lo-T_c)/T_lo fixes T_lo.
        q_in = K_hot * (T_hot - T_hi)
        s = q_in / T_hi
        T_lo = K_cold * T_cold / (K_cold - s)
        return q_in * (1 - T_lo / T_hi), 1 - T_lo / T_hi

    res = minimize_scalar(lambda T_hi: -power(T_hi)[0], bounds=(T_cold + 1, T_hot - 1), method="bounded", options={"xatol": 1e-9})
    assert np.isclose(power(res.x)[1], t.curzon_ahlborn_efficiency(T_hot, T_cold), rtol=1e-5)
    np.testing.assert_allclose(t.curzon_ahlborn_efficiency([800.0, 400.0], 300.0), 1 - np.sqrt([300 / 800, 300 / 400]))
    assert t.curzon_ahlborn_efficiency(T_hot, T_cold) < t.carnot_efficiency(T_hot, T_cold)
    with pytest.raises(ValueError):
        t.curzon_ahlborn_efficiency(300.0, 300.0)


def test_rankine_cycle_with_if97_matches_textbook_steam_table_values():
    # Cengel & Boles: 10 MPa, 500 °C to 10 kPa gives about 40.2% efficiency.
    h4 = t.steam_properties(10e3, x=0.0)
    h1 = h4.h + h4.v * (10e6 - 10e3)
    turbine_in = t.steam_properties(10e6, 773.15)
    turbine_out = t.steam_properties(10e3, s=turbine_in.s)
    eta = t.rankine_ideal_efficiency(h1, turbine_in.h, turbine_out.h, h4.h)
    assert round(eta, 3) == 0.402 and round(turbine_out.x, 3) == 0.793
    # The pump work v dp agrees with the region 1 isentropic enthalpy rise.
    pumped = t.steam_properties(10e6, s=h4.s)
    assert np.isclose(pumped.h, h1, rtol=1e-3)
    assert np.isclose(quad(lambda p: t.steam_properties(p, s=h4.s).v, 10e3, 10e6)[0], pumped.h - h4.h, rtol=1e-6)


def test_vectorized_otto_efficiency_and_result_container():
    np.testing.assert_allclose(t.otto_efficiency(np.array([4.0, 8.0])), [t.otto_efficiency(4.0), t.otto_efficiency(8.0)])
    st = t.steam_properties(1e5, 400.0)
    assert isinstance(st, t.SteamState) and set(st.summary()) == {"p", "T", "v", "h", "s", "u", "cp", "x", "region"}

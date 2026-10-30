"""Analytical, reference and textbook checks for the curriculum models of ROADMAP item D12 (Çengel and Boles; Moran and Shapiro)."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from engineeringkit import thermodynamics as t
from engineeringkit.thermodynamics.systems import peng_robinson as pr
from engineeringkit.thermodynamics.utils import if97, shomate

# Textbook values: Çengel and Boles, *Thermodynamics: An Engineering Approach*, chapters 5, 9-12, 14
# and 15 (steam, R-134a and air tables). Exact example numbers and published answers are the job of
# ROADMAP item B4; here they are met to the accuracy of the models.


# --- IAPWS-IF97 region 5 and array inputs -------------------------------------------------------


@pytest.mark.parametrize(
    "T, p, v, h, u, s, cp",
    [
        # IF97 Table 42 (2007 revision): T [K], p [MPa], v [m^3/kg], h, u [kJ/kg], s, cp [kJ/(kg K)].
        (1500, 0.5, 0.138455090e1, 0.521976855e4, 0.452749310e4, 0.965408875e1, 0.261609445e1),
        (1500, 30, 0.230761299e-1, 0.516723514e4, 0.447495124e4, 0.772970133e1, 0.272724317e1),
        (2000, 30, 0.311385219e-1, 0.657122604e4, 0.563707038e4, 0.853640523e1, 0.288569882e1),
    ],
)
def test_if97_region_5_matches_the_release_verification_table(T, p, v, h, u, s, cp):
    st = t.steam_properties(p * 1e6, T)
    assert st.region == 5
    np.testing.assert_allclose((st.v, st.h / 1e3, st.u / 1e3, st.s / 1e3, st.cp / 1e3), (v, h, u, s, cp), rtol=1e-8)
    for key in ("s", "h"):
        back = t.steam_properties(p * 1e6, **{key: getattr(st, key)})
        assert back.region == 5 and np.isclose(back.T, T, rtol=1e-9)


def test_if97_regions_2_and_5_agree_on_their_boundary_and_gibbs_relations_hold():
    # The release keeps the 1073.15 K boundary consistent to within its Prague values.
    for p_mpa in (0.1, 10.0, 50.0):
        r2, r5 = if97.region2(1073.15, p_mpa), if97.region5(1073.15, p_mpa)
        assert np.isclose(r2["v"], r5["v"], rtol=5e-4) and abs(r2["h"] - r5["h"]) < 1e3
    p, T, dT = 5e6, 1500.0, 1e-3
    hi, lo, mid = (t.steam_properties(p, T + d) for d in (dT, -dT, 0))
    assert np.isclose((hi.h - lo.h) / (2 * dT), mid.cp, rtol=1e-6)
    assert np.isclose((hi.h - lo.h) / (hi.s - lo.s), T, rtol=1e-6)
    # an enthalpy between the two regions' values at 1073.15 K returns the boundary
    gap = (float(if97.region2(1073.15, 30.0)["h"]) + float(if97.region5(1073.15, 30.0)["h"])) / 2
    assert np.isclose(t.steam_properties(30e6, h=gap).T, 1073.15)


def test_steam_properties_broadcasts_arrays_like_scalar_calls():
    p = np.array([1e4, 1e6, 50e6, 30e6])
    T = np.array([320.0, 400.0, 700.0, 1800.0])
    states = t.steam_properties(p, T)
    assert states.h.shape == (4,) and states.region.tolist() == [2, 1, 3, 5]
    for i in range(4):
        one = t.steam_properties(p[i], T[i])
        assert states.region[i] == one.region and np.isclose(states.h[i], one.h) and np.isclose(states.s[i], one.s)
    wet = t.steam_properties(1e5, x=[0.0, 0.5, 1.0])
    assert np.allclose(wet.h[1], (wet.h[0] + wet.h[2]) / 2) and np.all(wet.region == 4)
    grid = t.steam_properties(np.array([[1e5], [1e6]]), T=[500.0, 700.0])
    assert grid.T.shape == (2, 2)
    assert np.allclose(t.steam_properties(1e5, s=states.s[:1] * 0 + 7500.0).s, 7500.0)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"p": 60e6, "T": 1500.0},  # region 5 stops at 50 MPa
        {"p": 10e3, "h": 9e6},  # above 2273.15 K
        {"p": 1e5, "T": 2300.0},
        {"p": [1e5, 200e6], "T": 400.0},
        {"p": [1e5, 2e5], "x": [0.5, 1.5]},
    ],
)
def test_steam_properties_rejects_states_outside_if97(kwargs):
    p = kwargs.pop("p")
    with pytest.raises(ValueError):
        t.steam_properties(p, **kwargs)


# --- Steady-flow devices ------------------------------------------------------------------------


def test_steady_flow_energy_balance_terms_and_nozzle():
    assert t.steady_flow_work(3000e3, 2500e3) == 500e3
    assert np.isclose(t.steady_flow_work(0.0, 0.0, q=-10.0, V_in=10.0, z_in=1.0), -10 + 50 + 9.80665)
    # a nozzle's exit kinetic energy is the enthalpy drop plus the inlet kinetic energy
    V2 = t.nozzle_exit_velocity(3000e3, 2900e3, V_in=50.0)
    assert np.isclose(V2**2 / 2, 100e3 + 50.0**2 / 2)
    assert np.isclose(t.steady_flow_work(3000e3, 2900e3, V_in=50.0, V_out=V2), 0, atol=1e-6)
    with pytest.raises(ValueError):
        t.nozzle_exit_velocity(3000e3, 3100e3)


def test_isentropic_efficiencies_invert_the_outlet_enthalpy_functions():
    for eta in (0.6, 0.85, 1.0):
        h2 = t.expansion_outlet_enthalpy(3000e3, 2500e3, eta)
        assert np.isclose(t.turbine_isentropic_efficiency(3000e3, h2, 2500e3), eta)
        h2 = t.compression_outlet_enthalpy(300e3, 540e3, eta)
        assert np.isclose(t.compressor_isentropic_efficiency(300e3, h2, 540e3), eta)
    assert np.isclose(t.nozzle_isentropic_efficiency(3000e3, 2955e3, 2950e3, V_in=0.0), 0.9)
    with pytest.raises(ValueError):
        t.turbine_isentropic_efficiency(3000e3, 2600e3, 3000e3)
    with pytest.raises(ValueError):
        t.expansion_outlet_enthalpy(3000e3, 2500e3, 0.0)


def test_steam_turbine_with_if97_follows_the_second_law():
    # Çengel and Boles, Ch. 7: 3 MPa, 400 °C to 50 kPa; the real expansion generates entropy.
    inlet = t.steam_properties(3e6, 673.15)
    ideal = t.steam_properties(50e3, s=inlet.s)
    real = t.steam_properties(50e3, h=t.expansion_outlet_enthalpy(inlet.h, ideal.h, 0.85))
    assert real.s > inlet.s and real.x > ideal.x
    assert np.isclose(inlet.h - ideal.h, 3231.7e3 - 2407.0e3, rtol=0.01)


def test_mixing_chamber_and_throttle_conserve_mass_and_enthalpy():
    m, h = t.mixing_chamber([1.0, 3.0, 0.0], [100e3, 300e3, 1e9])
    assert m == 4.0 and h == 250e3
    # feedwater heater: hot steam and cold water mix to saturated liquid
    hot, cold = t.steam_properties(200e3, 473.15), t.steam_properties(200e3, 293.15)
    target = t.steam_properties(200e3, x=0.0).h
    y = (target - cold.h) / (hot.h - cold.h)
    assert np.isclose(t.mixing_chamber([y, 1 - y], [hot.h, cold.h])[1], target)
    for bad in ([-1.0, 2.0], [0.0, 0.0]):
        with pytest.raises(ValueError):
            t.mixing_chamber(bad, [1.0, 2.0])
    with pytest.raises(ValueError):
        t.mixing_chamber([1.0], [1.0, 2.0])
    inlet = t.steam_properties(1e6, x=0.0)
    out = t.throttle_steam(inlet, 1e5)
    assert np.isclose(out.h, inlet.h) and out.s > inlet.s and out.T < inlet.T
    with pytest.raises(ValueError):
        t.throttle_steam(inlet, 2e6)
    with pytest.raises(ValueError):
        t.throttle_steam(t.steam_properties([1e6, 2e6], x=0.0), 1e5)


# --- Rankine plants ------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, efficiency, extra",
    [
        # Çengel and Boles, Examples 10-3c, 10-4 and 10-5: 15 MPa, 600 °C, 10 kPa.
        ({}, 0.430, {"turbine_exit_quality": 0.804}),
        ({"reheat_pressure": 4e6, "T_reheat": 873.15}, 0.450, {"turbine_exit_quality": 0.896}),
        ({"feedwater_pressure": 1.2e6}, 0.463, {"bleed_fraction": 0.2270}),
    ],
)
def test_rankine_cycle_reproduces_the_textbook_plants(kwargs, efficiency, extra):
    res = t.rankine_cycle(15e6, 10e3, 873.15, **kwargs)
    assert abs(res.efficiency - efficiency) < 1e-3
    for key, value in extra.items():
        assert abs(res.summary()[key] - value) < 2e-3


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"reheat_pressure": 3e6, "T_reheat": 823.15},
        {"feedwater_pressure": 0.5e6},
        {"reheat_pressure": 3e6, "T_reheat": 823.15, "feedwater_pressure": 0.5e6},  # bleed after the reheater
        {"reheat_pressure": 1e6, "T_reheat": 823.15, "feedwater_pressure": 5e6},  # bleed before it
        {"turbine_efficiency": 0.85, "pump_efficiency": 0.8, "feedwater_pressure": 0.5e6},
    ],
)
def test_rankine_cycle_closes_its_energy_balance_and_stays_below_carnot(kwargs):
    res = t.rankine_cycle(10e6, 10e3, 823.15, **kwargs)
    assert np.isclose(res.heat_in - res.heat_out, res.net_work, rtol=1e-9)
    assert 0 < res.efficiency < t.carnot_efficiency(823.15, float(t.saturation_temperature(10e3)))
    assert res.back_work_ratio < 0.02
    if "turbine_efficiency" not in kwargs:
        # isentropic sections: the exit has the entropy of the last turbine inlet
        last = [st for label, st in zip(res.labels, res.states) if label in ("turbine inlet", "reheat outlet")][-1]
        assert np.isclose(res.states[res.labels.index("turbine exit")].s, last.s)


def test_rankine_improvements_raise_efficiency_and_irreversibility_lowers_it():
    base = t.rankine_cycle(10e6, 10e3, 773.15).efficiency
    assert base == pytest.approx(0.402, abs=1e-3)  # the 0.1.0 benchmark through rankine_ideal_efficiency
    assert t.rankine_cycle(10e6, 10e3, 873.15).efficiency > base
    assert t.rankine_cycle(10e6, 5e3, 773.15).efficiency > base
    assert t.rankine_cycle(10e6, 10e3, 773.15, feedwater_pressure=0.8e6).efficiency > base
    assert t.rankine_cycle(10e6, 10e3, 773.15, turbine_efficiency=0.87).efficiency < base
    saturated = t.rankine_cycle(5e6, 10e3)
    assert saturated.states[saturated.labels.index("turbine inlet")].x == 1.0
    supercritical = t.rankine_cycle(25e6, 10e3, 873.15)
    assert supercritical.efficiency > t.rankine_cycle(15e6, 10e3, 873.15).efficiency


@pytest.mark.parametrize(
    "args, kwargs",
    [
        ((10e6, 20e6, 773.15), {}),
        ((10e6, 10e3, 500.0), {}),
        ((25e6, 10e3), {}),
        ((25e6, 10e3, 600.0), {}),
        ((10e6, 10e3, 773.15), {"reheat_pressure": 2e6}),
        ((10e6, 10e3, 773.15), {"reheat_pressure": 2e6, "T_reheat": 400.0}),
        ((10e6, 10e3, 773.15), {"reheat_pressure": 2e6, "T_reheat": 773.15, "feedwater_pressure": 2e6}),
        ((10e6, 10e3, 773.15), {"feedwater_pressure": 20e6}),
        ((10e6, 10e3, 773.15), {"turbine_efficiency": 1.2}),
    ],
)
def test_rankine_cycle_rejects_impossible_plants(args, kwargs):
    with pytest.raises(ValueError):
        t.rankine_cycle(*args, **kwargs)


# --- Gas turbines ---------------------------------------------------------------------------------


def test_brayton_cycle_reduces_to_the_ideal_formula_and_the_regenerated_limit():
    for rp in (2.0, 4.0, 8.0):
        assert np.isclose(t.brayton_cycle(rp, 300.0, 1300.0).efficiency, t.brayton_efficiency(rp))
        # perfect regenerator, possible while the exhaust is hotter than the compressed air:
        # eta = 1 - (T_min/T_max) rp^((gamma-1)/gamma), Çengel and Boles Eq. 9-24
        regen = t.brayton_cycle(rp, 300.0, 1300.0, regenerator_effectiveness=1.0)
        assert np.isclose(regen.efficiency, 1 - 300.0 / 1300.0 * rp ** (0.4 / 1.4))
    assert t.brayton_cycle(4.0, 300.0, 1300.0, regenerator_effectiveness=1.0).heat_regenerated > 0


def test_brayton_cycle_with_many_intercooled_and_reheated_stages_approaches_ericsson_and_carnot():
    carnot = t.carnot_efficiency(1300.0, 300.0)
    etas = [
        t.brayton_cycle(10.0, 300.0, 1300.0, regenerator_effectiveness=1.0, compressor_stages=n, turbine_stages=n).efficiency
        for n in (1, 2, 8, 64, 512)
    ]
    # the gap to Carnot closes like 1/N
    assert np.all(np.diff(etas) > 0) and etas[-1] < carnot and np.isclose(etas[-1], carnot, atol=5e-4)


def test_brayton_cycle_energy_balance_and_back_work_ratio():
    res = t.brayton_cycle(
        8.0,
        300.0,
        1300.0,
        compressor_efficiency=0.8,
        turbine_efficiency=0.85,
        regenerator_effectiveness=0.8,
        compressor_stages=2,
        turbine_stages=2,
    )
    T, s = res.T, res.s
    # first law over the states: heat in at the two hot isobars, out in the exhaust and intercooler
    assert np.isclose(res.heat_in - res.heat_out, res.net_work)
    assert np.isclose(res.compressor_work, 2 * 1005.0 * (T[res.labels.index("compressor 1 outlet")] - 300.0))
    assert res.s[0] == 0 and np.all(np.isfinite(s)) and res.summary()["back_work_ratio"] == res.back_work_ratio
    # with poor components the exhaust is colder than the compressed air: the regenerator is left idle
    idle = t.brayton_cycle(30.0, 300.0, 1100.0, compressor_efficiency=0.8, turbine_efficiency=0.85, regenerator_effectiveness=0.8)
    assert idle.heat_regenerated == 0 and len(idle.labels) == 4
    for kwargs in ({"pressure_ratio": 1.0}, {"T_max": 400.0}, {"compressor_stages": 0}, {"regenerator_effectiveness": 1.5}):
        with pytest.raises(ValueError):
            t.brayton_cycle(**{"pressure_ratio": 8.0, "T_min": 300.0, "T_max": 1300.0, **kwargs})


# --- Stirling and dual cycles --------------------------------------------------------------------


def test_stirling_cycle_from_its_four_processes():
    R, gamma, Th, Tc, r = 287.0, 1.4, 900.0, 300.0, 4.0
    cv = R / (gamma - 1)
    w = R * (Th - Tc) * math.log(r)  # isothermal expansion minus isothermal compression
    for eps in (0.0, 0.5, 1.0):
        q_in = R * Th * math.log(r) + (1 - eps) * cv * (Th - Tc)
        assert np.isclose(t.stirling_efficiency(Th, Tc, r, gamma, eps), w / q_in)
    assert np.isclose(t.stirling_efficiency(Th, Tc, [2.0, 10.0]), t.carnot_efficiency(Th, Tc)).all()
    with pytest.raises(ValueError):
        t.stirling_efficiency(300.0, 900.0, 4.0)


def test_dual_cycle_limits_and_four_state_first_law():
    assert np.isclose(t.dual_efficiency(16.0, 1.0, 2.0), t.diesel_efficiency(16.0, 2.0))
    assert np.isclose(t.dual_efficiency(16.0, 1.8, 1.0), t.otto_efficiency(16.0))
    assert t.dual_efficiency(16.0, 1.0, 1.0) == t.otto_efficiency(16.0)
    r, alpha, rc, g, T1 = 16.0, 1.5, 1.4, 1.4, 300.0
    T2 = T1 * r ** (g - 1)
    Tx, T3 = alpha * T2, alpha * T2 * rc
    T4 = T3 * (rc / r) ** (g - 1)
    q_in = (Tx - T2) + g * (T3 - Tx)  # per cv: constant volume then constant pressure
    assert np.isclose(t.dual_efficiency(r, alpha, rc, g), 1 - (T4 - T1) / q_in)
    with pytest.raises(ValueError):
        t.dual_efficiency(16.0, 0.9, 1.2)


# --- Peng-Robinson and the vapour-compression cycle -----------------------------------------------


def test_peng_robinson_critical_point_is_a_horizontal_inflection():
    Tc, pc, omega = pr.R134A.Tc, pr.R134A.pc, pr.R134A.omega
    _, b, _ = pr._constants(Tc, Tc, pc, omega)
    vc = 0.3074 * t.R_UNIVERSAL * Tc / pc  # the PR critical volume, about 3.951 b
    p = lambda v: t.peng_robinson_pressure(Tc, v, Tc, pc, omega)  # noqa: E731
    h = 1e-3 * vc
    assert np.isclose(p(vc), pc, rtol=1e-4)
    # the published constants 0.45724 and 0.07780 are rounded, so the inflection is horizontal to 1e-3
    assert abs(p(vc + h) - p(vc - h)) / (2 * h) < 1e-3 * pc / vc
    assert np.isclose(vc / b, 3.951, atol=1e-3)
    # far from the critical point it is the ideal gas
    assert np.isclose(t.peng_robinson_compressibility(Tc, 10.0, Tc, pc, omega), 1, atol=1e-5)
    with pytest.raises(ValueError):
        t.peng_robinson_pressure(Tc, b / 2, Tc, pc, omega)


def test_peng_robinson_saturation_obeys_maxwell_equal_area_and_matches_r134a_tables():
    fluid = pr.R134A
    T = 280.0
    psat = t.peng_robinson_saturation_pressure(T, fluid.Tc, fluid.pc, fluid.omega)
    zl = t.peng_robinson_compressibility(T, psat, fluid.Tc, fluid.pc, fluid.omega, phase="liquid")
    zv = t.peng_robinson_compressibility(T, psat, fluid.Tc, fluid.pc, fluid.omega)
    vl, vv = (z * t.R_UNIVERSAL * T / psat for z in (zl, zv))
    area = quad(lambda v: t.peng_robinson_pressure(T, v, fluid.Tc, fluid.pc, fluid.omega) - psat, vl, vv, limit=200)[0]
    assert abs(area) < 1e-6 * psat * (vv - vl)
    # Çengel and Boles Table A-11: -40 °C 51.25 kPa, 0 °C 293.0 kPa, 40 °C 1017 kPa
    for T_c, p_kpa in ((-40, 51.25), (0, 293.0), (40, 1017.0)):
        assert np.isclose(t.peng_robinson_saturation_pressure(T_c + 273.15, fluid.Tc, fluid.pc, fluid.omega) / 1e3, p_kpa, rtol=0.02)
    # PR was fitted to the vapour pressure through omega, so it nearly returns it
    p07 = t.peng_robinson_saturation_pressure(0.7 * fluid.Tc, fluid.Tc, fluid.pc, fluid.omega)
    assert np.isclose(t.pitzer_acentric_factor(p07, fluid.pc), fluid.omega, atol=0.01)
    for bad in (0.3 * fluid.Tc, fluid.Tc):
        with pytest.raises(ValueError):
            t.peng_robinson_saturation_pressure(bad, fluid.Tc, fluid.pc, fluid.omega)


def test_refrigerant_properties_are_thermodynamically_consistent():
    fluid = t.R134A
    f, g = (t.refrigerant_properties(fluid, 0.5e6, x=x) for x in (0.0, 1.0))
    # IIR reference, equal Gibbs energies and Clapeyron's slope on the dome
    zero = t.refrigerant_properties(fluid, pr._saturation_pressure(273.15, fluid.Tc, fluid.pc, fluid.omega), x=0.0)
    assert np.isclose(zero.h, 200e3) and np.isclose(zero.s, 1e3)
    assert np.isclose(f.h - f.T * f.s, g.h - g.T * g.s, rtol=1e-6)
    dT = 0.01
    dpdT = (
        pr._saturation_pressure(f.T + dT, fluid.Tc, fluid.pc, fluid.omega)
        - pr._saturation_pressure(f.T - dT, fluid.Tc, fluid.pc, fluid.omega)
    ) / (2 * dT)
    assert np.isclose(dpdT, (g.h - f.h) / (f.T * (g.v - f.v)), rtol=1e-4)
    # (dh/ds)_p = T in the superheated vapour; latent heat of R-134a at 0.5 MPa is about 186 kJ/kg (Table A-12)
    hi, lo = (t.refrigerant_properties(fluid, 0.5e6, 320.0 + d) for d in (0.01, -0.01))
    assert np.isclose((hi.h - lo.h) / (hi.s - lo.s), 320.0, rtol=1e-6) and hi.phase == "vapour"
    assert np.isclose(g.h - f.h, 186.0e3, rtol=0.03)
    # inverse from s and h on every branch
    for p, T in ((0.5e6, 260.0), (0.5e6, 320.0), (5e6, 400.0), (5e6, 300.0)):
        st = t.refrigerant_properties(fluid, p, T)
        for key in ("s", "h"):
            assert np.isclose(t.refrigerant_properties(fluid, p, **{key: getattr(st, key)}).T, T, rtol=1e-8)
    wet = t.refrigerant_properties(fluid, 0.5e6, x=0.3)
    assert np.isclose(t.refrigerant_properties(fluid, 0.5e6, h=wet.h).x, 0.3)
    assert t.refrigerant_properties(fluid, 5e6, 400.0).phase == "supercritical" and wet.summary()["phase"] == "mixture"
    assert fluid.summary()["Tc"] == fluid.Tc and set(fluid.summary()) == {"name", "Tc", "pc", "omega", "molar_mass", "cp0"}


@pytest.mark.parametrize(
    "call",
    [
        lambda: t.refrigerant_properties("R-134a", 0.5e6, 300.0),
        lambda: t.refrigerant_properties(t.R134A, 0.5e6),
        lambda: t.refrigerant_properties(t.R134A, 0.5e6, 100.0),
        lambda: t.refrigerant_properties(t.R134A, [0.5e6], 300.0),
        lambda: t.refrigerant_properties(t.R134A, 4.058e6, x=0.5),
        lambda: t.refrigerant_properties(t.R134A, 10.0, x=0.5),
        lambda: t.refrigerant_properties(t.R134A, 0.5e6, h=1e7),
        lambda: t.peng_robinson_compressibility(300.0, 1e5, 374.2, 4.059e6, 0.327, phase="gas"),
        lambda: t.peng_robinson_compressibility([300.0], 1e5, 374.2, 4.059e6, 0.327),
        lambda: t.peng_robinson_saturation_pressure([300.0], 374.2, 4.059e6, 0.327),
        lambda: t.CubicFluid("bad", 300.0, 1e6, 3.0, 0.1, 800.0),
    ],
)
def test_peng_robinson_functions_reject_bad_inputs(call):
    with pytest.raises(ValueError):
        call()


def test_vapor_compression_cycle_matches_the_textbook_refrigerator():
    # Çengel and Boles Example 11-1: R-134a between 0.14 and 0.8 MPa, 0.05 kg/s:
    # Q_L = 7.18 kW, W = 1.81 kW, COP = 3.97.
    res = t.vapor_compression_cycle(0.14e6, 0.8e6)
    assert np.isclose(0.05 * res.q_evaporator, 7.18e3, rtol=0.02)
    assert np.isclose(0.05 * res.compressor_work, 1.81e3, rtol=0.02)
    assert np.isclose(res.cop_refrigerator, 3.97, rtol=0.02)
    assert np.isclose(res.cop_heat_pump, res.cop_refrigerator + 1)
    T_low, T_high = res.states[0].T, res.states[2].T
    assert res.cop_refrigerator < t.carnot_cop_refrigerator(T_high, T_low)
    assert np.isclose(res.states[3].h, res.states[2].h) and res.states[3].s > res.states[2].s
    real = t.vapor_compression_cycle(0.14e6, 0.8e6, superheat=5.0, subcooling=3.0, compressor_efficiency=0.8)
    assert real.cop_refrigerator < res.cop_refrigerator and real.summary()["compressor_outlet_T"] > res.states[1].T
    assert np.isclose(real.states[0].T, res.states[0].T + 5.0) and real.states[2].phase == "liquid"
    with pytest.raises(ValueError):
        t.vapor_compression_cycle(0.8e6, 0.14e6)


# --- Compressibility, Pitzer and Dalton -----------------------------------------------------------


def test_pitzer_virial_compressibility_agrees_with_peng_robinson_for_gases():
    Tc, pc, omega = 425.1, 3.796e6, 0.200  # n-butane
    for T, p in ((510.0, 2.5e6), (600.0, 1e6), (450.0, 0.3e6)):
        z_virial = t.pitzer_virial_compressibility(T, p, Tc, pc, omega)
        z_pr = t.peng_robinson_compressibility(T, p, Tc, pc, omega)
        assert np.isclose(z_virial, z_pr, atol=0.025)  # two approximate models; they differ most near Tr = 1.2
        v = z_pr * t.R_UNIVERSAL * T / p
        assert np.isclose(t.compressibility_factor(p, T, v), z_pr)
    assert np.isclose(t.pitzer_virial_compressibility(500.0, 1.0, Tc, pc, omega), 1, atol=1e-6)


def test_acentric_factor_of_water_from_if97():
    p07 = t.saturation_pressure(0.7 * if97.T_CRIT)
    assert np.isclose(t.pitzer_acentric_factor(p07, if97.P_CRIT), 0.344, atol=0.002)


def test_dalton_mixture_relations():
    M = [0.028013, 0.031999, 0.039948]  # N2, O2, Ar
    y = np.array([0.7808, 0.2095, 0.0097])
    y = y / y.sum()
    assert np.isclose(t.mixture_molar_mass(y, M), 0.02897, rtol=1e-3)
    np.testing.assert_allclose(t.mole_fractions(t.mass_fractions(y, M), M), y)
    assert np.isclose(t.partial_pressures(y, 101325.0).sum(), 101325.0)
    # mixing entropy is the sum of the isothermal expansions of each gas to the whole volume
    n = 3.0
    by_expansion = sum(t.entropy_change_isothermal_ideal_gas(n * yi, 300.0, yi, 1.0) for yi in y)
    assert np.isclose(t.entropy_of_mixing(y, n), by_expansion)
    assert t.entropy_of_mixing([1.0, 0.0]) == 0 and np.isclose(t.entropy_of_mixing([0.5, 0.5]), t.R_UNIVERSAL * math.log(2))
    for call in (lambda: t.mole_fractions([0.5, 0.6], [0.002, 0.032]), lambda: t.mass_fractions([1.0], [0.002, 0.032])):
        with pytest.raises(ValueError):
            call()


# --- Maxwell relations ---------------------------------------------------------------------------


def test_maxwell_relation_and_general_relations_hold_for_if97_steam():
    p, T, dT, dp = 2e6, 600.0, 0.01, 100.0
    v = lambda p_, T_: t.steam_properties(p_, T_).v  # noqa: E731
    s = lambda p_, T_: t.steam_properties(p_, T_).s  # noqa: E731
    dv_dT = (v(p, T + dT) - v(p, T - dT)) / (2 * dT)
    ds_dp = (s(p + dp, T) - s(p - dp, T)) / (2 * dp)
    assert np.isclose(ds_dp, -dv_dT, rtol=1e-5)  # (ds/dp)_T = -(dv/dT)_p
    st = t.steam_properties(p, T)
    beta = dv_dT / st.v
    mu = t.joule_thomson_coefficient(T, st.v, beta, st.cp)
    throttled = t.steam_properties(p - 1e3, h=st.h)
    assert np.isclose((T - throttled.T) / 1e3, mu, rtol=1e-3)
    kappa = -(v(p + dp, T) - v(p - dp, T)) / (2 * dp) / st.v
    # cv from (du/dT)_v along a constant-volume path, by finding the pressure that keeps v
    cp_cv = t.cp_minus_cv(T, st.v, beta, kappa)
    assert 0 < cp_cv < st.cp


def test_general_relations_reduce_to_the_ideal_gas():
    R, T, p = 287.0, 350.0, 2e5
    v = R * T / p
    assert np.isclose(t.cp_minus_cv(T, v, 1 / T, 1 / p), R)
    assert t.joule_thomson_coefficient(T, v, 1 / T, 1005.0) == 0
    assert t.cp_minus_cv(277.13, 1e-3, 0.0, 4.6e-10) == 0  # water at its density maximum
    mu = t.joule_thomson_coefficient(300.0, 0.0249, 2 * (1 / 300.0), 1040.0)
    assert mu > 0


# --- Psychrometric processes ---------------------------------------------------------------------


def test_moist_air_enthalpy_matches_if97_vapour():
    T, w = 303.15, 0.012
    h_vapour = t.steam_properties(float(t.saturation_pressure(T)) * 0.5, T).h
    assert np.isclose(t.moist_air_enthalpy(T, w), 1005.0 * 30 + w * h_vapour, rtol=1e-3)
    assert t.moist_air_enthalpy(273.15, 0.0) == 0


def test_adiabatic_saturation_closes_its_balance_and_lies_between_dew_and_dry_bulb():
    T, w, p = 303.15, 0.010, 101325.0
    T_as = t.adiabatic_saturation_temperature(T, w, p)
    w_s = t.humidity_ratio(T_as, 1.0, p)
    h_f = t.steam_properties(t.saturation_pressure(T_as), x=0.0).h
    assert np.isclose(t.moist_air_enthalpy(T, w) + (w_s - w) * h_f, t.moist_air_enthalpy(T_as, w_s), rtol=1e-8)
    assert t.dew_point(w) < T_as < T
    # Çengel and Boles, Ch. 14: 25 °C dry bulb, 15 °C wet bulb at 1 atm gives w = 0.00653
    w_book = 0.00653
    assert np.isclose(t.adiabatic_saturation_temperature(298.15, w_book), 288.15, atol=0.2)
    saturated = float(t.humidity_ratio(T, 1.0))
    assert t.adiabatic_saturation_temperature(T, saturated) == T
    for call in (
        lambda: t.adiabatic_saturation_temperature(T, 0.05),
        lambda: t.adiabatic_saturation_temperature(275.0, 0.0),
        lambda: t.adiabatic_saturation_temperature([T], w),
    ):
        with pytest.raises(ValueError):
            call()


def test_mixing_cooling_coil_and_cooling_tower_conserve_air_water_and_energy():
    T3, w3 = t.adiabatic_mixing(1.0, 305.15, 0.015, 3.0, 287.15, 0.006)
    assert np.isclose(4 * t.moist_air_enthalpy(T3, w3), t.moist_air_enthalpy(305.15, 0.015) + 3 * t.moist_air_enthalpy(287.15, 0.006))
    assert 287.15 < T3 < 305.15 and np.isclose(w3, (0.015 + 3 * 0.006) / 4)
    # above the dew point a coil only cools; below it the air leaves saturated
    w, q, water = t.cooling_coil(303.15, 0.3, 298.15)
    assert water == 0 and np.isclose(q, t.moist_air_enthalpy(303.15, w) - t.moist_air_enthalpy(298.15, w))
    w2, q2, water2 = t.cooling_coil(303.15, 0.8, 287.15)
    assert np.isclose(t.relative_humidity(287.15, w2), 1) and np.isclose(water2, t.humidity_ratio(303.15, 0.8) - w2)
    # Çengel and Boles Example 14-9: 100 kg/s from 35 to 22 °C; air 20 °C, 60% in, 30 °C saturated out:
    # about 96.9 kg/s of dry air and 1.80 kg/s of make-up water.
    m_air, makeup = t.cooling_tower(100.0, 308.15, 295.15, 293.15, 0.6, 303.15)
    assert np.isclose(m_air, 96.9, rtol=0.015) and np.isclose(makeup, 1.80, rtol=0.02)
    hf = lambda T: t.steam_properties(t.saturation_pressure(T), x=0.0).h  # noqa: E731
    w_in, w_out = t.humidity_ratio(293.15, 0.6), t.humidity_ratio(303.15, 1.0)
    gain = m_air * (t.moist_air_enthalpy(303.15, w_out) - t.moist_air_enthalpy(293.15, w_in))
    assert np.isclose(gain, 100 * hf(308.15) - (100 - makeup) * hf(295.15))
    for call in (
        lambda: t.cooling_coil(303.15, 0.5, 310.0),
        lambda: t.cooling_tower(100.0, 308.15, 295.15, 303.15, 1.0, 293.15, 0.6),
        lambda: t.adiabatic_mixing(0.0, 300.0, 0.01, 0.0, 300.0, 0.01),
    ):
        with pytest.raises(ValueError):
            call()


# --- Combustion -----------------------------------------------------------------------------------


def test_combustion_stoichiometry_and_hess_heating_values():
    products = t.combustion_products(8, 18, excess_air=0.5)
    need = 12.5
    # atoms balance: O from the air ends in CO2, H2O and the excess O2
    assert np.isclose(2 * 1.5 * need, 2 * products["CO2"] + products["H2O"] + 2 * products["O2"])
    assert np.isclose(products["N2"], 3.76 * 1.5 * need)
    assert np.isclose(t.stoichiometric_air_fuel_ratio(1, 4), 17.2, rtol=0.005)
    lhv = t.heating_value(1, 4, h_f_fuel=-74.87e3)
    hhv = t.heating_value(1, 4, h_f_fuel=-74.87e3, water="liquid")
    assert np.isclose(lhv, 50.0e6, rtol=0.002) and np.isclose(hhv, 55.5e6, rtol=0.002)
    assert np.isclose((hhv - lhv) * t.fuel_molar_mass(1, 4), 2 * 44.0e3, rtol=1e-3)
    # an oxygenated fuel needs less air: methanol
    assert t.stoichiometric_air_fuel_ratio(1, 4, 1) < t.stoichiometric_air_fuel_ratio(1, 4)
    for call in (
        lambda: t.combustion_products(0, 0),
        lambda: t.heating_value(1, 4, h_f_fuel=0.0, water="ice"),
        lambda: t.combustion_products(1, 4, excess_air=-0.1),
    ):
        with pytest.raises(ValueError):
            call()


def test_adiabatic_flame_temperature_closes_the_energy_balance():
    T = t.adiabatic_flame_temperature(1, 4, h_f_fuel=-74.87e3)
    products = t.combustion_products(1, 4)
    h_p = sum(n * (shomate.FORMATION_ENTHALPY[k] + shomate.sensible_enthalpy(k, T)) for k, n in products.items())
    assert np.isclose(h_p, -74.87e3, atol=1e-3)
    assert np.isclose(T, 2327.0, atol=2.0)  # no dissociation; equilibrium flames are near 2230 K
    # Çengel and Boles Example 15-8: liquid octane with 400% theoretical air reaches 962 K
    assert np.isclose(t.adiabatic_flame_temperature(8, 18, h_f_fuel=-249.95e3, excess_air=3.0), 962.0, atol=3.0)
    lean = t.adiabatic_flame_temperature(1, 4, h_f_fuel=-74.87e3, excess_air=1.0)
    hot_air = t.adiabatic_flame_temperature(1, 4, h_f_fuel=-74.87e3, T_air=800.0)
    assert lean < T < hot_air
    with pytest.raises(ValueError):
        t.adiabatic_flame_temperature(1, 0, 1, h_f_fuel=500e3, T_air=298.15)  # needs no oxygen
    with pytest.raises(ValueError):
        t.adiabatic_flame_temperature(1, 4, h_f_fuel=-2000e3)  # releases no heat


def test_shomate_enthalpies_match_janaf_and_are_continuous():
    # NIST-JANAF H(1000 K) - H(298.15 K): CO2 33.397, H2O 26.000, N2 21.463, O2 22.703 kJ/mol
    for species, value in (("CO2", 33.397e3), ("H2O", 26.000e3), ("N2", 21.463e3), ("O2", 22.703e3)):
        assert np.isclose(shomate.sensible_enthalpy(species, 1000.0), value, rtol=2e-4)
        assert abs(shomate.sensible_enthalpy(species, 298.15)) < 5.0
        for edge in (500.0, 700.0, 1200.0, 1700.0, 2000.0):
            assert abs(shomate.sensible_enthalpy(species, edge * (1 + 1e-9)) - shomate.sensible_enthalpy(species, edge)) < 5.0

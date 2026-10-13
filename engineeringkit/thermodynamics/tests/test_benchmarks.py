import numpy as np
import pytest
from scipy.integrate import quad

from engineeringkit import thermodynamics as t


def test_isothermal_work_equals_pdv_and_tds():
    n, T, V1, V2 = 2.0, 350.0, 0.01, 0.04
    integral = quad(lambda V: t.ideal_gas_pressure(n, T, V), V1, V2)[0]
    work = t.reversible_isothermal_work_ideal_gas(n, T, V1, V2)
    assert np.isclose(integral, work)
    assert np.isclose(work, T * t.entropy_change_isothermal_ideal_gas(n, T, V1, V2))
    assert np.isclose(t.ideal_gas_temperature(t.ideal_gas_pressure(n, T, V1), V1, n), T)


def test_carnot_energy_balance_and_cop_identity():
    eta = t.carnot_efficiency(600, 300)
    qin = 1000
    qout = (1 - eta) * qin
    assert np.isclose(t.heat_engine_work(qin, qout), eta * qin)
    assert np.isclose(t.carnot_cop_heat_pump(600, 300), t.carnot_cop_refrigerator(600, 300) + 1)
    assert t.clausius_inequality([qin / 600, -qout / 300])


def test_otto_efficiency_from_four_state_first_law():
    r, gamma = 8.0, 1.4
    T1, T3 = 300.0, 1800.0
    T2 = T1 * r ** (gamma - 1)
    T4 = T3 / r ** (gamma - 1)
    efficiency = 1 - (T4 - T1) / (T3 - T2)
    assert np.isclose(t.otto_efficiency(r, gamma), efficiency)


def test_diesel_otto_limit_and_compression_bounds():
    assert t.diesel_efficiency(8, 1) == t.otto_efficiency(8)
    assert np.isclose(t.diesel_efficiency(8, 1 + 1e-9), t.otto_efficiency(8), atol=1e-9)
    assert 0 < t.diesel_efficiency(8, 2) < t.otto_efficiency(8)
    assert t.brayton_efficiency(1) == 0


def test_rankine_first_law():
    # J/kg: pump outlet, boiler outlet, turbine outlet, condenser outlet.
    h1, h2, h3, h4 = 210e3, 3200e3, 2100e3, 200e3
    eta = t.rankine_ideal_efficiency(h1, h2, h3, h4)
    assert np.isclose(eta, 1 - (h3 - h4) / (h2 - h1))
    assert t.exergy_destruction(300, 2) == 600


@pytest.mark.parametrize("args", [(300, 300), (300, 400), (300, 0), (np.nan, 200)])
def test_invalid_carnot_reservoirs(args):
    for fn in (t.carnot_efficiency, t.carnot_cop_heat_pump, t.carnot_cop_refrigerator):
        with pytest.raises(ValueError):
            fn(*args)

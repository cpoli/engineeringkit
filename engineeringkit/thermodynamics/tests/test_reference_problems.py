"""Textbook reference problems (ROADMAP item B4): worked examples checked against the published answers.

Sources: Y. A. Çengel and M. A. Boles, *Thermodynamics: An Engineering Approach*, 9th ed.
(McGraw-Hill, 2019), and M. J. Moran, H. N. Shapiro, D. D. Boettner and M. B. Bailey,
*Fundamentals of Engineering Thermodynamics*, 9th ed. (Wiley, 2018). Each answer is met to the
digits the book prints (half a unit in the last place) unless a comment says otherwise; the
differences of method are listed in the guide, "Textbook reference problems".
"""

import pytest

from engineeringkit import thermodynamics as t


def printed(value, published, decimals):
    """True if ``value`` rounds to ``published`` at ``decimals`` places (with float slack)."""
    return abs(value - published) <= 0.5 * 10.0**-decimals + 1e-12


@pytest.mark.parametrize(
    "args, kwargs, efficiency, quality",
    [
        # Çengel and Boles, 9th ed., Example 10-1: 3 MPa, 350 °C, 75 kPa; eta = 0.260, x4 = 0.886.
        ((3e6, 75e3, 623.15), {}, 0.260, 0.886),
        # Example 10-3a: the condenser lowered to 10 kPa; eta = 0.334, x4 = 0.813.
        ((3e6, 10e3, 623.15), {}, 0.334, 0.813),
        # Example 10-3b: superheated to 600 °C; eta = 0.373, x4 = 0.915.
        ((3e6, 10e3, 873.15), {}, 0.373, 0.915),
        # Example 10-3c: the boiler raised to 15 MPa; eta = 0.430, x4 = 0.804.
        ((15e6, 10e3, 873.15), {}, 0.430, 0.804),
        # Example 10-4: reheated at 4 MPa to 600 °C; eta = 0.450, x6 = 0.896.
        ((15e6, 10e3, 873.15), {"reheat_pressure": 4e6, "T_reheat": 873.15}, 0.450, 0.896),
    ],
)
def test_cengel_boles_rankine_plants(args, kwargs, efficiency, quality):
    res = t.rankine_cycle(*args, **kwargs)
    assert printed(res.efficiency, efficiency, 3)
    assert printed(res.turbine_exit_quality, quality, 3)


def test_cengel_boles_example_10_5_open_feedwater_heater():
    # Çengel and Boles, 9th ed., Example 10-5: 15 MPa, 600 °C, 10 kPa, an open heater at 1.2 MPa;
    # y = 0.2270, eta = 0.463. The book's bleed enthalpy (2860.2 kJ/kg) is interpolated in the
    # superheat table, so y is met to three places, not four.
    res = t.rankine_cycle(15e6, 10e3, 873.15, feedwater_pressure=1.2e6)
    assert printed(res.efficiency, 0.463, 3)
    assert printed(res.bleed_fraction, 0.227, 3)


def test_moran_example_8_1_ideal_rankine_cycle():
    # Moran et al., 9th ed., Example 8.1: saturated vapour at 8.0 MPa, 0.008 MPa, 100 MW net;
    # eta = 37.1%, m = 3.77e5 kg/h, x2 = 0.6745. Moran's steam tables give s_g(8 MPa) = 5.7432
    # against IF97's 5.7450 kJ/(kg K), so x2 is met to three places, not four.
    res = t.rankine_cycle(8e6, 8e3)
    assert printed(res.efficiency, 0.371, 3)
    assert printed(res.turbine_exit_quality, 0.675, 3)
    assert printed(100e6 / res.net_work * 3600 / 1e5, 3.77, 2)


def test_moran_example_8_2_rankine_cycle_with_irreversibilities():
    # Moran et al., 9th ed., Example 8.2: Example 8.1 with 85% turbine and pump efficiencies;
    # eta = 31.4%.
    res = t.rankine_cycle(8e6, 8e3, turbine_efficiency=0.85, pump_efficiency=0.85)
    assert printed(res.efficiency, 0.314, 3)


def test_cengel_boles_example_9_5_cold_air_standard_brayton_cycle():
    # Çengel and Boles, 9th ed., Example 9-5: pressure ratio 8, 300 K to 1300 K. Under the
    # cold-air-standard assumptions (k = 1.4) the book gives eta = 0.448.
    res = t.brayton_cycle(8.0, 300.0, 1300.0)
    assert printed(res.efficiency, 0.448, 3)
    assert printed(t.brayton_efficiency(8.0), 0.448, 3)


def test_cengel_boles_example_15_8_adiabatic_flame_temperature():
    # Çengel and Boles, 9th ed., Example 15-8: liquid octane (h_f = -249,950 kJ/kmol) burned
    # steadily with 400% theoretical air, both at 25 °C; T_prod = 962 K.
    T = t.adiabatic_flame_temperature(8, 18, h_f_fuel=-249.95e3, excess_air=3.0)
    assert printed(T, 962.0, 0)


# --- Differences of method (see the guide, "Textbook reference problems") -------------------------


def test_cengel_boles_example_11_1_ideal_refrigerator_with_peng_robinson():
    # Çengel and Boles, 9th ed., Example 11-1: R-134a between 0.14 and 0.8 MPa, 0.05 kg/s;
    # Q_L = 7.18 kW, W = 1.81 kW, Q_H = 8.99 kW, COP_R = 3.97. The kit's Peng-Robinson R-134a
    # meets the work to the printed digits; its latent heat at 0.14 MPa is low, so Q_L and the
    # COP are 0.8% under the tables.
    res = t.vapor_compression_cycle(0.14e6, 0.8e6)
    assert printed(0.05 * res.compressor_work / 1e3, 1.81, 2)
    assert abs(0.05 * res.q_evaporator / 7.18e3 - 1) < 0.01
    assert abs(0.05 * res.q_condenser / 8.99e3 - 1) < 0.01
    assert abs(res.cop_refrigerator / 3.97 - 1) < 0.01


def test_cengel_boles_example_9_6_actual_gas_turbine_with_constant_specific_heats():
    # Çengel and Boles, 9th ed., Example 9-6: Example 9-5 with an 80% compressor and an 85%
    # turbine; from the air tables r_bw = 0.592, eta = 0.266, T4a = 853 K. With constant
    # cp = 1.005 kJ/(kg K) and k = 1.4, the room-temperature values, the kit gives 0.615, 0.274
    # and 805 K: the hot gas's larger cp and smaller k make the real expansion drop less in T.
    res = t.brayton_cycle(8.0, 300.0, 1300.0, compressor_efficiency=0.8, turbine_efficiency=0.85)
    assert printed(res.back_work_ratio, 0.615, 3) and 0 < res.back_work_ratio - 0.592 < 0.025
    assert printed(res.efficiency, 0.274, 3) and 0 < res.efficiency - 0.266 < 0.01
    T4a = res.T[res.labels.index("turbine 1 outlet")]
    assert printed(T4a, 805.0, 0) and 40 < 853.0 - T4a < 55

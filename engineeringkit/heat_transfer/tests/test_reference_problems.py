"""Worked examples from Incropera et al., checked to the precision the book prints.

F. P. Incropera, D. P. DeWitt, T. L. Bergman and A. S. Lavine, *Fundamentals
of Heat and Mass Transfer*, 6th ed. (Hoboken, NJ: Wiley, 2007). The book
works in rounded property values and gives answers to three significant
figures; each test uses the book's inputs and accepts half a unit in the
last printed digit.
"""

import numpy as np
import pytest

from engineeringkit import heat_transfer as h


def test_example_1_1_furnace_wall_conduction():
    # 0.15 m fireclay wall, k = 1.7 W/(m K), 0.5 m by 1.2 m, faces at 1400 K and 1150 K.
    q = h.fourier_heat_rate(1.7, 0.5 * 1.2, 1400 - 1150, 0.15)
    assert q / (0.5 * 1.2) == pytest.approx(2833, abs=0.5)
    assert q == pytest.approx(1700, abs=0.5)


def test_example_1_2_uninsulated_steam_pipe():
    # D = 70 mm, Ts = 200 C, air and surroundings at 25 C, eps = 0.8, h = 15 W/(m^2 K).
    Ts, T_inf, eps, A = 473.0, 298.0, 0.8, np.pi * 0.070
    assert eps * h.stefan_boltzmann_emissive_power(Ts) == pytest.approx(2270, abs=5)
    q_conv = h.convection_heat_rate(15, A, Ts, T_inf)
    q_rad = h.radiation_heat_rate(eps, A, Ts, T_inf)
    assert q_conv == pytest.approx(577, abs=0.5)
    assert q_rad == pytest.approx(421, abs=0.5)
    assert q_conv + q_rad == pytest.approx(998, abs=0.5)


def test_example_5_1_thermocouple_junction_time_constant():
    # Spherical junction, h = 400 W/(m^2 K), k = 20 W/(m K), c = 400 J/(kg K),
    # rho = 8500 kg/m^3, sized for a 1 s time constant; gas at 200 C, junction from 25 C.
    rho, c, hc, k = 8500, 400, 400, 20
    D = 6 * hc * 1.0 / (rho * c)
    assert D == pytest.approx(7.06e-4, abs=5e-7)
    assert h.biot_number(hc, D / 6, k) == pytest.approx(2.35e-3, abs=5e-6)
    A, V = np.pi * D**2, np.pi * D**3 / 6
    # One time constant leaves 1/e of the initial difference.
    assert h.lumped_temperature(298, 473, hc, A, rho, c, V, 1.0) == pytest.approx(473 - 175 / np.e, rel=1e-12)
    assert h.lumped_temperature(298, 473, hc, A, rho, c, V, 5.2) - 273 == pytest.approx(199, abs=0.5)


def test_example_7_4_cylinder_in_crossflow_churchill_bernstein():
    # D = 12.7 mm cylinder in 10 m/s air, film-temperature properties: Re_D = 6071, Pr = 0.700.
    assert h.churchill_bernstein_nusselt(6071, 0.700) == pytest.approx(40.6, abs=0.05)


def test_example_11_1_counterflow_oil_cooler():
    # Oil 0.1 kg/s, cp = 2131 J/(kg K), 100 -> 60 C; water 0.2 kg/s, cp = 4178 J/(kg K),
    # enters at 30 C in a 25 mm inner tube; h_i = 2250 and h_o = 38.8 W/(m^2 K).
    q = 0.1 * 2131 * (100 - 60)
    assert q == pytest.approx(8524, abs=0.5)
    Tc_out = 30 + q / (0.2 * 4178)
    assert Tc_out == pytest.approx(40.2, abs=0.05)
    dT_lm = h.log_mean_temperature_difference(100 - Tc_out, 60 - 30)
    assert dT_lm == pytest.approx(43.2, abs=0.05)
    U = h.overall_u(1 / 2250, 1 / 38.8)
    L = q / h.heat_exchanger_rate(U, np.pi * 0.025, dT_lm)
    assert L == pytest.approx(65.9, abs=0.05)

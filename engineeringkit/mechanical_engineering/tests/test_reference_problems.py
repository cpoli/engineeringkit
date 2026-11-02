"""Worked problems from Shigley and Rao, checked against the answers the books print (ROADMAP B4).

R. G. Budynas and J. K. Nisbett, *Shigley's Mechanical Engineering Design*, 10th ed.
(McGraw-Hill, 2015), Example 14-1 and the worked solutions to end-of-chapter problems in the
publisher's solutions manual; S. S. Rao, *Mechanical Vibrations*, 5th ed. (Prentice Hall,
2011). Each tolerance is half a unit in the last printed digit; the differences of method are in
the guide. US customary data are converted to SI here and the answers converted back.
"""

import numpy as np
import pytest

from engineeringkit import mechanical_engineering as m

IN = 0.0254
LBF = 4.4482216152605
PSI = LBF / IN**2
HP = 550 * LBF * 0.3048


def test_shigley_problem_10_3_music_wire_spring_rate_and_load_to_close_solid():
    # Shigley, 10th ed., Problem 10-3: music wire, d = 2.5 mm, OD = 31 mm, plain ground ends,
    # Nt = 14 (Na = 13), G = 81.0 GPa: k = 1.314 N/mm. Sut = 2211/d^0.145 = 1936 MPa,
    # Ssy = 0.45 Sut = 871.2 MPa, K_B = 1.117 and the load at yield is F_s = 167.9 N.
    d, D = 2.5e-3, 31e-3 - 2.5e-3
    assert m.spring_index(D, d) == pytest.approx(11.4)
    assert m.spring_rate(d, D, 13, 81.0e9) / 1e3 == pytest.approx(1.314, abs=5e-4)
    assert m.bergstrasser_factor(11.4) == pytest.approx(1.117, abs=5e-4)
    Ssy = 0.45 * 2211e6 / 2.5**0.145
    assert Ssy / 1e6 == pytest.approx(871.2, abs=0.05)
    assert Ssy / m.spring_shear_stress(1.0, D, d, correction="bergstrasser") == pytest.approx(167.9, abs=0.05)


def test_shigley_problem_10_5_oil_tempered_spring_closed_solid():
    # Shigley, 10th ed., Problem 10-5: oil-tempered wire, d = 0.2 in, D = 2 in, Nt = 12,
    # L0 = 5 in, squared ends (Ls = 2.6 in, Na = 10), G = 11.2 Mpsi: k = 28 lbf/in,
    # F_s = 67.2 lbf, K_B = 1.135, tau_s = 48.56 kpsi, Sut = 198.6 kpsi and n_s = 2.04 on
    # Ssy = 0.50 Sut.
    k = m.spring_rate(0.2 * IN, 2 * IN, 10, 11.2e6 * PSI) / (LBF / IN)
    assert k == pytest.approx(28.0, abs=0.5)
    Fs = k * (5 - 2.6)
    assert Fs == pytest.approx(67.2, abs=0.05)
    assert m.bergstrasser_factor(10.0) == pytest.approx(1.135, abs=5e-4)
    tau = m.spring_shear_stress(Fs * LBF, 2 * IN, 0.2 * IN, correction="bergstrasser") / PSI / 1e3
    assert tau == pytest.approx(48.56, abs=5e-3)
    Sut = 147 / 0.2**0.187
    assert Sut == pytest.approx(198.6, abs=0.05)
    # Unrounded n_s = 2.0450 sits on the rounding boundary; Shigley divides the printed 99.3 by
    # the printed 48.56.
    assert 0.50 * Sut / tau == pytest.approx(2.045, abs=5e-4)
    assert 99.3 / 48.56 == pytest.approx(2.04, abs=5e-3)


def test_shigley_problem_7_1_de_goodman_shaft_diameter():
    # Shigley, 10th ed., Problem 7-1: Kf = 2.2, Kfs = 1.8, Ma = 70 N m, Ta = 45 N m,
    # Mm = 55 N m, Tm = 35 N m, Se = 210 MPa, Sut = 700 MPa, n = 2: DE-Goodman d = 27.27 mm.
    d = m.shaft_diameter_goodman(210e6, 700e6, Ma=70.0, Tm=35.0, Mm=55.0, Ta=45.0, safety_factor=2.0, Kf=2.2, Kfs=1.8)
    assert d * 1e3 == pytest.approx(27.27, abs=5e-3)


@pytest.mark.parametrize(
    ("load", "life", "rating"),
    [
        # Shigley, 10th ed., Problems 11-8, 11-9 and 11-10, ball bearings at R = 0.90 by
        # Eq. (11-3): C10 = 20 kN, 5050 lbf and 24.9 kN.
        (2e3, 1e9, (20e3, 0.5e3)),
        (800 * LBF, 12_000 * 350 * 60, (5050 * LBF, 5 * LBF)),
        (4e3, 8_000 * 500 * 60, (24.9e3, 0.05e3)),
    ],
    ids=["11-8", "11-9", "11-10"],
)
def test_shigley_problems_11_8_to_11_10_catalog_rating_for_a_design_life(load, life, rating):
    # bearing_l10_life gives the life for a rating; the printed rating is right to its last
    # digit when the design life lies between the lives of the ratings half a digit either side.
    C, half_digit = rating
    assert m.bearing_l10_life(C - half_digit, load) <= life <= m.bearing_l10_life(C + half_digit, load)


def test_shigley_example_14_1_lewis_rating_of_a_stock_spur_gear():
    # Shigley, 10th ed., Example 14-1: P = 8 teeth/in, F = 1.5 in, N = 16, 20 deg full depth,
    # AISI 1020 as rolled (Sy = 30 kpsi), n_d = 3, 1200 rev/min: Y = 0.296, V = 628.3 ft/min,
    # cut-profile K_v = 1.52, W_t = 365 lbf and H = 6.95 hp. Shigley rounds K_v before dividing;
    # the unrounded K_v gives 364 lbf (see the guide).
    Y = m.lewis_form_factor(16)
    assert Y == pytest.approx(0.296)
    V = m.belt_speed(2 * IN, 1200.0)
    assert V / 0.00508 == pytest.approx(628.3, abs=0.05)
    Kv = m.barth_velocity_factor(V)
    assert Kv == pytest.approx(1.52, abs=5e-3)
    allowable = 30e3 / 3 * PSI
    Wt = allowable / m.lewis_bending_stress(1.0, 1.5 * IN, IN / 8, Y, velocity_factor=1.52)
    assert Wt / LBF == pytest.approx(365, abs=0.5)
    assert Wt * V / HP == pytest.approx(6.95, abs=5e-3)
    assert allowable / m.lewis_bending_stress(1.0, 1.5 * IN, IN / 8, Y, velocity_factor=Kv) / LBF == pytest.approx(364, abs=0.5)


def test_rao_example_6_11_three_equal_masses_on_a_spring_chain():
    # Rao, 5th ed., Example 6.11: m1 = m2 = m3 = m and k1 = k2 = k3 = k, fixed at one end.
    # Rao prints omega = 0.44504, 1.2471 and 1.8025 sqrt(k/m) (lambda = 0.19806, 1.5553,
    # 3.2490). The first matches to the printed digits; Rao's roots of the cubic
    # lambda^3 - 5 lambda^2 + 6 lambda - 1 = 0 are low-precision for the other two, whose exact
    # values are 1.24698 and 1.80194 (see the guide).
    K = [[2.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 1.0]]
    w = m.normal_modes(np.eye(3), K).frequencies
    assert w[0] == pytest.approx(0.44504, abs=5e-6)
    assert w[0] ** 2 == pytest.approx(0.19806, abs=5e-6)
    np.testing.assert_allclose(w[1:], [1.2471, 1.8025], atol=6e-4)
    lam = w**2
    np.testing.assert_allclose(lam**3 - 5 * lam**2 + 6 * lam - 1, 0.0, atol=1e-12)
    assert abs(1.8025**2 - lam[2]) > 1e-3

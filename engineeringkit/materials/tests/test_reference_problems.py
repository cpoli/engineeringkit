"""Textbook reference problems for materials (ROADMAP item B4).

Each test is one worked example from Callister and Rethwisch, *Materials Science and
Engineering: An Introduction*, 9th edition (Wiley, 2014), cited by its example number. The
given data are the book's, and each answer is checked to the precision the book prints:
the tolerance is half a unit in its last printed digit.
"""

import numpy as np

from engineeringkit import materials as m

AVOGADRO = 6.022e23


def test_callister_example_3_3_copper_density():
    # FCC copper, R = 0.128 nm, A = 63.5 g/mol: rho = 8.89 g/cm^3.
    a = m.cubic_lattice_parameter(0.128e-9, "fcc")
    rho = m.theoretical_density(4, 63.5e-3, a**3)
    assert abs(rho - 8.89e3) < 5.0


def test_callister_example_5_1_steady_carbon_flux():
    # Iron plate at 700 C, 1.2 and 0.8 kg/m^3 at 5 and 10 mm, D = 3e-11 m^2/s: J = 2.4e-9 kg/(m^2 s).
    J = m.fick_flux(3e-11, (0.8 - 1.2) / (10e-3 - 5e-3))
    assert abs(J - 2.4e-9) < 0.05e-9


def test_callister_example_5_2_carburizing_time():
    # 0.25 wt% C steel, surface 1.20 wt% C, 0.80 wt% C at 0.5 mm, D = 1.6e-11 m^2/s at 950 C:
    # t = 25,400 s = 7.1 h. The book interpolates its erf table to z = 0.392; the exact inverse
    # error function gives z = 0.3926 and 25,370 s, which still rounds to the printed 25,400 s.
    # The hours are the book's rounding of 25,400 s, not a separate answer (25,370 s is 7.05 h).
    t = m.diffusion_time(0.5e-3, 0.80, 1.6e-11, 0.25, 1.20)
    assert abs(t - 25400) < 50


def test_callister_example_5_4_copper_in_aluminum_diffusivity():
    # Cu in Al at 500 C, D0 = 6.5e-5 m^2/s and Qd = 136 kJ/mol (Table 5.2): D = 4.2e-14 m^2/s.
    D = m.diffusion_coefficient(6.5e-5, 136e3, 500 + 273.0)
    assert abs(D - 4.2e-14) < 0.05e-14


def test_callister_example_6_1_copper_elastic_elongation():
    # Copper, 276 MPa on a 305 mm specimen, E = 110 GPa: delta_l = 0.77 mm.
    delta_l = 305e-3 * 276e6 / m.hooke_stress(110e9, 1.0)
    assert abs(delta_l - 0.77e-3) < 0.005e-3


def test_callister_example_6_2_brass_rod_load_for_diameter_change():
    # Brass rod, d0 = 10 mm, E = 97 GPa, nu = 0.34, diameter reduced by 2.5e-3 mm:
    # eps_x = -2.5e-4, eps_z = 7.35e-4, F = 5600 N.
    eps_x = m.normal_strain(-2.5e-6, 10e-3)
    eps_z = -eps_x / 0.34
    assert abs(eps_x + 2.5e-4) < 0.05e-4
    assert abs(eps_z - 7.35e-4) < 0.005e-4
    sigma = m.hooke_stress(97e9, eps_z)
    F = sigma / m.normal_stress(1.0, np.pi * (10e-3) ** 2 / 4)
    assert abs(F - 5600) < 50


def test_callister_example_8_1_glass_maximum_flaw_size():
    # Glass at 40 MPa, E = 69 GPa, gamma_s = 0.3 J/m^2: largest surface flaw a = 8.2 um.
    # Griffith's stress scales as a^-1/2, so the flaw at 40 MPa is (sigma_f(1 m) / 40 MPa)^2 m.
    a = (m.griffith_fracture_stress(69e9, 0.3, 1.0) / 40e6) ** 2
    assert abs(a - 8.2e-6) < 0.05e-6


def test_callister_design_example_8_2_rupture_lifetime():
    # S-590 at 140 MPa and 800 C, where Figure 8.32 gives T(20 + log t_r) = 24.0e3: t_r = 233 h.
    t_r = m.larson_miller_rupture_time(24.0e3, 800 + 273.0) / 3600
    assert abs(t_r - 233) < 0.5


def test_callister_example_9_2_lead_tin_phase_fractions():
    # 40 wt% Sn-60 wt% Pb at 150 C, C_alpha = 11 and C_beta = 98 wt% Sn (Example 9.1):
    # W_alpha = 0.67, W_beta = 0.33.
    W_alpha, W_beta = m.lever_rule(40.0, 11.0, 98.0)
    assert abs(W_alpha - 0.67) < 0.005
    assert abs(W_beta - 0.33) < 0.005


def test_callister_example_19_1_brass_thermal_stress_on_heating():
    # Brass rod held between rigid supports, stress-free at 20 C, alpha = 20.0e-6 /C and
    # E = 100 GPa, compressive stress at most 172 MPa: delta_T = 86 C, T_f = 106 C.
    dT = -172e6 / m.constrained_thermal_stress(100e9, 20.0e-6, 1.0)
    assert abs(dT - 86) < 0.5
    assert abs(20 + dT - 106) < 0.5

"""Textbook reference problems for structural_engineering (ROADMAP item B4).

Each test is one worked example from R. C. Hibbeler, *Mechanics of Materials*, 10th edition
(Pearson, 2017), or J. M. Gere and B. J. Goodno, *Mechanics of Materials*, 7th edition
(Cengage, 2009), cited by its example number. The given data are the book's, converted to SI
where the book uses US customary units, and each answer is checked to the precision the book
prints: the tolerance is half a unit in its last printed digit.
"""

import numpy as np

from engineeringkit import structural_engineering as s

PSI = 6894.757293168  # Pa
INCH = 0.0254  # m
LBF = 4.4482216152605  # N


def test_hibbeler_example_8_1_allowable_vessel_pressure():
    # Cylinder of 4 ft inner diameter and 1/2 in wall, neither stress above 20 ksi: p = 417 psi;
    # a sphere of the same diameter and wall: p = 833 psi. Membrane stresses are linear in p.
    hoop, _ = s.thin_wall_vessel_stresses(1.0, 24 * INCH, 0.5 * INCH)
    sphere, _ = s.thin_wall_vessel_stresses(1.0, 24 * INCH, 0.5 * INCH, "sphere")
    assert abs(20e3 / hoop - 417) < 0.5
    assert abs(20e3 / sphere - 833) < 0.5


def test_hibbeler_example_9_2_stresses_on_a_rotated_element():
    # sigma_x = -80 MPa, sigma_y = 50 MPa, tau_xy = -25 MPa, element rotated 30 deg clockwise:
    # sigma_x' = -25.8 MPa, sigma_y' = -4.15 MPa, tau_x'y' = -68.8 MPa.
    sx, sy, txy = s.stress_transformation(-80e6, 50e6, -25e6, np.radians(-30))
    assert abs(sx / 1e6 + 25.8) < 0.05
    assert abs(sy / 1e6 + 4.15) < 0.005
    assert abs(txy / 1e6 + 68.8) < 0.05


def test_hibbeler_example_13_1_w200x46_pinned_column():
    # A992 W200x46, pinned, 4 m long, I_y = 15.3e6 mm^4, A = 5890 mm^2, E = 200 GPa,
    # sigma_Y = 345 MPa: buckling governs at P_cr = 1887.6 kN, sigma_cr = 320 MPa.
    P = s.euler_buckling_load(200e9, 15.3e-6, 4.0)
    sigma = P / 5890e-6
    assert abs(P / 1e6 - 1.89) < 0.005
    assert abs(sigma / 1e6 - 320) < 0.5
    assert sigma < 345e6


def test_gere_example_3_1_solid_bar_in_torsion():
    # d = 1.5 in, L = 54 in, G = 11.5e6 psi, T = 250 lb-ft: tau_max = 4530 psi,
    # phi = 0.02834 rad = 1.62 deg. With tau_allow = 6000 psi and phi_allow = 2.5 deg the
    # torques are 3980 and 4620 lb-in, so T_max = 3980 lb-in.
    d, L, G = 1.5 * INCH, 54 * INCH, 11.5e6 * PSI
    J = np.pi * d**4 / 32
    T = 250 * 12 * LBF * INCH
    tau = s.torsion_shear_stress(T, d / 2, J)
    phi = s.torsion_angle(T, L, G, J)
    assert abs(tau / PSI - 4530) < 5
    assert abs(phi - 0.02834) < 0.000005
    assert abs(np.degrees(phi) - 1.62) < 0.005
    lb_in = LBF * INCH
    T_stress = T * (6000 * PSI) / tau
    T_twist = T * np.radians(2.5) / phi
    assert abs(T_stress / lb_in - 3980) < 5
    assert abs(T_twist / lb_in - 4620) < 5
    assert abs(min(T_stress, T_twist) / lb_in - 3980) < 5


def test_gere_example_7_2_element_inclined_at_45_degrees():
    # sigma_x = 16,000 psi, sigma_y = 6,000 psi, tau_xy = 4,000 psi, theta = 45 deg:
    # sigma_x1 = 15,000 psi, sigma_y1 = 7,000 psi, tau_x1y1 = -5,000 psi.
    sx, sy, txy = s.stress_transformation(16000 * PSI, 6000 * PSI, 4000 * PSI, np.radians(45))
    assert abs(sx / PSI - 15000) < 500
    assert abs(sy / PSI - 7000) < 500
    assert abs(txy / PSI + 5000) < 500


def test_gere_example_7_3_principal_stresses():
    # sigma_x = 12,300 psi, sigma_y = -4,200 psi, tau_xy = -4,700 psi: sigma_1 = 13,500 psi at
    # theta_p1 = -14.8 deg, sigma_2 = -5,440 psi, tau_max = 9,490 psi.
    _, radius, s1, s2 = s.mohr_circle(12300 * PSI, -4200 * PSI, -4700 * PSI)
    assert abs(s1 / PSI - 13500) < 50
    assert abs(s2 / PSI + 5440) < 5
    assert abs(radius / PSI - 9490) < 5
    theta_p1 = np.radians(-14.8)
    sx1, _, txy1 = s.stress_transformation(12300 * PSI, -4200 * PSI, -4700 * PSI, theta_p1)
    assert abs(sx1 / PSI - 13500) < 50
    assert abs(txy1) < 0.005 * radius

"""Textbook reference problems (ROADMAP B4.fluid_engineering).

Each test works a published example with the package's own functions and
checks the book's answer to the precision it prints: a printed 4.84 is
met within 0.005. The books round intermediate values and read friction
factors off the Moody chart, which the guide's "Reference problems"
section covers where it matters.

F. M. White, *Fluid Mechanics*, 5th ed. (New York: McGraw-Hill, 2003).
B. R. Munson, D. F. Young, T. H. Okiishi and W. W. Huebsch,
*Fundamentals of Fluid Mechanics*, 6th ed. (Hoboken: Wiley, 2009).
"""

import numpy as np
from scipy.optimize import brentq

from engineeringkit import fluid_engineering as f


def _pipe_head_loss(Q, d, L, nu, roughness):
    v = Q / (np.pi * d * d / 4)
    factor = f.colebrook_friction_factor(v * d / nu, roughness / d)
    return f.darcy_head_loss(factor, L, d, v)


def test_white_example_6_7_oil_in_sloping_cast_iron_pipe():
    # Oil, rho = 900 kg/m^3 and nu = 1e-5 m^2/s, at 0.2 m^3/s through 500 m of
    # 200-mm cast-iron pipe (epsilon = 0.26 mm) sloping down 10 degrees.
    rho, nu, Q, L, d = 900.0, 1e-5, 0.2, 500.0, 0.2
    v = Q / (np.pi * d * d / 4)
    assert abs(v - 6.37) < 0.005
    assert abs(v * d / nu - 127_000) < 500
    head_loss = _pipe_head_loss(Q, d, L, nu, 0.26e-3)
    assert abs(head_loss - 117) < 0.5
    # The book subtracts its rounded 117 m from the 86.8 m drop, with a
    # chart-read f, and prints about 265 kPa; Colebrook's f gives 270 kPa.
    delta_p = rho * f.G * (head_loss - L * np.sin(np.radians(10)))
    assert abs(delta_p - 270e3) < 1e3


def test_white_example_6_9_flow_rate_from_head_loss():
    # Oil, nu = 2e-5 m^2/s, in a 30-cm pipe 100 m long with an 8 m head loss
    # and epsilon/d = 0.0002: find the velocity and the flow rate.
    nu, L, d, head_loss = 2e-5, 100.0, 0.3, 8.0
    Q = brentq(lambda Q: _pipe_head_loss(Q, d, L, nu, 0.0002 * d) - head_loss, 0.01, 2.0)
    assert abs(Q / (np.pi * d * d / 4) - 4.84) < 0.005
    assert abs(Q - 0.342) < 0.0005


def test_white_example_6_10_diameter_from_head_loss():
    # Example 6.9 worked backward: Q = 0.342 m^3/s and epsilon = 0.06 mm are
    # known, the diameter is not.
    nu, L, Q, head_loss = 2e-5, 100.0, 0.342, 8.0
    d = brentq(lambda d: _pipe_head_loss(Q, d, L, nu, 0.06e-3) - head_loss, 0.1, 1.0)
    assert abs(d - 0.300) < 0.0005


def test_white_example_10_8_hydraulic_jump():
    # A jump in a wide channel at q = 10 m^3/(s m) from y1 = 1.25 m.
    q, y1 = 10.0, 1.25
    v1 = q / y1
    assert abs(f.froude_number(v1, y1) - 2.285) < 0.0005
    y2 = f.sequent_depth(y1, v1)
    v2 = q / y2
    assert abs(y2 - 3.46) < 0.005
    assert abs(v2 - 2.89) < 0.005
    assert abs(f.froude_number(v2, y2) - 0.496) < 0.0005
    head_loss = f.hydraulic_jump_head_loss(y1, y2)
    assert abs(head_loss - 0.625) < 0.0005
    assert abs(100 * head_loss / f.specific_energy(y1, q) - 14) < 0.5
    power = 998.0 * f.G * q * head_loss
    assert abs(power - 61e3) < 0.5e3
    assert abs(power / (q * 998.0 * 4200) - 0.0015) < 0.00005


def test_munson_example_2_6_inclined_circular_gate():
    # A 4-m circular gate in a wall inclined at 60 degrees, its horizontal
    # shaft 10 m below the surface, gamma = 9.80 kN/m^3.
    area, second_moment = np.pi * 2.0**2, np.pi * 2.0**4 / 4
    force = f.hydrostatic_force_plane(1000.0, 10.0, area, g=9.80)
    assert abs(force - 1.23e6) < 0.005e6
    y_c = 10.0 / np.sin(np.radians(60))
    y_r = f.center_of_pressure_depth(y_c, area, second_moment)
    assert abs(y_r - 11.6) < 0.05
    assert abs(force * (y_r - y_c) - 1.07e5) < 0.005e5

import numpy as np
import pytest
from scipy.integrate import quad

from engineeringkit import fluid_engineering as f


def test_poiseuille_profile_integrates_to_volume_flow():
    dp, r, L, mu = 100.0, 0.002, 1.0, 0.001
    profile_flow = quad(lambda x: 2 * np.pi * x * dp * (r * r - x * x) / (4 * mu * L), 0, r)[0]
    assert np.isclose(f.poiseuille_flow(dp, r, L, mu), profile_flow)
    v = profile_flow / (np.pi * r * r)
    Re = f.reynolds_number(1000, v, 2 * r, mu)
    loss = f.darcy_head_loss(f.laminar_friction_factor(Re), L, 2 * r, v)
    assert np.isclose(1000 * f.G * loss, dp)


def test_bernoulli_venturi_recovery():
    rho, A1, A2, Q = 1000.0, 0.02, 0.01, 0.02
    v1, v2 = Q / A1, Q / A2
    p2 = f.bernoulli_pressure(1e5, rho, v1, v2)
    assert np.isclose(f.venturi_flow(1e5 - p2, rho, A1, A2), Q)
    assert np.isclose(f.continuity_velocity(A1, v1, A2), v2)


@pytest.mark.parametrize("Re,rough", [(4000, 0), (1e5, 1e-4), (1e8, 0.05)])
def test_colebrook_satisfies_implicit_equation(Re, rough):
    factor = f.colebrook_friction_factor(Re, rough)
    residual = 1 / np.sqrt(factor) + 2 * np.log10(rough / 3.7 + 2.51 / (Re * np.sqrt(factor)))
    assert abs(residual) < 1e-8


def test_pump_power_and_hydrostatic_invariants():
    assert np.isclose(f.pump_input_power(1000, 0.01, 10, 0.8) * 0.8, f.pump_hydraulic_power(1000, 0.01, 10))
    assert f.buoyancy_force(1000, 0.01) == f.hydrostatic_pressure(1000, 1) * 0.01


def test_transition_not_silently_turbulent():
    with pytest.raises(ValueError, match="transition"):
        f.colebrook_friction_factor(3000, 1e-4)

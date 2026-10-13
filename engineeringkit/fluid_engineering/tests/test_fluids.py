import numpy as np

from engineeringkit.fluid_engineering import *


def test_reynolds():
    assert reynolds_number(1000, 1, 0.1, 0.001) == 100000


def test_laminar_friction():
    assert np.isclose(laminar_friction_factor(1000), 0.064)


def test_poiseuille_scaling():
    assert np.isclose(poiseuille_flow(2, 1, 1, 1), 2 * poiseuille_flow(1, 1, 1, 1))


def test_heads():
    assert np.isclose(velocity_head(np.sqrt(2 * G)), 1)


def test_dimensionless():
    assert np.isclose(froude_number(np.sqrt(G), 1), 1)
    assert weber_number(1000, 1, 0.01, 0.1) == 100


def test_water_hammer():
    assert water_hammer_pressure(1000, 1200, 1) == 1.2e6

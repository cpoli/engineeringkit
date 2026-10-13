import numpy as np

from engineeringkit.thermodynamics import *


def test_carnot():
    assert np.isclose(carnot_efficiency(600, 300), 0.5)


def test_otto():
    assert 0 < otto_efficiency(10) < 1


def test_isothermal_work():
    assert np.isclose(reversible_isothermal_work_ideal_gas(1, 300, 1, 2), 8.31446261815324 * 300 * np.log(2))


def test_clausius_and_exergy():
    assert clausius_inequality([-2, 1])
    assert exergy_destruction(300, 2) == 600


def test_rankine_positive():
    assert 0 < rankine_ideal_efficiency(200, 3200, 2400, 190) < 1

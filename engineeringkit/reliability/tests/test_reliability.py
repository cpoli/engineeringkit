import numpy as np

import engineeringkit.reliability as m


def test_system_reliability():
    assert np.isclose(m.series_reliability(0.9, 0.8), 0.72)
    assert np.isclose(m.parallel_reliability(0.9, 0.8), 0.98)


def test_weibull_exponential_limit():
    t = 2.0
    assert np.isclose(m.weibull_reliability(t, eta=5, beta=1), m.exponential_reliability(t, 1 / 5))

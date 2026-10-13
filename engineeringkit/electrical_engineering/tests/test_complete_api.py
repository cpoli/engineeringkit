import numpy as np
import pytest

from engineeringkit import electrical_engineering as e


def test_dc_sources_and_ac_power():
    assert e.ohms_law(current=2, resistance=5) == 10
    assert e.ohms_law(voltage=10, current=2) == 5
    assert np.isclose(e.rms_sine(np.sqrt(2)), 1)
    assert e.real_power(10, 2, 0.8) == 16
    assert np.isclose(e.three_phase_power(10, 2, 0.8), 16 * np.sqrt(3))
    assert e.transformer_secondary(240, 100, 50) == 120
    assert e.thevenin_from_source(10, 5) == (10, 5)
    assert e.norton_from_thevenin(10, 5) == (2, 5)
    assert e.wheatstone_balanced(1, 2, 2, 4)
    assert not e.wheatstone_balanced(1, 2, 3, 4)
    assert e.rl_time_constant(2, 4) == 0.5
    assert e.power_factor(8, 10) == 0.8
    assert np.isclose(e.reactive_power(10, 2, np.pi / 2), 20)


@pytest.mark.parametrize("kwargs", [{}, {"voltage": 1}, {"voltage": 1, "current": 1, "resistance": 1}, {"voltage": 1, "current": 0}])
def test_ohm_incomplete_or_degenerate(kwargs):
    with pytest.raises(ValueError):
        e.ohms_law(**kwargs)


def test_circuit_plot_and_summary():
    import matplotlib.pyplot as plt

    r = e.SeriesRLC(2, 1, 0.25).step_response()
    assert r.summary()["peak"] == r.response.max()
    fig, ax = r.plot()
    assert ax.get_ylabel() == "current [A]"
    plt.close(fig)
    with pytest.raises(ValueError):
        e.SeriesRLC(2, 1, 0.25).step_response([0, 0])

import numpy as np
import pytest

from engineeringkit import heat_transfer as h


def test_resistance_conductance_and_convection():
    assert np.isclose(h.sphere_resistance(1, 2, 1), 1 / (8 * np.pi))
    assert h.convection_heat_rate(10, 2, 400, 300) == 2000
    assert h.convection_resistance(10, 2) == 0.05
    assert h.overall_u(0.1, 0.1) == 5
    assert h.critical_insulation_radius_cylinder(2, 10) == 0.2
    assert h.critical_insulation_radius_sphere(2, 10) == 0.4
    assert h.biot_number(10, 0.01, 1) == 0.1
    assert h.fourier_number(0.01, 2, 0.1) == pytest.approx(2)


def test_radiation_fin_and_diffusion_boundary():
    hr = h.radiation_coefficient(0.8, 400, 300)
    assert np.isclose(hr * 2 * 100, h.radiation_heat_rate(0.8, 2, 400, 300))
    assert h.semi_infinite_temperature(400, 300, 0, 0.01, 1) == 300
    assert h.fin_parameter(1, 1, 1, 1) == 1
    assert np.isclose(h.fin_heat_rate_adiabatic(1, 1, 1, 1, 1, 400, 300), 100 * np.tanh(1))
    assert h.heat_exchanger_rate(10, 2, 20) == 400
    assert np.isclose(h.blackbody_peak_wavelength(1000), 2.897771955e-6)
    assert np.isclose(h.stefan_boltzmann_emissive_power(300), 5.670374419e-8 * 300**4)


def test_wall_plot_summary_and_empty_layers():
    import matplotlib.pyplot as plt

    r = h.composite_wall(400, 300, [(0.1, 1)])
    assert r.summary() == {"heat_rate_W": 1000, "interfaces_K": [400, 300]}
    fig, ax = r.plot()
    assert ax.get_ylabel() == "temperature [K]"
    plt.close(fig)
    with pytest.raises(ValueError):
        h.composite_wall(400, 300, [])

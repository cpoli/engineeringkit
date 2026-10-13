import numpy as np
import pytest

from engineeringkit.heat_transfer import fin_efficiency_adiabatic


def test_fin_limit():
    assert fin_efficiency_adiabatic(0, 1) == 1
    assert fin_efficiency_adiabatic(1, 0) == 1
    assert np.isclose(fin_efficiency_adiabatic(1e-8, 1), 1)
    with pytest.raises(ValueError):
        fin_efficiency_adiabatic(-1, 1)

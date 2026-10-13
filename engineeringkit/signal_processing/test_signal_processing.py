import numpy as np

import engineeringkit.signal_processing as m


def test_fft_roundtrip():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert np.allclose(m.idft(m.dft(x)).real, x)


def test_alias_and_nyquist():
    assert m.nyquist_frequency(100) == 50
    assert np.isclose(m.alias_frequency(70, 100), 30)

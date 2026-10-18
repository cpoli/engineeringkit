import numpy as np
import pytest

from engineeringkit import signal_processing as s


def test_discrete_transforms_and_filters():
    assert s.sampling_period(10) == 0.1
    assert s.db20(10) == 20
    np.testing.assert_allclose(s.convolution([1, 2], [1, 1]), [1, 3, 2])
    np.testing.assert_allclose(s.correlation([1, 2], [1, 2]), [2, 5, 2])
    np.testing.assert_allclose(s.moving_average([1, 2, 3], 2), [1.5, 2.5])
    np.testing.assert_allclose(s.fft_frequencies(4, 4), [0, 1, -2, -1])
    f, magnitude = s.windowed_spectrum(np.ones(4), 4, window="boxcar")
    np.testing.assert_allclose(f, [0, 1, 2])
    np.testing.assert_allclose(magnitude, [4, 0, 0])
    b, a = s.butterworth_highpass(1, 1, 10)
    _, gain = s.frequency_response(b, a, worN=[0], fs=10)
    assert abs(gain[0]) < 1e-15
    np.testing.assert_allclose(s.filter_signal([1, 1], [1], [1, 0, 0]), [1, 1, 0])
    with pytest.raises(ValueError):
        s.quantization_step(1, 1, 8)
    with pytest.raises(ValueError):
        s.quantize([np.nan], 0, 1, 8)


def test_spectrum_plot():
    import matplotlib.pyplot as plt

    spec = s.Spectrum.from_signal([1, 1, 1, 1], 4)
    np.testing.assert_allclose(spec.magnitude, [4, 0, 0])
    fig, ax = spec.plot()
    assert ax.get_xlabel() == "frequency [Hz]"
    plt.close(fig)

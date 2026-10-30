"""Result containers for engineeringkit.signal_processing.

Most models return plain arrays. A one-sided spectrum, a short-time Fourier
transform, an adaptive filter run and an autoregressive model carry several
related arrays, so each gets a frozen dataclass with ``summary()`` and
``plot()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive

__all__ = ["Spectrum", "STFTResult", "LMSResult", "ARModel"]


@dataclass(frozen=True)
class Spectrum:
    """One-sided, unnormalized FFT magnitude of a real signal."""

    frequency: np.ndarray
    """ndarray: Bin frequencies [Hz], from 0 to the Nyquist frequency."""

    magnitude: np.ndarray
    """ndarray: :math:`|X_k|` of the unnormalized real FFT [signal units]."""

    @classmethod
    def from_signal(cls, x, fs):
        """Magnitude spectrum of real samples.

        Parameters
        ----------
        x : array_like
            Real, finite samples.
        fs : float
            Sampling rate [Hz], positive.

        Returns
        -------
        Spectrum
            One-sided frequencies [Hz] and unnormalized FFT magnitudes.
        """
        _positive(fs=fs)
        _finite(x=x)
        x = np.asarray(x)
        return cls(np.fft.rfftfreq(len(x), 1 / fs), np.abs(np.fft.rfft(x)))

    def summary(self) -> dict:
        """Number of bins, bin spacing [Hz] and frequency of the largest bin [Hz]."""
        spacing = float(self.frequency[1] - self.frequency[0]) if self.frequency.size > 1 else 0.0
        return {"bins": int(self.frequency.size), "resolution": spacing, "peak_frequency": float(self.frequency[np.argmax(self.magnitude)])}

    def plot(self, ax=None):
        """Magnitude against frequency; returns ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots() if ax is None else (ax.figure, ax)
        ax.plot(self.frequency, self.magnitude)
        ax.set(xlabel="frequency [Hz]", ylabel="magnitude")
        return fig, ax


@dataclass(frozen=True)
class STFTResult:
    """Output of :func:`engineeringkit.signal_processing.stft`."""

    time: np.ndarray
    """ndarray: Centre time of each frame [s]."""

    frequency: np.ndarray
    """ndarray: One-sided bin frequencies [Hz]."""

    values: np.ndarray
    """ndarray of complex: Coefficients, shape ``(frequency.size, time.size)`` [signal units]."""

    def summary(self) -> dict:
        """Frame and bin counts, and the time step and frequency resolution of the grid."""
        return {
            "frames": int(self.time.size),
            "bins": int(self.frequency.size),
            "time_step": float(self.time[1] - self.time[0]) if self.time.size > 1 else 0.0,
            "frequency_resolution": float(self.frequency[1] - self.frequency[0]),
        }

    def plot(self, ax=None):
        """Spectrogram :math:`|X|^2` in decibels relative to its peak; returns ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots() if ax is None else (ax.figure, ax)
        level = np.abs(self.values) ** 2
        level_db = 10 * np.log10(np.maximum(level / level.max(), 1e-12))
        mesh = ax.pcolormesh(self.time, self.frequency, level_db, shading="nearest", vmin=-60, vmax=0)
        fig.colorbar(mesh, ax=ax, label="level [dB]")
        ax.set(xlabel="time [s]", ylabel="frequency [Hz]")
        return fig, ax


@dataclass(frozen=True)
class LMSResult:
    """Output of :func:`engineeringkit.signal_processing.lms_filter`."""

    output: np.ndarray
    """ndarray: Filter output :math:`y_n = w_n^T u_n` [signal units]."""

    error: np.ndarray
    """ndarray: Error :math:`e_n = d_n - y_n` [signal units]."""

    weights: np.ndarray
    """ndarray: Weight vector after the last sample [-]."""

    weight_history: np.ndarray
    """ndarray: Weights before each update, shape ``(n_samples, n_taps)`` [-]."""

    def summary(self) -> dict:
        """Number of taps and mean squared error over the first and last tenth of the run."""
        tenth = max(1, self.error.size // 10)
        return {
            "taps": int(self.weights.size),
            "initial_mse": float(np.mean(np.abs(self.error[:tenth]) ** 2)),
            "final_mse": float(np.mean(np.abs(self.error[-tenth:]) ** 2)),
        }

    def plot(self, ax=None):
        """Learning curve, the squared error against sample index on a log axis; returns ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots() if ax is None else (ax.figure, ax)
        ax.semilogy(np.maximum(np.abs(self.error) ** 2, 1e-300))
        ax.set(xlabel="sample", ylabel="squared error")
        return fig, ax


@dataclass(frozen=True)
class ARModel:
    r"""Autoregressive model :math:`x_n + \sum_{k=1}^{p} a_k x_{n-k} = e_n`.

    Returned by :func:`engineeringkit.signal_processing.levinson_durbin` and
    :func:`engineeringkit.signal_processing.yule_walker`.
    """

    coefficients: np.ndarray
    """ndarray: :math:`[1, a_1, \\ldots, a_p]`, the prediction-error filter :math:`A(z)` [-]."""

    noise_variance: float
    """float: Variance :math:`\\sigma^2 = E_p` of the white input :math:`e_n` [signal units squared]."""

    reflection_coefficients: np.ndarray
    """ndarray: :math:`k_1..k_p`, all inside (-1, 1) for a stable model [-]."""

    prediction_errors: np.ndarray
    """ndarray: Mean squared prediction error :math:`E_0..E_p` of each order [signal units squared]."""

    @property
    def order(self) -> int:
        """int: Model order :math:`p`."""
        return int(self.coefficients.size - 1)

    def psd(self, fs=1.0, n=512):
        r"""One-sided power spectral density :math:`P(f) = 2\sigma^2/(f_s|A(e^{j2\pi f/f_s})|^2)`.

        Parameters
        ----------
        fs : float, optional
            Sampling rate [Hz], positive.
        n : int, optional
            Number of frequencies from 0 to :math:`f_s/2`, at least 2.

        Returns
        -------
        tuple of ndarray
            ``(frequency, psd)`` [Hz] and [signal units squared per Hz]; its
            integral over :math:`[0, f_s/2]` is the model's variance.
        """
        _positive(fs=fs)
        if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2:
            raise ValueError("n must be an integer >= 2")
        f = np.linspace(0, fs / 2, n)
        A = np.exp(-2j * np.pi * np.outer(f / fs, np.arange(self.coefficients.size))) @ self.coefficients
        return f, 2 * self.noise_variance / (fs * np.abs(A) ** 2)

    def summary(self) -> dict:
        """Order, noise variance, prediction gain :math:`E_0/E_p` and whether every :math:`|k_m| < 1`."""
        return {
            "order": self.order,
            "noise_variance": float(self.noise_variance),
            "prediction_gain": float(self.prediction_errors[0] / self.prediction_errors[-1]),
            "stable": bool(np.all(np.abs(self.reflection_coefficients) < 1)),
        }

    def plot(self, ax=None, fs=1.0):
        """Model spectrum in decibels against frequency; returns ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots() if ax is None else (ax.figure, ax)
        f, P = self.psd(fs)
        ax.plot(f, 10 * np.log10(P))
        ax.set(xlabel="frequency [Hz]", ylabel="PSD [dB]")
        return fig, ax

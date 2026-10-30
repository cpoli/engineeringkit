r"""Thermal noise, the noise figure of cascaded amplifiers, and the noise performance of AM and FM.

Johnson measured, and Nyquist explained, the noise power :math:`kTB`
delivered by any resistor; Friis showed that in a chain of amplifiers the
first stage dominates the total noise factor. Crosby (1937) showed that
above a threshold wideband FM improves the output signal-to-noise ratio by
:math:`3\beta^2` times the message power, while envelope-detected AM loses
the power spent on its carrier.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.communications.utils.decibels import db

__all__ = [
    "BOLTZMANN",
    "thermal_noise_power",
    "noise_figure",
    "cascaded_noise_factor",
    "baseband_snr",
    "am_output_snr",
    "fm_output_snr",
    "fm_threshold_snr",
]

BOLTZMANN = 1.380649e-23
"""Boltzmann constant [J/K] (exact)."""


def thermal_noise_power(B, T=290.0, k=BOLTZMANN):
    r"""Available thermal (Johnson-Nyquist) noise power :math:`N = kTB`.

    Parameters
    ----------
    B : float or array_like
        Noise bandwidth [Hz], nonnegative.
    T : float or array_like, optional
        Absolute noise temperature [K], nonnegative; 290 K is the standard
        reference.
    k : float, optional
        Boltzmann constant [J/K].

    Returns
    -------
    float or ndarray
        Noise power [W] (-174 dBm in 1 Hz at 290 K).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    J. B. Johnson, "Thermal agitation of electricity in conductors," *Phys.
    Rev.* 32 (1928) 97-109; H. Nyquist, *Phys. Rev.* 32 (1928) 110-113.

    Examples
    --------
    >>> round(float(10 * np.log10(thermal_noise_power(1.0) / 1e-3)), 2)
    -173.98
    """
    _nonnegative(B=B, T=T)
    _positive(k=k)
    return k * np.asarray(T) * B


def noise_figure(F):
    r"""Noise figure :math:`NF = 10\log_{10} F` from a linear noise factor.

    Parameters
    ----------
    F : float or array_like
        Noise factor [-], at least one (one for a noiseless device).

    Returns
    -------
    float or ndarray
        Noise figure [dB].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(noise_figure(2.0)) > 3.0
    True
    """
    _interval("noise factor", F, 1, np.inf, closed="left")
    return db(F)


def cascaded_noise_factor(gains, noise_factors):
    r"""Friis's formula for the noise factor of amplifiers in cascade.

    .. math::

        F = F_1 + \frac{F_2 - 1}{G_1} + \frac{F_3 - 1}{G_1 G_2} + \cdots

    Parameters
    ----------
    gains : array_like
        Linear available power gain of each stage [-], positive. The last
        stage's gain does not affect the result.
    noise_factors : array_like
        Linear noise factor of each stage [-], at least one.

    Returns
    -------
    float
        Overall noise factor [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    H. T. Friis, "Noise figures of radio receivers," *Proc. IRE* 32 (1944)
    419-422.

    Examples
    --------
    >>> round(cascaded_noise_factor([10.0, 10.0], [2.0, 3.0]), 6)
    2.2
    """
    gains = _vector("gains", gains)
    noise_factors = _vector("noise_factors", noise_factors)
    if gains.shape != noise_factors.shape:
        raise ValueError("one gain and noise factor are required per stage")
    _positive(gains=gains)
    _interval("noise factors", noise_factors, 1, np.inf, closed="left")
    upstream_gain = np.concatenate(([1.0], np.cumprod(gains[:-1])))
    return float(noise_factors[0] + np.sum((noise_factors[1:] - 1) / upstream_gain[1:]))


def baseband_snr(received_power, N0, message_bandwidth):
    r"""Reference SNR :math:`\gamma = P_R/(N_0 W)` against which analog modulations are compared.

    It is the output SNR of baseband transmission, and also of coherent
    DSB-SC and SSB, with the same received power.

    Parameters
    ----------
    received_power : float or array_like
        Received signal power :math:`P_R` [W], nonnegative.
    N0 : float or array_like
        One-sided noise power spectral density [W/Hz], positive.
    message_bandwidth : float or array_like
        Message bandwidth :math:`W` [Hz], positive.

    Returns
    -------
    float or ndarray
        :math:`\gamma` [-] (linear).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(baseband_snr(1e-9, 1e-15, 1e4))
    100.0
    """
    _nonnegative(received_power=received_power)
    _positive(N0=N0, message_bandwidth=message_bandwidth)
    return np.asarray(received_power, dtype=float) / (np.asarray(N0, dtype=float) * message_bandwidth)


def am_output_snr(snr_baseband, mu=1.0, message_power=0.5):
    r"""Output SNR of conventional AM with an envelope detector above threshold, :math:`\dfrac{\mu^2 P_m}{1 + \mu^2 P_m}\gamma`.

    Parameters
    ----------
    snr_baseband : float or array_like
        Reference SNR :math:`\gamma` from :func:`baseband_snr` [-], nonnegative.
    mu : float or array_like, optional
        Modulation index [-] in ``[0, 1]``.
    message_power : float or array_like, optional
        Mean square of the normalised message, :math:`P_m` [-] in ``[0, 1]``;
        0.5 for a full-scale tone.

    Returns
    -------
    float or ndarray
        Output SNR [-] (linear): at best one third of :math:`\gamma`, for a
        tone at :math:`\mu = 1`.

    Notes
    -----
    Valid when the carrier is well above the noise; below that the envelope
    detector's output collapses. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/communications`.

    Examples
    --------
    >>> round(float(am_output_snr(300.0)), 6)
    100.0
    """
    _nonnegative(snr_baseband=snr_baseband)
    _interval("mu", mu, 0, 1)
    _interval("message_power", message_power, 0, 1)
    efficiency = np.asarray(mu, dtype=float) ** 2 * message_power
    return efficiency / (1 + efficiency) * np.asarray(snr_baseband, dtype=float)


def fm_output_snr(snr_baseband, beta, message_power=0.5):
    r"""Output SNR of FM with a limiter-discriminator above threshold, :math:`3\beta^2 P_m\gamma`.

    Parameters
    ----------
    snr_baseband : float or array_like
        Reference SNR :math:`\gamma` from :func:`baseband_snr` [-], nonnegative.
    beta : float or array_like
        Modulation index :math:`\beta = \Delta f/W` [-], nonnegative.
    message_power : float or array_like, optional
        Mean square of the message normalised to :math:`|m| \le 1`,
        :math:`P_m` [-] in ``[0, 1]``; 0.5 for a full-scale tone.

    Returns
    -------
    float or ndarray
        Output SNR [-] (linear). The improvement grows as the square of the
        deviation, so FM trades bandwidth (:func:`carson_bandwidth`) for
        noise immunity. No pre-emphasis is assumed.

    Notes
    -----
    Valid only above the threshold of :func:`fm_threshold_snr`. Inputs are
    raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    M. G. Crosby, "Frequency modulation noise characteristics," *Proc. IRE*
    25 (1937) 472-514.

    Examples
    --------
    >>> float(fm_output_snr(100.0, 5.0))
    3750.0
    """
    _nonnegative(snr_baseband=snr_baseband, beta=beta)
    _interval("message_power", message_power, 0, 1)
    return 3 * np.asarray(beta, dtype=float) ** 2 * message_power * np.asarray(snr_baseband, dtype=float)


def fm_threshold_snr(beta, carrier_to_noise=10.0):
    r"""Reference SNR :math:`\gamma` at the FM threshold, :math:`2(\beta + 1)\,\mathrm{CNR}_{th}`.

    Below threshold, noise spikes ("clicks") from phase wrap-arounds swamp
    the discriminator output and :func:`fm_output_snr` no longer holds.

    Parameters
    ----------
    beta : float or array_like
        Modulation index [-], nonnegative.
    carrier_to_noise : float, optional
        Carrier-to-noise ratio in the Carson bandwidth at threshold [-]
        (linear), positive; about 10 (10 dB) is customary.

    Returns
    -------
    float or ndarray
        Threshold value of :math:`\gamma = P_R/(N_0 W)` [-] (linear).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> float(fm_threshold_snr(5.0))
    120.0
    """
    _nonnegative(beta=beta)
    _positive(carrier_to_noise=carrier_to_noise)
    return 2 * (np.asarray(beta, dtype=float) + 1) * carrier_to_noise

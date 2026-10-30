r"""Rayleigh fading and diversity combining.

Rayleigh (1880) showed that the sum of many vibrations of equal frequency
and random phase has an amplitude with density :math:`2r/\Omega\,e^{-r^2/\Omega}`.
A radio signal reaching a receiver by many scattered paths fades the same
way, and the bit-error rate falls only as :math:`1/\bar\gamma` instead of
exponentially. Combining several independently fading branches (Kahn's
maximal-ratio combiner, Brennan's comparison of combiners) restores a steep
error curve, with slope :math:`\bar\gamma^{-L}` for :math:`L` branches.
"""

from __future__ import annotations

import numpy as np
from scipy import special

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "rayleigh_pdf",
    "ber_bpsk_rayleigh",
    "ber_bpsk_mrc_rayleigh",
    "selection_outage_probability",
    "maximal_ratio_combine",
]


def rayleigh_pdf(r, mean_power=1.0):
    r"""Rayleigh density of a fading envelope, :math:`f(r) = \dfrac{2r}{\Omega}e^{-r^2/\Omega}`.

    Parameters
    ----------
    r : float or array_like
        Envelope amplitude [V or any amplitude unit], nonnegative.
    mean_power : float, optional
        Mean square envelope :math:`\Omega = E[r^2]` [amplitude unit squared],
        positive.

    Returns
    -------
    float or ndarray
        Probability density [1/amplitude unit]. The power :math:`r^2` is
        exponentially distributed with mean :math:`\Omega`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    References
    ----------
    Lord Rayleigh, "On the resultant of a large number of vibrations of the
    same pitch and of arbitrary phase," *Phil. Mag.* 10 (1880) 73-78.

    Examples
    --------
    >>> round(float(rayleigh_pdf(np.sqrt(0.5))), 6)
    0.857764
    """
    _nonnegative(r=r)
    _positive(mean_power=mean_power)
    r = np.asarray(r, dtype=float)
    return 2 * r / mean_power * np.exp(-(r**2) / mean_power)


def _mu(mean_ebn0):
    _nonnegative(mean_ebn0=mean_ebn0)
    g = np.asarray(mean_ebn0, dtype=float)
    return np.sqrt(g / (1 + g))


def ber_bpsk_rayleigh(mean_ebn0):
    r"""Average bit-error rate of coherent BPSK in flat Rayleigh fading, :math:`\tfrac12\left(1 - \sqrt{\bar\gamma/(1+\bar\gamma)}\right)`.

    Parameters
    ----------
    mean_ebn0 : float or array_like
        Mean :math:`E_b/N_0` [-] (linear), nonnegative.

    Returns
    -------
    float or ndarray
        Bit-error probability [-], about :math:`1/(4\bar\gamma)` at high SNR.

    Notes
    -----
    The channel gain is assumed known at the receiver and constant over a
    bit. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> round(float(ber_bpsk_rayleigh(100.0)), 6)
    0.002481
    """
    return 0.5 * (1 - _mu(mean_ebn0))


def ber_bpsk_mrc_rayleigh(mean_ebn0, branches):
    r"""Average BPSK bit-error rate with ``L``-branch maximal-ratio combining in independent Rayleigh fading.

    .. math::

        P_b = \left(\frac{1-\mu}{2}\right)^L \sum_{k=0}^{L-1}\binom{L-1+k}{k}\left(\frac{1+\mu}{2}\right)^k,
        \qquad \mu = \sqrt{\frac{\bar\gamma}{1+\bar\gamma}}

    Parameters
    ----------
    mean_ebn0 : float or array_like
        Mean :math:`E_b/N_0` per branch [-] (linear), nonnegative.
    branches : int
        Number of independently fading branches :math:`L`, at least 1.

    Returns
    -------
    float or ndarray
        Bit-error probability [-], falling as :math:`\bar\gamma^{-L}`.

    Notes
    -----
    ``branches=1`` is :func:`ber_bpsk_rayleigh`. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/communications`.

    References
    ----------
    D. G. Brennan, "Linear diversity combining techniques," *Proc. IRE* 47
    (1959) 1075-1102.

    Examples
    --------
    >>> bool(np.isclose(ber_bpsk_mrc_rayleigh(10.0, 1), ber_bpsk_rayleigh(10.0)))
    True
    """
    _integer("branches", branches, 1)
    mu = _mu(mean_ebn0)
    total = sum(special.comb(branches - 1 + k, k) * ((1 + mu) / 2) ** k for k in range(branches))
    return ((1 - mu) / 2) ** branches * total


def selection_outage_probability(threshold_snr, mean_snr, branches=1):
    r"""Probability that the best of ``L`` Rayleigh branches is below a threshold, :math:`(1 - e^{-\gamma_t/\bar\gamma})^L`.

    Parameters
    ----------
    threshold_snr : float or array_like
        SNR the receiver needs, :math:`\gamma_t` [-] (linear), nonnegative.
    mean_snr : float or array_like
        Mean SNR of each branch :math:`\bar\gamma` [-] (linear), positive.
    branches : int, optional
        Number of independently fading branches :math:`L` for selection
        combining, at least 1.

    Returns
    -------
    float or ndarray
        Outage probability [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/communications`.

    Examples
    --------
    >>> round(float(selection_outage_probability(1.0, 10.0, 2)), 6)
    0.009056
    """
    _nonnegative(threshold_snr=threshold_snr)
    _positive(mean_snr=mean_snr)
    _integer("branches", branches, 1)
    return (-np.expm1(-np.asarray(threshold_snr, dtype=float) / mean_snr)) ** branches


def maximal_ratio_combine(received, gains):
    r"""Maximal-ratio combining of diversity branches, :math:`\hat s = \sum_l h_l^* r_l / \sum_l |h_l|^2`.

    Each branch is weighted by its conjugate channel gain, so strong branches
    count for more and the output SNR is the sum of the branch SNRs.

    Parameters
    ----------
    received : array_like
        Received samples, shape ``(L, N)``: one row per branch.
    gains : array_like
        Complex channel gains :math:`h_l`, shape ``(L, N)`` or ``(L,)``.

    Returns
    -------
    ndarray of complex
        Combined estimate of the ``N`` transmitted symbols (unit gain on the
        signal).

    Notes
    -----
    Every column needs at least one nonzero gain. See
    :doc:`/guides/communications`.

    References
    ----------
    L. R. Kahn, "Ratio squarer," *Proc. IRE* 42 (1954) 1704.

    Examples
    --------
    >>> maximal_ratio_combine([[2.0], [1j]], [2.0, 1j]).tolist()
    [(1+0j)]
    """
    r = np.asarray(received)
    h = np.asarray(gains)
    if r.ndim != 2 or not r.size:
        raise ValueError("received must be a nonempty (branches, samples) array")
    h = h[:, None] if h.ndim == 1 else h
    if h.shape[0] != r.shape[0] or h.shape[1] not in (1, r.shape[1]):
        raise ValueError("gains must have one entry per branch (and per sample)")
    _finite(received=r, gains=h)
    power = np.sum(np.abs(h) ** 2, axis=0)
    if np.any(power == 0):
        raise ValueError("every sample needs at least one nonzero branch gain")
    return np.sum(np.conj(h) * r, axis=0) / power

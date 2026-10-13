r"""Orthogonal frequency-division multiplexing by the discrete Fourier transform.

Chang (1966) showed that overlapping subcarriers spaced by the symbol rate
remain orthogonal; Weinstein and Ebert (1971) generated them all with one
inverse DFT and added a cyclic prefix that turns a dispersive channel's
linear convolution into a circular one, so each subcarrier sees a single
complex gain.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer

__all__ = ["ofdm_modulate", "ofdm_demodulate"]


def ofdm_modulate(symbols, n_subcarriers, cyclic_prefix=0):
    r"""Map complex symbols onto OFDM time-domain blocks with an inverse DFT.

    Parameters
    ----------
    symbols : array_like
        Complex constellation symbols; the length must be a multiple of
        ``n_subcarriers``.
    n_subcarriers : int
        Number of subcarriers (IDFT size), at least 2.
    cyclic_prefix : int, optional
        Number of samples copied from the end of each block to its front,
        nonnegative and less than ``n_subcarriers``.

    Returns
    -------
    ndarray of complex
        Baseband samples, ``n_subcarriers + cyclic_prefix`` per block,
        normalised so average power equals average symbol energy.

    Notes
    -----
    See :doc:`/guides/communications`.

    References
    ----------
    S. B. Weinstein and P. M. Ebert, "Data transmission by frequency-division
    multiplexing using the discrete Fourier transform," *IEEE Trans. Commun.
    Technol.* 19 (1971) 628-634.

    Examples
    --------
    >>> x = ofdm_modulate([1, 1, 1, 1], 4, cyclic_prefix=1)
    >>> np.round(x.real, 12).tolist()
    [0.0, 2.0, 0.0, 0.0, 0.0]
    """
    _integer("n_subcarriers", n_subcarriers, 2)
    _integer("cyclic_prefix", cyclic_prefix, 0)
    if cyclic_prefix >= n_subcarriers:
        raise ValueError("cyclic_prefix must be shorter than a block")
    s = np.asarray(symbols, dtype=complex)
    _finite(symbols=s)
    if s.ndim != 1 or s.size % n_subcarriers:
        raise ValueError("symbols must be a vector whose length is a multiple of n_subcarriers")
    blocks = np.fft.ifft(s.reshape(-1, n_subcarriers), axis=1) * np.sqrt(n_subcarriers)
    if cyclic_prefix:
        blocks = np.hstack((blocks[:, -cyclic_prefix:], blocks))
    return blocks.reshape(-1)


def ofdm_demodulate(samples, n_subcarriers, cyclic_prefix=0):
    """Recover subcarrier symbols from OFDM samples (inverse of :func:`ofdm_modulate`).

    Parameters
    ----------
    samples : array_like
        Received baseband samples; the length must be a multiple of
        ``n_subcarriers + cyclic_prefix``.
    n_subcarriers : int
        Number of subcarriers, at least 2.
    cyclic_prefix : int, optional
        Cyclic-prefix length used by the transmitter.

    Returns
    -------
    ndarray of complex
        One complex value per subcarrier per block. After a channel shorter
        than the cyclic prefix, each equals the transmitted symbol times the
        channel's frequency response at that subcarrier.

    Notes
    -----
    See :doc:`/guides/communications`.

    Examples
    --------
    >>> rx = ofdm_demodulate(ofdm_modulate([1, -1, 1j, -1j], 4, 2), 4, 2)
    >>> np.allclose(rx, [1, -1, 1j, -1j])
    True
    """
    _integer("n_subcarriers", n_subcarriers, 2)
    _integer("cyclic_prefix", cyclic_prefix, 0)
    if cyclic_prefix >= n_subcarriers:
        raise ValueError("cyclic_prefix must be shorter than a block")
    x = np.asarray(samples, dtype=complex)
    _finite(samples=x)
    block = n_subcarriers + cyclic_prefix
    if x.ndim != 1 or x.size % block:
        raise ValueError("samples must be a vector of whole OFDM blocks")
    blocks = x.reshape(-1, block)[:, cyclic_prefix:]
    return (np.fft.fft(blocks, axis=1) / np.sqrt(n_subcarriers)).reshape(-1)

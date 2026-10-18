r"""Phasor helpers: polar and rectangular forms of complex amplitudes.

A sinusoid :math:`A\cos(\omega t + \phi)` is represented by the complex
number :math:`A e^{j\phi}`. Circuit textbooks quote phasors in polar form, so
these two helpers convert between the forms.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative

__all__ = ["phasor", "polar"]


def phasor(magnitude, phase):
    r"""Complex phasor :math:`A e^{j\phi}` from a magnitude and a phase.

    Parameters
    ----------
    magnitude : float or array_like
        Amplitude or RMS value [V, A or other], nonnegative.
    phase : float or array_like
        Phase angle [rad].

    Returns
    -------
    complex or ndarray of complex
        Phasor in rectangular form, same units as ``magnitude``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> complex(np.round(phasor(2.0, np.pi / 2), 12))
    2j
    """
    _nonnegative(magnitude=magnitude)
    _finite(phase=phase)
    return magnitude * np.exp(1j * np.asarray(phase, dtype=float))


def polar(z):
    r"""Magnitude and phase of a phasor, the inverse of :func:`phasor`.

    Parameters
    ----------
    z : complex or array_like
        Phasor [V, A, ohm or other].

    Returns
    -------
    magnitude : float or ndarray
        :math:`|z|`, same units as ``z``.
    phase : float or ndarray
        :math:`\arg z` in :math:`(-\pi, \pi]` [rad].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> polar(1j)
    (1.0, 1.5707963267948966)
    """
    _finite(z=z)
    z = np.asarray(z, dtype=complex)
    magnitude, phase = np.abs(z), np.angle(z)
    if z.ndim == 0:
        return float(magnitude), float(phase)
    return magnitude, phase

r"""Evans's root locus: closed-loop poles as the loop gain varies.

For unity negative feedback around :math:`K\,N(s)/D(s)`, the closed-loop
poles are the roots of :math:`D(s) + K N(s) = 0`. Evans (1948) gave rules
for sketching how they move from the open-loop poles (:math:`K = 0`) toward
the open-loop zeros or infinity.
"""

from __future__ import annotations

import numpy as np

__all__ = ["root_locus"]


def root_locus(num, den, gains):
    """Closed-loop characteristic roots of :math:`D(s) + K N(s)` for each gain.

    Parameters
    ----------
    num, den : array_like
        Open-loop numerator and denominator, descending powers of :math:`s`.
        The open loop must be strictly proper so the number of roots is
        constant.
    gains : array_like
        Nonempty vector of loop gains :math:`K` [-].

    Returns
    -------
    ndarray of complex
        Shape ``(len(gains), len(den) - 1)``; row ``i`` holds the roots for
        ``gains[i]``. Roots within a row are not tracked along branches.

    Notes
    -----
    See :doc:`/guides/control_systems`.

    References
    ----------
    W. R. Evans, "Graphical analysis of control systems," *Trans. AIEE* 67
    (1948) 547-551.

    Examples
    --------
    >>> np.sort(root_locus([1], [1, 3, 2], [0.25])[0].real).tolist()
    [-1.5, -1.5]
    """
    num, den = np.asarray(num, float), np.asarray(den, float)
    gains = np.asarray(gains, float)
    if (
        num.ndim != 1
        or den.ndim != 1
        or not num.size
        or len(den) < 2
        or len(num) >= len(den)
        or den[0] == 0
        or gains.ndim != 1
        or not gains.size
        or not all(np.all(np.isfinite(v)) for v in (num, den, gains))
    ):
        raise ValueError("requires finite strictly proper polynomials and a nonempty gain vector")
    return np.asarray([np.roots(np.polyadd(den, k * num)) for k in gains])

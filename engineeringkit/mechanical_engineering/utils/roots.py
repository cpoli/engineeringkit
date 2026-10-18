"""Root bracketing for the characteristic functions of mechanical models.

Holzer's residual torque and the Sommerfeld relation are smooth scalar
functions whose roots are wanted, not minimised. Scanning a grid for sign
changes and refining each bracket with Brent's method finds every root the
grid resolves.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

__all__ = ["bracketed_roots"]


def bracketed_roots(f, grid):
    """Roots of a scalar function between consecutive grid points where it changes sign.

    Parameters
    ----------
    f : callable
        ``f(x) -> float``, vectorized over a 1-D array.
    grid : array_like
        Increasing sample points; roots closer together than the grid
        spacing may be missed.

    Returns
    -------
    ndarray
        Roots in increasing order. A grid point where ``f`` is exactly zero
        is returned as a root.
    """
    x = np.asarray(grid, dtype=float)
    y = np.asarray(f(x), dtype=float)
    roots = list(x[y == 0])
    for i in np.flatnonzero(y[:-1] * y[1:] < 0):
        roots.append(brentq(lambda s: float(f(np.array([s]))[0]), x[i], x[i + 1], xtol=1e-14, rtol=1e-13))
    return np.sort(np.array(roots, dtype=float))

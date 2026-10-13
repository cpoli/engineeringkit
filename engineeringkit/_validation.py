"""Input checks shared by every engineeringkit domain.

Each helper raises :class:`ValueError` with the offending argument's name, so
a public function states its physical contract in one line at the top of its
body instead of relying on downstream NumPy/SciPy errors. All checks are
elementwise and therefore work for scalars and arrays alike.
"""

from __future__ import annotations

import numpy as np

__all__ = ["finite", "positive", "nonnegative", "interval", "integer", "vector"]


def finite(**values):
    """Reject empty, non-numeric, or nonfinite values.

    Parameters
    ----------
    **values
        ``name=value`` pairs; ``name`` appears in the error message.

    Raises
    ------
    ValueError
        If any value is empty, not numeric, a boolean, or contains NaN/inf.
    """
    for name, value in values.items():
        a = np.asarray(value)
        if not a.size or a.dtype == bool or not np.issubdtype(a.dtype, np.number) or not np.all(np.isfinite(a)):
            raise ValueError(f"{name} must contain finite numeric values")


def positive(**values):
    """Require finite, real, strictly positive values."""
    finite(**values)
    for name, value in values.items():
        if np.iscomplexobj(value) or np.any(np.asarray(value) <= 0):
            raise ValueError(f"{name} must be positive")


def nonnegative(**values):
    """Require finite, real, nonnegative values."""
    finite(**values)
    for name, value in values.items():
        if np.iscomplexobj(value) or np.any(np.asarray(value) < 0):
            raise ValueError(f"{name} must be nonnegative")


def interval(name, value, low, high, *, closed: bool | str = True):
    """Require every element of ``value`` to lie between ``low`` and ``high``.

    Parameters
    ----------
    name : str
        Argument name used in the error message.
    value : float or array_like
    low, high : float
        Bounds; either may be infinite.
    closed : bool or {"left", "right"}
        ``True`` for ``[low, high]``, ``False`` for ``(low, high)``, ``"left"``
        for ``[low, high)`` and ``"right"`` for ``(low, high]``.
    """
    finite(**{name: value})
    a = np.asarray(value)
    if np.iscomplexobj(a):
        raise ValueError(f"{name} must be real")
    left = closed is True or closed == "left"
    right = closed is True or closed == "right"
    too_low = a < low if left else a <= low
    too_high = a > high if right else a >= high
    if np.any(too_low | too_high):
        raise ValueError(f"{name} must lie in {'[' if left else '('}{low}, {high}{']' if right else ')'}")


def integer(name, value, minimum=0):
    """Require a Python or NumPy integer (not a bool) that is at least ``minimum``."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def vector(name, value, *, min_length=1):
    """Return ``value`` as a finite 1-D float array with at least ``min_length`` entries."""
    a = np.asarray(value)
    if a.ndim != 1 or len(a) < min_length:
        raise ValueError(f"{name} must be a 1-D array with at least {min_length} entries")
    finite(**{name: a})
    if np.iscomplexobj(a):
        raise ValueError(f"{name} must be real")
    return a.astype(float)

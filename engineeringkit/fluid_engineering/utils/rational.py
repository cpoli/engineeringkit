"""Exact rational linear algebra for dimensional analysis.

Dimension matrices hold small integers, and the exponents of the Pi groups
should come out as small integers too. Floating-point SVD null spaces mix the
groups and leave rounding noise, so the reduction here uses
:class:`fractions.Fraction` and is written out step by step.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, lcm

import numpy as np

__all__ = ["rational_rref", "integer_null_space"]


def rational_rref(matrix):
    """Reduced row echelon form of an integer matrix, computed exactly.

    Parameters
    ----------
    matrix : array_like of int
        Matrix of shape ``(m, n)``.

    Returns
    -------
    rref : list of list of Fraction
        The reduced matrix, ``m`` rows of ``n`` entries.
    pivots : list of int
        Column index of each pivot, in order; ``len(pivots)`` is the rank.
    """
    rows = [[Fraction(int(v)) for v in row] for row in np.asarray(matrix)]
    n = len(rows[0]) if rows else 0
    pivots: list[int] = []
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        lead = rows[r][c]
        rows[r] = [v / lead for v in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                factor = rows[i][c]
                rows[i] = [a - factor * b for a, b in zip(rows[i], rows[r])]
        pivots.append(c)
        r += 1
        if r == len(rows):
            break
    return rows, pivots


def integer_null_space(matrix):
    """Integer basis of the null space of an integer matrix.

    One basis vector is built per free (non-pivot) column: that column's
    entry is set to one and the pivot entries are read from the reduced
    matrix. Each vector is then scaled to coprime integers.

    Parameters
    ----------
    matrix : array_like of int
        Matrix of shape ``(m, n)``.

    Returns
    -------
    ndarray of int
        Basis vectors as rows, shape ``(n - rank, n)``; the free column of
        each row is positive.
    """
    rref, pivots = rational_rref(matrix)
    n = np.asarray(matrix).shape[1]
    basis = []
    for free in (c for c in range(n) if c not in pivots):
        v = [Fraction(0)] * n
        v[free] = Fraction(1)
        for row, p in zip(rref, pivots):
            v[p] = -row[free]
        scale = lcm(*(x.denominator for x in v))
        ints = [int(x * scale) for x in v]
        divisor = gcd(*ints)
        basis.append([x // divisor for x in ints])
    return np.array(basis, dtype=int).reshape(len(basis), n)

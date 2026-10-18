r"""Pipe networks: Hardy Cross's method of balancing heads.

In a looped water-distribution network the flows must satisfy continuity
at every junction and add up to zero head loss round every loop. Cross
(1936) started from any flows that satisfy continuity and corrected one
loop at a time by the Newton step

.. math::

    \Delta Q = -\frac{\sum s_i r_i Q_i |Q_i|^{n-1}}{n \sum r_i |Q_i|^{n-1}},

which leaves the junction balances untouched, because a loop correction
adds as much flow to each junction of the loop as it removes.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.fluid_engineering.core.base import NetworkResult

__all__ = ["hardy_cross"]


def hardy_cross(resistances, loops, initial_flows, n=2.0, tol=1e-10, max_iter=200):
    r"""Balance the flows of a looped pipe network by Hardy Cross's method.

    Parameters
    ----------
    resistances : array_like
        Resistance coefficient :math:`r_i` of each pipe, so that
        :math:`h_i = r_i Q_i |Q_i|^{n-1}` [s^n/m^(3n-1)], positive; for
        Darcy-Weisbach use :func:`pipe_resistance`.
    loops : array_like
        Loop matrix of shape ``(loops, pipes)``: +1 where the pipe's
        reference direction runs clockwise round the loop, -1 where it runs
        anticlockwise, 0 where the pipe is not in the loop.
    initial_flows : array_like
        Starting flow in each pipe [m^3/s], signed along its reference
        direction. It must satisfy continuity at every junction; the
        corrections preserve it.
    n : float, optional
        Head-loss exponent, at least 1: 2 for Darcy-Weisbach with fixed
        friction factor, 1.852 for Hazen-Williams.
    tol : float, optional
        Stop when the largest loop correction falls below this [m^3/s].
    max_iter : int, optional
        Largest number of sweeps over all loops.

    Returns
    -------
    NetworkResult
        Flows, head losses, loop residuals and the correction of each sweep.

    Raises
    ------
    ValueError
        If a loop has no pipes, or all its starting flows are zero so that
        the first correction is undefined.
    RuntimeError
        If the corrections do not fall below ``tol`` within ``max_iter`` sweeps.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Loops are corrected one after another
    (Gauss-Seidel order), as Cross did by hand; a pipe shared by two loops
    receives both corrections, with opposite signs.

    References
    ----------
    H. Cross, "Analysis of flow in networks of conduits or conductors,"
    *Univ. Illinois Eng. Exp. Stn. Bull.* 286 (1936).

    Examples
    --------
    One loop of two parallel pipes carrying 0.3 m³/s between two junctions:

    >>> result = hardy_cross([1.0, 4.0], [[1, -1]], [0.15, 0.15])
    >>> [round(float(q), 4) for q in result.flows]
    [0.2, 0.1]
    """
    r = _vector("resistances", resistances)
    _positive(resistances=r)
    Q = _vector("initial_flows", initial_flows).copy()
    A = np.atleast_2d(np.asarray(loops))
    _finite(loops=A)
    if A.shape[1:] != (r.size,) or Q.size != r.size:
        raise ValueError("loops must have one column per pipe, and initial_flows one value per pipe")
    if not np.all(np.isin(A, (-1, 0, 1))) or np.any(np.all(A == 0, axis=1)):
        raise ValueError("loops entries must be -1, 0 or 1, with at least one pipe per loop")
    _interval("n", n, 1, np.inf)
    _positive(tol=tol)
    _integer("max_iter", max_iter, 1)

    corrections = []
    for _ in range(max_iter):
        largest = 0.0
        for row in A:
            mask = row != 0
            scale = r[mask] * np.abs(Q[mask]) ** (n - 1)
            denominator = n * scale.sum()
            if denominator == 0:
                raise ValueError("every starting flow in a loop is zero; give the loop a nonzero guess")
            dQ = -np.sum(row[mask] * scale * Q[mask]) / denominator
            Q[mask] += row[mask] * dQ
            largest = max(largest, abs(dQ))
        corrections.append(largest)
        if largest < tol:
            h = r * Q * np.abs(Q) ** (n - 1)
            return NetworkResult(flows=Q, head_losses=h, loop_residuals=A @ h, corrections=np.array(corrections))
    raise RuntimeError(f"Hardy Cross did not converge in {max_iter} sweeps")

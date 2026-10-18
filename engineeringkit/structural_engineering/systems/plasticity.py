r"""Plastic collapse of beams by the mechanism (upper-bound) method.

When a section reaches its plastic moment :math:`M_p` it rotates freely as
a *plastic hinge*. A beam collapses when enough hinges form to turn it into
a mechanism; equating the work of the loads to the work dissipated in the
hinges, :math:`\sum P\delta = \sum M_p\theta`, gives the collapse load.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import positive as _positive

__all__ = ["plastic_collapse_load"]

_SUPPORTS = ("simple", "propped", "fixed")
_LOADS = ("point", "uniform")


def plastic_collapse_load(Mp, L, supports="fixed", load="point", position=None):
    r"""Collapse load of a prismatic beam with plastic moment :math:`M_p`.

    - ``"simple"``: :math:`P_c = M_pL/(ab)` for a point load at :math:`a`
      (:math:`b = L - a`), :math:`w_c = 8M_p/L^2` for a uniform load.
    - ``"propped"``: :math:`P_c = M_p(2/a + 1/b)`,
      :math:`w_c = 2(3 + 2\sqrt{2})M_p/L^2 \approx 11.66M_p/L^2`.
    - ``"fixed"``: :math:`P_c = 2M_pL/(ab)`, :math:`w_c = 16M_p/L^2`.

    A propped cantilever is fixed at the left end (:math:`x = 0`) and simply
    supported at the right.

    Parameters
    ----------
    Mp : float or array_like
        Plastic moment of the section [N m], positive (:math:`\sigma_y Z_p`).
    L : float or array_like
        Span [m], positive.
    supports : {"fixed", "propped", "simple"}, optional
        End conditions.
    load : {"point", "uniform"}, optional
        A single point load or a uniform load over the whole span.
    position : float, optional
        Distance of the point load from the left end [m], strictly inside
        the span; midspan by default.

    Returns
    -------
    float or ndarray
        Collapse load [N] for a point load, or collapse intensity [N/m] for a
        uniform load.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Rigid-plastic behaviour is assumed: no strain
    hardening, no instability and no reduction of :math:`M_p` by shear or
    axial force.

    References
    ----------
    G. Kazinczy, "Kísérletek befalazott tartókkal" (Experiments with clamped
    girders), *Betonszemle* 2 (1914), nos. 4-6; J. F. Baker, M. R. Horne and
    J. Heyman, *The Steel Skeleton*, vol. 2 (Cambridge, 1956).
    """
    _positive(Mp=Mp, L=L)
    if supports not in _SUPPORTS:
        raise ValueError(f"supports must be one of {_SUPPORTS}")
    if load not in _LOADS:
        raise ValueError(f"load must be one of {_LOADS}")
    Mp, L = np.asarray(Mp, float), np.asarray(L, float)
    if load == "uniform":
        if position is not None:
            raise ValueError("position applies to point loads only")
        factor = {"simple": 8.0, "propped": 2 * (3 + 2 * np.sqrt(2)), "fixed": 16.0}[supports]
        return factor * Mp / L**2
    a = L / 2 if position is None else np.asarray(position, float)
    _positive(position=a, remainder=L - a)
    b = L - a
    if supports == "simple":
        return Mp * L / (a * b)
    if supports == "propped":
        return Mp * (2 / a + 1 / b)
    return 2 * Mp * L / (a * b)

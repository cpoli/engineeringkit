r"""Cables and arches: the parabolic and catenary cable, and the three-hinged arch.

A cable carries load in tension only, so it takes the shape of its own
moment diagram: under a load uniform along the span it hangs in a parabola
with horizontal tension :math:`H = wL^2/(8f)`, and under its own weight in a
catenary. Inverted, the same shape is an arch in pure compression (Hooke's
"as hangs the flexible line, so but inverted will stand the rigid arch").
A three-hinged arch is statically determinate: the crown hinge gives the
extra equation for the thrust, :math:`H = M_\mathrm{beam}(L/2)/f`, and the
bending moment is the simple-beam moment less :math:`Hy`.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit.structural_engineering.core.base import ArchResult

__all__ = ["parabolic_cable", "catenary_cable", "three_hinged_arch"]


def parabolic_cable(w, span, sag):
    r"""Tension and length of a cable carrying a load uniform along the span.

    With supports at the same level and the sag :math:`f` at midspan, the
    cable is the parabola :math:`y = 4fx(L - x)/L^2` and

    .. math::

        H = \frac{wL^2}{8f},\qquad T_\max = \sqrt{H^2 + (wL/2)^2},\qquad
        s = \frac{L}{2}\sqrt{1 + 16n^2} + \frac{L}{8n}\sinh^{-1}4n,\quad n = f/L.

    Parameters
    ----------
    w : float or array_like
        Load per unit horizontal length [N/m], positive (a suspension
        bridge deck).
    span : float or array_like
        Horizontal distance between the supports [m], positive.
    sag : float or array_like
        Midspan sag [m], positive.

    Returns
    -------
    tuple of float or ndarray
        ``(H, T_max, length)``: horizontal tension [N], tension at the
        supports [N] and cable length [m].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. The cable is perfectly flexible and inextensible.
    """
    _positive(w=w, span=span, sag=sag)
    w, L, f = (np.asarray(v, float) for v in (w, span, sag))
    H = w * L**2 / (8 * f)
    n = f / L
    length = L / 2 * np.sqrt(1 + 16 * n**2) + L / (8 * n) * np.arcsinh(4 * n)
    return H, np.hypot(H, w * L / 2), length


def catenary_cable(w, span, sag):
    r"""Tension and length of a cable hanging under its own weight, the catenary.

    The cable is :math:`y = a(\cosh(x/a) - 1)` from its lowest point, with
    :math:`a = H/w`. The sag fixes :math:`a` through
    :math:`f = a(\cosh(L/2a) - 1)`, solved by Brent's method; then

    .. math::

        T_\max = H + wf,\qquad s = 2a\sinh\frac{L}{2a}.

    Parameters
    ----------
    w : float
        Weight per unit length of cable [N/m], positive.
    span : float
        Horizontal distance between the supports, at the same level [m], positive.
    sag : float
        Midspan sag [m], positive.

    Returns
    -------
    tuple of float
        ``(H, T_max, length)``: horizontal tension [N], tension at the
        supports [N] and cable length [m]. For a shallow cable they approach
        :func:`parabolic_cable` with the same ``w``.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Scalars only.

    References
    ----------
    The problem was posed by Jacob Bernoulli in *Acta Eruditorum* (1690)
    and solved there by Johann Bernoulli, G. W. Leibniz and C. Huygens
    (1691).
    """
    _positive(w=w, span=span, sag=sag)
    if np.ndim(w) or np.ndim(span) or np.ndim(sag):
        raise ValueError("catenary_cable takes scalars")
    ratio = sag / span
    # u = L/(2a) solves sinh(u/2)^2 / u = f/L; the left side exceeds u/4, so u < 4f/L
    upper = min(4 * ratio, 1400.0)
    u = brentq(lambda u: np.sinh(u / 2) ** 2 / u - ratio, 1e-9 * upper, upper, xtol=1e-15, rtol=1e-14)
    a = span / (2 * u)
    H = w * a
    return float(H), float(H + w * sag), float(2 * a * np.sinh(u))


def three_hinged_arch(span, rise, w=0.0, point_loads=(), n=401):
    r"""Bending moment and thrust of a parabolic three-hinged arch.

    The arch axis is :math:`y = 4fx(L - x)/L^2`, pinned at both springings
    and hinged at the crown. With :math:`M_0(x)` the moment of a simply
    supported beam of the same span under the same loads,

    .. math::

        H = \frac{M_0(L/2)}{f},\qquad M(x) = M_0(x) - Hy(x),

    and the axial force is :math:`N = -(H\cos\varphi + V_0\sin\varphi)`,
    with :math:`\varphi` the slope of the axis and :math:`V_0` the beam shear.
    A load uniform along the span gives :math:`M \equiv 0`: the parabola
    is its funicular shape.

    Parameters
    ----------
    span : float
        Distance between the springings, at the same level [m], positive.
    rise : float
        Height of the crown hinge above the springings [m], positive.
    w : float, optional
        Load per unit horizontal length over the whole span [N/m], downward
        positive.
    point_loads : sequence of (float, float), optional
        ``(position [m], P [N])``, downward positive.
    n : int, optional
        Number of stations, at least 2.

    Returns
    -------
    ArchResult

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Small deflections: the shape is not updated
    under load.

    References
    ----------
    R. Hooke, *A Description of Helioscopes, and Some Other Instruments*
    (London: John Martyn, 1676), postscript anagram; R. C. Hibbeler,
    *Structural Analysis*, 10th ed. (Pearson, 2018), section 5.5.
    """
    _positive(span=span, rise=rise)
    _finite(w=w)
    _integer("n", n, 2)
    L, f = float(span), float(rise)
    loads = []
    for load in point_loads:
        a, P = load
        _interval("load position", a, 0, L)
        _finite(P=P)
        loads.append((float(a), float(P)))
    left = w * L / 2 + sum(P * (L - a) / L for a, P in loads)
    right = w * L / 2 + sum(P * a / L for a, P in loads)

    def beam_moment(x):
        return left * x - w * x**2 / 2 - sum(P * np.where(x >= a, x - a, 0.0) for a, P in loads)

    H = float(beam_moment(L / 2)) / f
    x = np.linspace(0, L, n)
    y = 4 * f * x * (L - x) / L**2
    V0 = left - w * x - sum(P * (x >= a) for a, P in loads)
    phi = np.arctan(4 * f * (L - 2 * x) / L**2)
    moment = beam_moment(x) - H * y
    axial = -(H * np.cos(phi) + V0 * np.sin(phi))
    return ArchResult(x, y, moment, axial, (float(left), float(right)), H)

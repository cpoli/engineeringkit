r"""Beams: Navier's bending theory, Macaulay's singularity functions, Clapeyron's continuous beams and Winkler's influence lines.

Sign convention: loads act downward when positive, reactions upward,
applied couples clockwise, sagging bending moment is positive and downward
deflection is negative. Euler-Bernoulli beam theory, :math:`EIy'' = M`, is
assumed throughout. Unsymmetrical bending uses Hibbeler's and Gere's section
axes: :math:`y` up, :math:`z` horizontal, beam axis :math:`x`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.structural_engineering.core.base import BeamResult

__all__ = [
    "bending_stress",
    "unsymmetrical_bending_stress",
    "neutral_axis_angle",
    "SimplySupportedBeam",
    "Beam",
    "superpose_beam_results",
    "three_moment_equation",
    "influence_line",
]

_INFLUENCE_QUANTITIES = ("left_reaction", "right_reaction", "shear", "moment")
_SUPPORT_KINDS = ("pin", "roller", "fixed")


def bending_stress(M, y, I):
    r"""Flexure formula, :math:`\sigma = My/I`.

    Parameters
    ----------
    M : float or array_like
        Bending moment [N m], sagging positive.
    y : float or array_like
        Distance from the neutral axis [m]. With a sagging-positive ``M``,
        measure ``y`` downward so that tension comes out positive; pass the
        extreme-fibre distance ``c`` for the peak stress.
    I : float or array_like
        Second moment of area about the neutral axis [m^4], positive.

    Returns
    -------
    float or ndarray
        Normal stress [Pa].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    C.-L. Navier, *Résumé des leçons données à l'École des ponts et
    chaussées sur l'application de la mécanique à l'établissement des
    constructions et des machines* (Paris, 1826).
    """
    _finite(M=M, y=y, I=I)
    _positive(I=I)
    return np.asarray(M) * y / I


def _bending_coefficients(My, Mz, Iy, Iz, Iyz):
    """Return the stress gradients ``(d sigma/dy, d sigma/dz)`` [Pa/m] of the generalized flexure formula."""
    _finite(My=My, Mz=Mz, Iyz=Iyz)
    _positive(Iy=Iy, Iz=Iz)
    det = np.asarray(Iy) * Iz - np.asarray(Iyz) ** 2
    if np.any(det <= 0):
        raise ValueError("Iy * Iz - Iyz**2 must be positive")
    return -(np.asarray(Mz) * Iy + np.asarray(My) * Iyz) / det, (np.asarray(My) * Iz + np.asarray(Mz) * Iyz) / det


def unsymmetrical_bending_stress(My, Mz, y, z, Iy, Iz, Iyz=0.0):
    r"""Normal stress under bending about two axes, the generalized flexure formula.

    .. math::

        \sigma = -\frac{M_zI_y + M_yI_{yz}}{I_yI_z - I_{yz}^2}\,y
                 + \frac{M_yI_z + M_zI_{yz}}{I_yI_z - I_{yz}^2}\,z,

    which reduces to :math:`\sigma = -M_zy/I_z + M_yz/I_y` on principal
    axes. Positive :math:`M_z` compresses the fibres at positive :math:`y`
    (for a beam along :math:`x` loaded downward, it is the sagging moment).

    Parameters
    ----------
    My, Mz : float or array_like
        Moment components about the centroidal :math:`y` (vertical) and
        :math:`z` (horizontal) axes [N m], right-hand rule.
    y, z : float or array_like
        Coordinates of the point from the centroid [m], :math:`y` up.
    Iy, Iz : float or array_like
        Centroidal second moments of area [m^4], positive.
    Iyz : float or array_like, optional
        Product of inertia :math:`\int yz\,dA` [m^4]; zero on principal axes.

    Returns
    -------
    float or ndarray
        Normal stress [Pa], tension positive.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Loads must pass through the shear centre, or the
    beam also twists.

    References
    ----------
    R. C. Hibbeler, *Mechanics of Materials*, 10th ed. (Pearson, 2017),
    section 6.5; J. M. Gere and B. J. Goodno, *Mechanics of Materials*,
    9th ed. (Cengage, 2018), sections 6.4-6.5.
    """
    _finite(y=y, z=z)
    dy, dz = _bending_coefficients(My, Mz, Iy, Iz, Iyz)
    return dy * np.asarray(y) + dz * np.asarray(z)


def neutral_axis_angle(My, Mz, Iy, Iz, Iyz=0.0):
    r"""Inclination of the neutral axis under bending about two axes.

    The neutral axis is the line of zero :func:`unsymmetrical_bending_stress`
    through the centroid. On principal axes, with the moment vector at
    :math:`\theta` from the :math:`z` axis, :math:`\tan\alpha =
    (I_z/I_y)\tan\theta`: unless :math:`I_y = I_z` the neutral axis is not
    perpendicular to the plane of the loads.

    Parameters
    ----------
    My, Mz : float or array_like
        Moment components [N m], not both zero.
    Iy, Iz : float or array_like
        Centroidal second moments of area [m^4], positive.
    Iyz : float or array_like, optional
        Product of inertia [m^4].

    Returns
    -------
    float or ndarray
        Angle :math:`\alpha` of the neutral axis from the :math:`z` axis [rad],
        counterclockwise positive, in :math:`(-\pi/2, \pi/2]`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    dy, dz = _bending_coefficients(My, Mz, Iy, Iz, Iyz)
    if np.any((dy == 0) & (dz == 0)):
        raise ValueError("My and Mz must not both be zero")
    alpha = np.arctan2(-dz, dy)
    return np.where(alpha > np.pi / 2, alpha - np.pi, np.where(alpha <= -np.pi / 2, alpha + np.pi, alpha))


@dataclass(frozen=True)
class SimplySupportedBeam:
    """Prismatic Euler-Bernoulli beam on a pin and a roller.

    Parameters
    ----------
    length : float
        Span [m], positive.
    young_modulus : float
        Young's modulus [Pa], positive.
    second_moment : float
        Second moment of area [m^4], positive.
    """

    length: float
    young_modulus: float
    second_moment: float

    def __post_init__(self):
        _positive(length=self.length, young_modulus=self.young_modulus, second_moment=self.second_moment)

    def solve_point_load(self, load: float, position: float | None = None, n: int = 401):
        """Shear, moment and exact deflection under one point load.

        Parameters
        ----------
        load : float
            Point load [N], downward positive.
        position : float, optional
            Distance of the load from the left support [m]; midspan by default.
        n : int, optional
            Number of stations, at least 2.

        Returns
        -------
        BeamResult
        """
        _finite(load=load)
        _integer("n", n, 2)
        L = self.length
        a = L / 2 if position is None else position
        if not 0 <= a <= L:
            raise ValueError("position must lie on beam")
        b = L - a
        ra = load * b / L
        rb = load * a / L
        x = np.linspace(0, L, n)
        V = np.where(x < a, ra, ra - load)
        M = np.where(x <= a, ra * x, rb * (L - x))
        # exact Euler-Bernoulli deflection, downward negative
        E, I = self.young_modulus, self.second_moment
        y = np.empty_like(x)
        left = x <= a
        y[left] = -load * b * x[left] * (L**2 - b**2 - x[left] ** 2) / (6 * L * E * I)
        xr = L - x[~left]
        y[~left] = -load * a * xr * (L**2 - a**2 - xr**2) / (6 * L * E * I)
        return BeamResult(x, V, M, y, (ra, rb))

    def solve_uniform_load(self, w: float, n: int = 401):
        """Shear, moment and exact deflection under a uniform load.

        Parameters
        ----------
        w : float
            Load intensity [N/m], downward positive.
        n : int, optional
            Number of stations, at least 2.

        Returns
        -------
        BeamResult
        """
        _finite(w=w)
        _integer("n", n, 2)
        L = self.length
        x = np.linspace(0, L, n)
        R = w * L / 2
        V = R - w * x
        M = R * x - w * x * x / 2
        y = -w * x * (L**3 - 2 * L * x * x + x**3) / (24 * self.young_modulus * self.second_moment)
        return BeamResult(x, V, M, y, (R, R))


def _bracket(x, a, p):
    """Macaulay bracket :math:`\\langle x - a\\rangle^p`, zero for ``x < a`` (``p >= 0``)."""
    return np.where(x >= a, (x - a) ** p, 0.0)


def _integrate_terms(terms, x):
    """Return ``(M, V, EI slope, EI deflection)`` at ``x`` for ``coef <x - a>^p`` moment terms, without constants."""
    x = np.asarray(x, float)
    M, V, slope, deflection = np.zeros((4,) + x.shape)
    for coef, a, p in terms:
        M += coef * _bracket(x, a, p)
        if p:
            V += coef * p * _bracket(x, a, p - 1)
        slope += coef * _bracket(x, a, p + 1) / (p + 1)
        deflection += coef * _bracket(x, a, p + 2) / ((p + 1) * (p + 2))
    return M, V, slope, deflection


@dataclass(frozen=True)
class Beam:
    r"""Prismatic Euler-Bernoulli beam on any supports, solved by Macaulay's method.

    The bending moment is written once for the whole beam with
    singularity functions :math:`\langle x - a\rangle^n`, which are zero
    to the left of :math:`a`:

    .. math::

        M(x) = \sum R_k\langle x - a_k\rangle - \sum P\langle x - a\rangle
               + \sum C\langle x - a\rangle^0 - \tfrac{w}{2}\langle x - a\rangle^2 + \cdots

    and integrated twice, :math:`EIy = \iint M\,dx\,dx + C_1x + C_2`. The
    unknown reactions and the two constants follow from vertical and
    moment equilibrium and the support conditions (:math:`y = 0` at every
    support, :math:`y' = 0` at a fixed one), so statically indeterminate
    beams are solved the same way as determinate ones.

    Parameters
    ----------
    length : float
        Length [m], positive; the beam runs from :math:`x = 0` to ``length``.
    young_modulus : float
        Young's modulus [Pa], positive.
    second_moment : float
        Second moment of area [m^4], positive.
    supports : sequence of (float, str)
        ``(position [m], kind)`` for each support, with ``kind`` one of
        ``"pin"``, ``"roller"`` (both restrain deflection only) or
        ``"fixed"`` (restrains deflection and slope).

    References
    ----------
    W. H. Macaulay, "Note on the deflection of beams," *Messenger of
    Mathematics* 48 (1919) 129-130; A. Clebsch, *Theorie der Elasticität
    fester Körper* (Leipzig: Teubner, 1862).
    """

    length: float
    young_modulus: float
    second_moment: float
    supports: tuple

    def __post_init__(self):
        _positive(length=self.length, young_modulus=self.young_modulus, second_moment=self.second_moment)
        if not len(self.supports):
            raise ValueError("at least one support is required")
        for support in self.supports:
            if len(support) != 2 or support[1] not in _SUPPORT_KINDS:
                raise ValueError(f"each support must be (position, kind) with kind one of {_SUPPORT_KINDS}")
            _interval("support position", support[0], 0, self.length)

    def solve(self, point_loads=(), moments=(), distributed_loads=(), n: int = 401):
        """Reactions, shear, moment, slope and deflection under any combination of loads.

        Parameters
        ----------
        point_loads : sequence of (float, float), optional
            ``(position [m], P [N])``, downward positive.
        moments : sequence of (float, float), optional
            ``(position [m], C [N m])`` applied couples, clockwise positive.
        distributed_loads : sequence of tuple, optional
            ``(start [m], end [m], w_start [N/m])`` for a uniform load or
            ``(start, end, w_start, w_end)`` for a linearly varying one,
            downward positive.
        n : int, optional
            Number of equally spaced stations, at least 2.

        Returns
        -------
        BeamResult
            With one reaction per support, the reaction moments of the
            fixed supports and the slope. Shear at a point load or support
            takes the value just to its right, except at ``x = length``,
            where it takes the value just to its left.

        Raises
        ------
        ValueError
            If the supports leave the beam free to move as a mechanism.
        """
        _integer("n", n, 2)
        L = self.length
        terms = []
        for load in point_loads:
            a, P = load
            _interval("load position", a, 0, L)
            _finite(P=P)
            terms.append((-float(P), float(a), 1))
        for couple in moments:
            a, C = couple
            _interval("moment position", a, 0, L)
            _finite(C=C)
            terms.append((float(C), float(a), 0))
        for load in distributed_loads:
            if len(load) not in (3, 4):
                raise ValueError("a distributed load must be (start, end, w_start) or (start, end, w_start, w_end)")
            a, b, w1 = load[:3]
            w2 = load[3] if len(load) == 4 else w1
            _interval("load start", a, 0, L)
            _interval("load end", b, a, L, closed="right")
            _finite(w_start=w1, w_end=w2)
            k = (w2 - w1) / (b - a)
            terms += [(-w1 / 2, a, 2), (-k / 6, a, 3), (w2 / 2, b, 2), (k / 6, b, 3)]
        unknowns = [(float(s[0]), 1) for s in self.supports] + [(float(s[0]), 0) for s in self.supports if s[1] == "fixed"]
        size = len(unknowns) + 2
        A = np.zeros((size, size))
        rhs = np.zeros(size)
        rows = [("shear", L), ("moment", L)] + [("deflection", float(s[0])) for s in self.supports]
        rows += [("slope", float(s[0])) for s in self.supports if s[1] == "fixed"]
        for r, (quantity, at) in enumerate(rows):
            index = ("moment", "shear", "slope", "deflection").index(quantity)
            for c, (a, p) in enumerate(unknowns):
                A[r, c] = float(_integrate_terms([(1.0, a, p)], at)[index])
            A[r, -2:] = {"shear": (0, 0), "moment": (0, 0), "slope": (1, 0), "deflection": (at, 1)}[quantity]
            rhs[r] = -float(_integrate_terms(terms, at)[index])
        row_scale = np.max(np.abs(A), axis=1, keepdims=True)
        col_scale = np.max(np.abs(A), axis=0, keepdims=True)
        scaled = A / row_scale / col_scale
        if np.linalg.matrix_rank(scaled) < size:
            raise ValueError("the supports leave the beam free to move as a mechanism")
        solution = np.linalg.solve(scaled, rhs / row_scale[:, 0]) / col_scale[0]
        values, (C1, C2) = solution[:-2], solution[-2:]
        all_terms = terms + [(value, a, p) for value, (a, p) in zip(values, unknowns)]
        x = np.linspace(0, L, n)
        M, V, slope, deflection = _integrate_terms(all_terms, x)
        # as in SimplySupportedBeam, the end station shows the shear just left of a reaction or load at x = L
        V[-1] = _integrate_terms([term for term in all_terms if term[1] < L], L)[1]
        EI = self.young_modulus * self.second_moment
        nsup = len(self.supports)
        return BeamResult(
            x,
            V,
            M,
            (deflection + C1 * x + C2) / EI,
            tuple(float(v) for v in values[:nsup]),
            tuple(float(v) for v in values[nsup:]),
            (slope + C1) / EI,
        )


def superpose_beam_results(*results):
    """Add the responses of one beam to several load cases (the principle of superposition).

    Shear, moment, deflection, slope and reactions of a linear elastic beam
    are linear in the loads, so the response to a combination of loads is
    the sum of the responses to each.

    Parameters
    ----------
    *results : BeamResult
        Responses of the same beam, computed at the same stations and with
        the same supports.

    Returns
    -------
    BeamResult
        The summed response. The slope is summed when every result has one,
        and is ``None`` otherwise.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    if not results:
        raise ValueError("at least one result is required")
    first = results[0]
    for other in results[1:]:
        if not np.array_equal(other.x, first.x):
            raise ValueError("all results must use the same stations")
        if len(other.reactions) != len(first.reactions) or len(other.support_moments) != len(first.support_moments):
            raise ValueError("all results must come from the same supports")
    slopes = [r.slope for r in results]
    return BeamResult(
        first.x,
        sum(r.shear for r in results),
        sum(r.moment for r in results),
        sum(r.deflection for r in results),
        tuple(float(v) for v in np.sum([r.reactions for r in results], axis=0)),
        tuple(float(v) for v in np.sum([r.support_moments for r in results], axis=0).reshape(-1)),
        None if any(s is None for s in slopes) else sum(slopes),
    )


def three_moment_equation(spans, w=0.0, I=1.0):
    r"""Support moments of a continuous beam by Clapeyron's three-moment equation.

    For the supports :math:`i-1, i, i+1` around spans :math:`L_i` and
    :math:`L_{i+1}` carrying uniform loads :math:`w_i, w_{i+1}`,

    .. math::

        M_{i-1}\frac{L_i}{I_i} + 2M_i\left(\frac{L_i}{I_i} + \frac{L_{i+1}}{I_{i+1}}\right)
        + M_{i+1}\frac{L_{i+1}}{I_{i+1}} = -\frac{w_iL_i^3}{4I_i} - \frac{w_{i+1}L_{i+1}^3}{4I_{i+1}},

    one equation per interior support, with :math:`M = 0` at the two
    simply supported ends. The modulus is taken to be the same in every span.

    Parameters
    ----------
    spans : array_like
        Span lengths [m], positive, left to right.
    w : float or array_like, optional
        Uniform load on each span [N/m], downward positive.
    I : float or array_like, optional
        Second moment of area of each span [m^4], positive. Only ratios
        between spans matter.

    Returns
    -------
    ndarray
        Bending moment at each of the ``len(spans) + 1`` supports [N m],
        sagging positive (interior supports of a loaded beam hog, so they
        come out negative).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    É. Clapeyron, "Calcul d'une poutre élastique reposant librement sur des
    appuis inégalement espacés," *C. R. Acad. Sci. Paris* 45 (1857)
    1076-1080.

    Examples
    --------
    Two equal spans: the middle support moment is :math:`-wL^2/8`.

    >>> three_moment_equation([4.0, 4.0], w=1000.0)
    array([    0., -2000.,     0.])
    """
    L = _vector("spans", spans)
    _positive(spans=L)
    _finite(w=w, I=I)
    _positive(I=I)
    n = L.size
    wk = np.broadcast_to(np.asarray(w, float), (n,))
    Ik = np.broadcast_to(np.asarray(I, float), (n,))
    moments = np.zeros(n + 1)
    if n == 1:
        return moments
    a = L / Ik
    load = wk * L**3 / (4 * Ik)
    A = np.zeros((n - 1, n - 1))
    rhs = np.zeros(n - 1)
    for row in range(n - 1):  # interior support row + 1, between spans row and row + 1
        A[row, row] = 2 * (a[row] + a[row + 1])
        if row > 0:
            A[row, row - 1] = a[row]
        if row < n - 2:
            A[row, row + 1] = a[row + 1]
        rhs[row] = -(load[row] + load[row + 1])
    moments[1:-1] = np.linalg.solve(A, rhs)
    return moments


def influence_line(length, positions, quantity="moment", at=None):
    r"""Influence line of a simply supported beam for a moving unit load.

    Gives the value of a response quantity when a unit downward load
    stands at each of ``positions``. Multiplying by a load and summing over
    several loads gives the response to any load train, and the position of
    the train that maximises it.

    Parameters
    ----------
    length : float
        Span [m], positive.
    positions : float or array_like
        Positions of the unit load from the left support [m], in ``[0, length]``.
    quantity : {"moment", "shear", "left_reaction", "right_reaction"}, optional
        Response quantity.
    at : float, optional
        Section [m] at which ``"moment"`` or ``"shear"`` is measured,
        in ``[0, length]``. Not used for reactions.

    Returns
    -------
    float or ndarray
        Response per unit load: reactions [-] (N/N), shear [-] (N/N) or
        moment [m] (N m/N). Shear follows :class:`SimplySupportedBeam`: with
        the load at or left of the section, the shear there is
        :math:`-\xi/L`; otherwise :math:`1 - \xi/L`.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    E. Winkler, *Die Lehre von der Elasticität und Festigkeit* (Prague:
    Dominicus, 1867); O. Mohr, "Beitrag zur Theorie der Holz- und
    Eisen-Constructionen," *Z. Architekten- und Ingenieur-Vereins zu
    Hannover* 14 (1868) 19-51.
    """
    _positive(length=length)
    _interval("positions", positions, 0, length)
    xi = np.asarray(positions, float)
    if quantity not in _INFLUENCE_QUANTITIES:
        raise ValueError(f"quantity must be one of {_INFLUENCE_QUANTITIES}")
    if quantity == "left_reaction":
        return 1 - xi / length
    if quantity == "right_reaction":
        return xi / length
    if at is None:
        raise ValueError(f"at is required for the {quantity} influence line")
    _interval("at", at, 0, length)
    if quantity == "shear":
        return np.where(xi <= at, -xi / length, 1 - xi / length)
    return np.where(xi <= at, xi * (length - at) / length, at * (length - xi) / length)

r"""Beams: Navier's bending theory, Clapeyron's continuous beams and Winkler's influence lines.

Sign convention: loads act downward when positive, reactions upward,
sagging bending moment is positive and downward deflection is negative.
Euler-Bernoulli beam theory, :math:`EIy'' = M`, is assumed throughout.
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

__all__ = ["bending_stress", "SimplySupportedBeam", "three_moment_equation", "influence_line"]

_INFLUENCE_QUANTITIES = ("left_reaction", "right_reaction", "shear", "moment")


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

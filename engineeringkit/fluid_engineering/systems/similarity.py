r"""Similarity: dimensionless groups and Buckingham's Pi theorem.

Two flows behave alike when their dimensionless groups match. Froude (1868)
scaled ship models by :math:`v/\sqrt{gL}`; Buckingham (1914) proved that a
relation among :math:`n` variables built from :math:`r` independent
dimensions can always be rewritten among :math:`n - r` dimensionless
groups, and the groups are the null space of the dimension matrix.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive
from engineeringkit.fluid_engineering.core.base import PiGroups
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY as G
from engineeringkit.fluid_engineering.utils.rational import integer_null_space, rational_rref

__all__ = ["froude_number", "weber_number", "euler_number", "mach_number", "buckingham_pi"]

_BASE_ORDER = ("M", "L", "T", "Θ", "I", "N", "J")


def froude_number(v, length, g=G):
    r"""Froude number :math:`Fr = v/\sqrt{gL}`.

    Parameters
    ----------
    v : float or array_like
        Speed [m/s].
    length : float or array_like
        Characteristic length [m], positive (waterline length of a hull,
        depth of a channel).
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Froude number [-], the ratio of inertia to gravity.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    References
    ----------
    W. Froude, "Observations and suggestions on the subject of determining
    by experiment the resistance of ships," memorandum to the Admiralty
    (1868), in *The Papers of William Froude* (London: INA, 1955).

    Examples
    --------
    >>> float(froude_number(3.0, 0.9, g=10.0))
    1.0
    """
    _finite(v=v, length=length, g=g)
    _positive(length=length, g=g)
    return v / np.sqrt(g * length)


def weber_number(rho, v, length, surface_tension):
    r"""Weber number :math:`We = \rho v^2 L/\sigma`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    v : float or array_like
        Speed [m/s].
    length : float or array_like
        Characteristic length, such as a droplet diameter [m], positive.
    surface_tension : float or array_like
        Surface tension [N/m], positive.

    Returns
    -------
    float or ndarray
        Weber number [-], the ratio of inertia to surface tension.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> weber_number(1000.0, 1.0, 0.01, 0.1)
    100.0
    """
    _finite(rho=rho, v=v, length=length, surface_tension=surface_tension)
    _positive(rho=rho, length=length, surface_tension=surface_tension)
    return rho * v * v * length / surface_tension


def euler_number(delta_p, rho, v):
    r"""Euler number :math:`Eu = \Delta p/(\rho v^2)`.

    Parameters
    ----------
    delta_p : float or array_like
        Pressure difference [Pa].
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    v : float or array_like
        Speed [m/s], nonzero.

    Returns
    -------
    float or ndarray
        Euler number [-], the ratio of pressure forces to inertia.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Some texts use
    :math:`\Delta p/(\tfrac12\rho v^2)`, twice this value.

    Examples
    --------
    >>> euler_number(4000.0, 1000.0, 2.0)
    1.0
    """
    _finite(delta_p=delta_p, rho=rho, v=v)
    _positive(rho=rho)
    if np.any(np.asarray(v) == 0):
        raise ValueError("velocity must be nonzero")
    return delta_p / (rho * v * v)


def mach_number(v, speed_of_sound):
    r"""Mach number :math:`Ma = v/a`.

    Parameters
    ----------
    v : float or array_like
        Speed [m/s].
    speed_of_sound : float or array_like
        Local speed of sound [m/s], positive.

    Returns
    -------
    float or ndarray
        Mach number [-]; compressibility matters above about 0.3.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> mach_number(170.0, 340.0)
    0.5
    """
    _finite(v=v, speed_of_sound=speed_of_sound)
    _positive(speed_of_sound=speed_of_sound)
    return v / speed_of_sound


def buckingham_pi(variables, repeating=None):
    r"""Dimensionless groups of a set of variables by Buckingham's Pi theorem.

    The variables' dimensions form a matrix :math:`D` (base dimensions by
    variables). Every vector :math:`x` with :math:`Dx = 0` lists the
    exponents of one dimensionless product, and a basis of that null space
    has :math:`n - \operatorname{rank} D` vectors. The repeating variables
    are placed first, so the exact row reduction pivots on them and each
    group contains exactly one other variable, raised to a positive power.

    Parameters
    ----------
    variables : mapping
        Variable name to its dimensions, itself a mapping from a base
        dimension (``"M"``, ``"L"``, ``"T"``, ``"Θ"``, ...) to an integer
        exponent; for example ``{"rho": {"M": 1, "L": -3}}``. Dimensionless
        variables map to ``{}``.
    repeating : sequence of str, optional
        Names of the repeating variables, as many as the rank of the
        dimension matrix and dimensionally independent. Defaults to the
        first independent variables in ``variables`` order.

    Returns
    -------
    PiGroups
        Variables, dimension matrix and the integer exponents of each group.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The exponents come from exact rational
    arithmetic (:func:`~engineeringkit.fluid_engineering.utils.integer_null_space`).
    The theorem fixes how many groups there are, not which ones: another
    choice of repeating variables gives products of these.

    References
    ----------
    E. Buckingham, "On physically similar systems; illustrations of the use
    of dimensional equations," *Phys. Rev.* 4 (1914) 345-376.

    Examples
    --------
    >>> pi = buckingham_pi(
    ...     {
    ...         "rho": {"M": 1, "L": -3},
    ...         "V": {"L": 1, "T": -1},
    ...         "D": {"L": 1},
    ...         "mu": {"M": 1, "L": -1, "T": -1},
    ...     }
    ... )
    >>> pi.groups()
    ['rho^-1 V^-1 D^-1 mu^1']
    """
    if not isinstance(variables, Mapping) or len(variables) < 2:
        raise ValueError("variables must map at least two names to their dimensions")
    for name, dims in variables.items():
        if not isinstance(dims, Mapping) or any(isinstance(e, bool) or not isinstance(e, (int, np.integer)) for e in dims.values()):
            raise ValueError(f"dimensions of {name!r} must map base dimensions to integer exponents")
    used = {b for dims in variables.values() for b, e in dims.items() if e}
    bases = tuple(b for b in _BASE_ORDER if b in used) + tuple(sorted(used - set(_BASE_ORDER)))
    names = list(variables)

    def matrix(order):
        return np.array([[int(variables[n].get(b, 0)) for n in order] for b in bases], dtype=int).reshape(len(bases), len(order))

    _, pivots = rational_rref(matrix(names)) if bases else ([], [])
    rank = len(pivots)
    if repeating is None:
        repeating = [names[p] for p in pivots]
    repeating = list(repeating)
    if len(set(repeating)) != len(repeating) or not set(repeating) <= set(names):
        raise ValueError("repeating variables must be distinct names from variables")
    if len(repeating) != rank or (rank and len(rational_rref(matrix(repeating))[1]) != rank):
        raise ValueError(f"need {rank} dimensionally independent repeating variables")
    order = repeating + [n for n in names if n not in repeating]
    D = matrix(order)
    exponents = integer_null_space(D) if bases else np.eye(len(order), dtype=int)
    return PiGroups(
        names=tuple(order),
        bases=bases,
        dimension_matrix=D,
        exponents=exponents,
        repeating=tuple(repeating),
    )

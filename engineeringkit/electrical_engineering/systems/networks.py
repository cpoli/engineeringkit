r"""Equivalent two-terminal networks and maximum power transfer.

Seen from two terminals, any linear network of sources and resistors acts
as one source behind one resistance: a voltage source in series (Thévenin)
or a current source in parallel (Norton). The load that draws the most
power from it matches that resistance (Jacobi).
"""

from __future__ import annotations

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["thevenin_from_source", "norton_from_thevenin", "maximum_power_transfer"]


def thevenin_from_source(vs, rs):
    r"""Thévenin equivalent of a real voltage source with internal resistance.

    A source :math:`V_s` behind :math:`R_s` already is a Thévenin
    equivalent; the function validates and returns it as a pair, so it can
    be passed to :func:`norton_from_thevenin` or
    :func:`maximum_power_transfer`. To reduce a whole network, use
    :func:`~engineeringkit.electrical_engineering.nodal_analysis` for the
    open-circuit voltage and the short-circuit current.

    Parameters
    ----------
    vs : float
        Open-circuit source voltage [V].
    rs : float
        Internal resistance [ohm], nonnegative.

    Returns
    -------
    tuple of float
        ``(V_th [V], R_th [ohm])``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    L. Thévenin, "Extension de la loi d'Ohm aux circuits électromoteurs
    complexes," *Annales Télégraphiques* 10 (1883) 222-224.

    Examples
    --------
    >>> thevenin_from_source(9, 2)
    (9.0, 2.0)
    """
    _finite(vs=vs, rs=rs)
    _nonnegative(rs=rs)
    return float(vs), float(rs)


def norton_from_thevenin(vth, rth):
    r"""Norton equivalent :math:`I_N = V_{th}/R_{th}` in parallel with :math:`R_{th}`.

    Parameters
    ----------
    vth : float
        Thévenin (open-circuit) voltage [V].
    rth : float
        Thévenin resistance [ohm], positive (an ideal voltage source has no
        Norton form).

    Returns
    -------
    tuple of float
        ``(I_N [A], R_N [ohm])``; :math:`I_N` is the short-circuit current.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    E. L. Norton, "Design of finite networks for uniform frequency
    characteristic," Bell Laboratories Technical Report TM26-0-1860 (1926).

    Examples
    --------
    >>> norton_from_thevenin(10, 5)
    (2.0, 5)
    """
    _finite(vth=vth, rth=rth)
    _positive(rth=rth)
    return vth / rth, rth


def maximum_power_transfer(vth, rth):
    r"""Matched load :math:`R_L = R_{th}` and the power :math:`V_{th}^2/4R_{th}` it receives.

    Parameters
    ----------
    vth : float
        Thévenin (open-circuit) voltage [V].
    rth : float
        Thévenin resistance [ohm], positive.

    Returns
    -------
    tuple of float
        ``(R_L [ohm], P_max [W])``. At the match the efficiency is 50%.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> maximum_power_transfer(12.0, 6.0)
    (6.0, 6.0)
    """
    _finite(vth=vth, rth=rth)
    _positive(rth=rth)
    return rth, vth * vth / (4 * rth)

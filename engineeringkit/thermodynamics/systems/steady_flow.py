r"""Steady-flow devices: the first law for an open system and isentropic efficiencies.

Per unit mass flowing through a control volume in steady state,

.. math:: q - w = (h_{out} - h_{in}) + \tfrac12(V_{out}^2 - V_{in}^2) + g(z_{out} - z_{in}),

so a nozzle turns enthalpy into kinetic energy, a turbine into shaft work, a
compressor or pump does the reverse, a throttle keeps its enthalpy and a
mixing chamber averages the enthalpies of its inlets. Real turbines,
compressors and nozzles are compared with the isentropic device between the
same pressures.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.constants import STANDARD_GRAVITY
from engineeringkit.thermodynamics.core.base import SteamState
from engineeringkit.thermodynamics.systems.steam import steam_properties

__all__ = [
    "steady_flow_work",
    "nozzle_exit_velocity",
    "turbine_isentropic_efficiency",
    "compressor_isentropic_efficiency",
    "nozzle_isentropic_efficiency",
    "expansion_outlet_enthalpy",
    "compression_outlet_enthalpy",
    "mixing_chamber",
    "throttle_steam",
]


def steady_flow_work(h_in, h_out, *, q=0.0, V_in=0.0, V_out=0.0, z_in=0.0, z_out=0.0):
    r"""Shaft work per unit mass of a steady-flow device from its energy balance.

    .. math:: w = q + (h_{in} - h_{out}) + \tfrac12(V_{in}^2 - V_{out}^2) + g(z_{in} - z_{out})

    Parameters
    ----------
    h_in, h_out : float or array_like
        Specific enthalpy at the inlet and the outlet [J/kg].
    q : float or array_like, optional
        Heat added per unit mass [J/kg]; negative for heat lost.
    V_in, V_out : float or array_like, optional
        Mean flow velocity at the inlet and the outlet [m/s], nonnegative.
    z_in, z_out : float or array_like, optional
        Elevation of the inlet and the outlet [m].

    Returns
    -------
    float or ndarray
        Work done *by* the fluid per unit mass [J/kg]: positive for a
        turbine, negative for a compressor or pump.

    Notes
    -----
    One inlet, one outlet, steady state. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(steady_flow_work(3000e3, 2500e3)), 6)
    500000.0
    """
    _finite(h_in=h_in, h_out=h_out, q=q, z_in=z_in, z_out=z_out)
    _nonnegative(V_in=V_in, V_out=V_out)
    kinetic = 0.5 * (np.asarray(V_in, dtype=float) ** 2 - np.asarray(V_out, dtype=float) ** 2)
    return (q + (np.asarray(h_in, dtype=float) - h_out) + kinetic + STANDARD_GRAVITY * (np.asarray(z_in, dtype=float) - z_out))[()]


def nozzle_exit_velocity(h_in, h_out, V_in=0.0):
    r"""Exit velocity of an adiabatic nozzle, :math:`V_2 = \sqrt{V_1^2 + 2(h_1 - h_2)}`.

    Parameters
    ----------
    h_in, h_out : float or array_like
        Specific enthalpy at the inlet and the exit [J/kg].
    V_in : float or array_like, optional
        Inlet velocity [m/s], nonnegative.

    Returns
    -------
    float or ndarray
        Exit velocity [m/s].

    Notes
    -----
    No work and no heat; raises ``ValueError`` if the enthalpy rise exceeds
    the inlet kinetic energy. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(nozzle_exit_velocity(3000e3, 2950e3)), 6)
    316.227766
    """
    _finite(h_in=h_in, h_out=h_out)
    _nonnegative(V_in=V_in)
    square = np.asarray(V_in, dtype=float) ** 2 + 2 * (np.asarray(h_in, dtype=float) - h_out)
    _nonnegative(**{"V_in**2 + 2 (h_in - h_out)": square})
    return np.sqrt(square)[()]


def _ratio(numerator, denominator, name):
    _positive(**{name: denominator})
    return (np.asarray(numerator, dtype=float) / denominator)[()]


def turbine_isentropic_efficiency(h_in, h_out, h_out_s):
    r"""Isentropic efficiency of a turbine, :math:`(h_1 - h_2)/(h_1 - h_{2s})`.

    Parameters
    ----------
    h_in : float or array_like
        Inlet enthalpy [J/kg].
    h_out : float or array_like
        Actual outlet enthalpy [J/kg].
    h_out_s : float or array_like
        Outlet enthalpy of an isentropic expansion to the same pressure
        [J/kg], below ``h_in``.

    Returns
    -------
    float or ndarray
        Actual work over isentropic work [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(turbine_isentropic_efficiency(3000e3, 2600e3, 2500e3)), 6)
    0.8
    """
    _finite(h_in=h_in, h_out=h_out, h_out_s=h_out_s)
    return _ratio(np.asarray(h_in, dtype=float) - h_out, np.asarray(h_in, dtype=float) - h_out_s, "h_in - h_out_s")


def compressor_isentropic_efficiency(h_in, h_out, h_out_s):
    r"""Isentropic efficiency of a compressor or pump, :math:`(h_{2s} - h_1)/(h_2 - h_1)`.

    Parameters
    ----------
    h_in : float or array_like
        Inlet enthalpy [J/kg].
    h_out : float or array_like
        Actual outlet enthalpy [J/kg], above ``h_in``.
    h_out_s : float or array_like
        Outlet enthalpy of an isentropic compression to the same pressure
        [J/kg].

    Returns
    -------
    float or ndarray
        Isentropic work over actual work [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(compressor_isentropic_efficiency(300e3, 600e3, 540e3)), 6)
    0.8
    """
    _finite(h_in=h_in, h_out=h_out, h_out_s=h_out_s)
    return _ratio(np.asarray(h_out_s, dtype=float) - h_in, np.asarray(h_out, dtype=float) - h_in, "h_out - h_in")


def nozzle_isentropic_efficiency(h_in, h_out, h_out_s, V_in=0.0):
    r"""Isentropic efficiency of a nozzle, the ratio of exit kinetic energies :math:`V_2^2/V_{2s}^2`.

    Parameters
    ----------
    h_in : float or array_like
        Inlet enthalpy [J/kg].
    h_out : float or array_like
        Actual exit enthalpy [J/kg].
    h_out_s : float or array_like
        Exit enthalpy of an isentropic expansion to the same pressure
        [J/kg].
    V_in : float or array_like, optional
        Inlet velocity [m/s], nonnegative.

    Returns
    -------
    float or ndarray
        Actual over isentropic exit kinetic energy [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(nozzle_isentropic_efficiency(3000e3, 2955e3, 2950e3)), 6)
    0.9
    """
    actual = nozzle_exit_velocity(h_in, h_out, V_in) ** 2
    ideal = nozzle_exit_velocity(h_in, h_out_s, V_in) ** 2
    return _ratio(actual, ideal, "isentropic exit kinetic energy")


def expansion_outlet_enthalpy(h_in, h_out_s, efficiency):
    r"""Actual outlet enthalpy of a turbine, :math:`h_2 = h_1 - \eta(h_1 - h_{2s})`.

    Parameters
    ----------
    h_in : float or array_like
        Inlet enthalpy [J/kg].
    h_out_s : float or array_like
        Isentropic outlet enthalpy [J/kg].
    efficiency : float or array_like
        Isentropic efficiency [-] in (0, 1].

    Returns
    -------
    float or ndarray
        Actual outlet enthalpy [J/kg].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(expansion_outlet_enthalpy(3000e3, 2500e3, 0.8)), 6)
    2600000.0
    """
    _finite(h_in=h_in, h_out_s=h_out_s)
    _interval("efficiency", efficiency, 0, 1, closed="right")
    h_in = np.asarray(h_in, dtype=float)
    return (h_in - efficiency * (h_in - h_out_s))[()]


def compression_outlet_enthalpy(h_in, h_out_s, efficiency):
    r"""Actual outlet enthalpy of a compressor or pump, :math:`h_2 = h_1 + (h_{2s} - h_1)/\eta`.

    Parameters
    ----------
    h_in : float or array_like
        Inlet enthalpy [J/kg].
    h_out_s : float or array_like
        Isentropic outlet enthalpy [J/kg].
    efficiency : float or array_like
        Isentropic efficiency [-] in (0, 1].

    Returns
    -------
    float or ndarray
        Actual outlet enthalpy [J/kg].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> round(float(compression_outlet_enthalpy(300e3, 540e3, 0.8)), 6)
    600000.0
    """
    _finite(h_in=h_in, h_out_s=h_out_s)
    _interval("efficiency", efficiency, 0, 1, closed="right")
    h_in = np.asarray(h_in, dtype=float)
    return (h_in + (np.asarray(h_out_s, dtype=float) - h_in) / efficiency)[()]


def mixing_chamber(mass_flows, enthalpies):
    r"""Adiabatic mixing of several streams: total mass flow and mixed enthalpy.

    Parameters
    ----------
    mass_flows : array_like
        Mass flow of each inlet stream [kg/s], nonnegative, not all zero.
    enthalpies : array_like
        Specific enthalpy of each inlet stream [J/kg], same length.

    Returns
    -------
    tuple of float
        Outlet mass flow :math:`\sum \dot m_i` [kg/s] and outlet enthalpy
        :math:`\sum \dot m_i h_i / \sum \dot m_i` [J/kg].

    Notes
    -----
    No work, no heat and negligible kinetic and potential energy. Inputs
    are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> mixing_chamber([1.0, 3.0], [100e3, 300e3])
    (4.0, 250000.0)
    """
    m = _vector("mass_flows", mass_flows)
    h = _vector("enthalpies", enthalpies)
    _nonnegative(mass_flows=m)
    if m.shape != h.shape:
        raise ValueError("mass_flows and enthalpies must have the same length")
    total = float(m.sum())
    _positive(**{"total mass flow": total})
    return total, float(m @ h / total)


def throttle_steam(state, p_out):
    r"""State of water or steam after a throttling valve: same enthalpy, lower pressure.

    Parameters
    ----------
    state : SteamState
        Inlet state, e.g. from :func:`steam_properties`.
    p_out : float
        Outlet pressure [Pa], positive and at most the inlet pressure.

    Returns
    -------
    SteamState
        Outlet state with :math:`h_2 = h_1`; its temperature drop is the
        Joule-Thomson effect and its entropy rise the lost work.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/thermodynamics`.

    Examples
    --------
    >>> out = throttle_steam(steam_properties(1e6, x=0.0), 1e5)
    >>> round(out.x, 4)
    0.1529
    """
    if not isinstance(state, SteamState) or np.ndim(state.h):
        raise ValueError("state must be a scalar SteamState")
    _positive(p_out=p_out)
    if p_out > state.p:
        raise ValueError("a throttle cannot raise the pressure")
    return steam_properties(p_out, h=state.h)

r"""Pipe flow: Reynolds's regimes, laminar Hagen-Poiseuille flow and turbulent friction.

Hagen (1839) and Poiseuille (1840) measured laminar flow in small tubes,
:math:`Q \propto \Delta p\,R^4`. Weisbach (1845) and Darcy (1857) wrote the
head loss of any pipe as :math:`h_f = f\,(L/D)\,v^2/(2g)`, and Reynolds
(1883) showed that the regime, and so the friction factor, depends on
:math:`Re = \rho v D/\mu`. Nikuradse (1932-1933) measured the smooth and
fully rough limits, Colebrook (1939) joined them for commercial pipes, Moody
(1944) drew the chart, and Swamee and Jain (1976) made it explicit.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.fluid_engineering.utils.constants import LAMINAR_LIMIT, TURBULENT_LIMIT
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY as G

__all__ = [
    "reynolds_number",
    "hydraulic_diameter",
    "poiseuille_flow",
    "poiseuille_velocity",
    "laminar_friction_factor",
    "darcy_head_loss",
    "minor_head_loss",
    "pipe_resistance",
    "smooth_pipe_friction_factor",
    "nikuradse_friction_factor",
    "colebrook_friction_factor",
    "swamee_jain_friction_factor",
]


def _elementwise(solve, *args):
    """Apply a scalar root solve to broadcast arrays; return a float for scalar input."""
    arrays = np.broadcast_arrays(*(np.asarray(a, dtype=float) for a in args))
    out = np.array([solve(*values) for values in zip(*(a.ravel() for a in arrays))]).reshape(arrays[0].shape)
    return float(out) if out.ndim == 0 else out


def reynolds_number(rho, v, D, mu):
    r"""Reynolds number :math:`Re = \rho |v| D/\mu`.

    Parameters
    ----------
    rho : float or array_like
        Fluid density [kg/m^3], positive.
    v : float or array_like
        Mean flow velocity [m/s]; its magnitude is used.
    D : float or array_like
        Pipe (or hydraulic) diameter [m], positive.
    mu : float or array_like
        Dynamic viscosity [Pa s], positive.

    Returns
    -------
    float or ndarray
        Reynolds number [-]. Pipe flow is laminar below about 2300.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    References
    ----------
    O. Reynolds, "An experimental investigation of the circumstances which
    determine whether the motion of water shall be direct or sinuous, and of
    the law of resistance in parallel channels," *Phil. Trans. R. Soc.* 174
    (1883) 935-982.

    Examples
    --------
    >>> reynolds_number(1000.0, 1.0, 0.1, 0.001)
    100000.0
    """
    _finite(rho=rho, v=v, D=D, mu=mu)
    _positive(rho=rho, D=D, mu=mu)
    return rho * abs(v) * D / mu


def hydraulic_diameter(area, wetted_perimeter):
    r"""Hydraulic diameter :math:`D_h = 4A/P` of a non-circular duct.

    Parameters
    ----------
    area : float or array_like
        Flow area [m^2], positive.
    wetted_perimeter : float or array_like
        Wetted perimeter [m], positive.

    Returns
    -------
    float or ndarray
        Hydraulic diameter [m]; equal to the diameter for a full circular pipe.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> hydraulic_diameter(1.0, 4.0)
    1.0
    """
    _finite(area=area, wetted_perimeter=wetted_perimeter)
    _positive(area=area, wetted_perimeter=wetted_perimeter)
    return 4 * area / wetted_perimeter


def poiseuille_flow(delta_p, radius, length, mu):
    r"""Hagen-Poiseuille volume flow, :math:`Q = \pi R^4 \Delta p/(8\mu L)`.

    Parameters
    ----------
    delta_p : float or array_like
        Pressure drop along the pipe [Pa]; the sign sets the direction.
    radius : float or array_like
        Pipe radius [m], positive.
    length : float or array_like
        Pipe length [m], positive.
    mu : float or array_like
        Dynamic viscosity [Pa s], positive.

    Returns
    -------
    float or ndarray
        Volume flow [m^3/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Fully developed laminar flow of a
    Newtonian fluid with no slip at the wall.

    References
    ----------
    G. Hagen, "Über die Bewegung des Wassers in engen cylindrischen Röhren,"
    *Ann. Phys. Chem.* 46 (1839) 423-442; J. L. M. Poiseuille, "Recherches
    expérimentales sur le mouvement des liquides dans les tubes de très
    petits diamètres," *C. R. Acad. Sci.* 11 (1840) 961-967, 1041-1048.

    Examples
    --------
    >>> f"{poiseuille_flow(100.0, 0.001, 1.0, 0.001):.4e}"
    '3.9270e-08'
    """
    _finite(delta_p=delta_p, radius=radius, length=length, mu=mu)
    _positive(radius=radius, length=length, mu=mu)
    return np.pi * radius**4 * delta_p / (8 * mu * length)


def poiseuille_velocity(r, delta_p, radius, length, mu):
    r"""Parabolic Hagen-Poiseuille velocity profile, :math:`u(r) = \Delta p\,(R^2 - r^2)/(4\mu L)`.

    Parameters
    ----------
    r : float or array_like
        Radial position [m], between 0 and ``radius``.
    delta_p : float or array_like
        Pressure drop along the pipe [Pa].
    radius : float or array_like
        Pipe radius [m], positive.
    length : float or array_like
        Pipe length [m], positive.
    mu : float or array_like
        Dynamic viscosity [Pa s], positive.

    Returns
    -------
    float or ndarray
        Axial velocity [m/s]; the centreline value is twice the mean.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> float(poiseuille_velocity(0.0, 4.0, 1.0, 1.0, 1.0))
    1.0
    """
    _finite(r=r, delta_p=delta_p, radius=radius, length=length, mu=mu)
    _positive(radius=radius, length=length, mu=mu)
    _nonnegative(r=r)
    if np.any(np.asarray(r) > np.asarray(radius)):
        raise ValueError("r must not exceed the pipe radius")
    return delta_p * (radius * radius - np.asarray(r) ** 2) / (4 * mu * length)


def laminar_friction_factor(Re):
    r"""Darcy friction factor of laminar pipe flow, :math:`f = 64/Re`.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number [-], positive and below 2300.

    Returns
    -------
    float or ndarray
        Darcy friction factor [-] (four times the Fanning factor).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Follows from the Hagen-Poiseuille
    profile and the Darcy-Weisbach equation.

    Examples
    --------
    >>> laminar_friction_factor(1000)
    0.064
    """
    _finite(Re=Re)
    _positive(Re=Re)
    if np.any(np.asarray(Re) >= LAMINAR_LIMIT):
        raise ValueError("laminar model requires Re < 2300")
    return 64 / Re


def darcy_head_loss(f, L, D, v, g=G):
    r"""Darcy-Weisbach friction head loss, :math:`h_f = f\,(L/D)\,v^2/(2g)`.

    Parameters
    ----------
    f : float or array_like
        Darcy friction factor [-], nonnegative.
    L : float or array_like
        Pipe length [m], nonnegative.
    D : float or array_like
        Pipe diameter [m], positive.
    v : float or array_like
        Mean velocity [m/s].
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Head loss [m]; multiply by :math:`\rho g` for the pressure drop.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    References
    ----------
    J. Weisbach, *Lehrbuch der Ingenieur- und Maschinen-Mechanik*
    (Braunschweig: Vieweg, 1845); H. Darcy, *Recherches expérimentales
    relatives au mouvement de l'eau dans les tuyaux* (Paris:
    Mallet-Bachelier, 1857).

    Examples
    --------
    >>> darcy_head_loss(0.02, 100.0, 0.1, 2.0, g=10.0)
    4.0
    """
    _finite(f=f, L=L, D=D, v=v, g=g)
    _nonnegative(f=f, L=L)
    _positive(D=D, g=g)
    return f * (L / D) * v * v / (2 * g)


def minor_head_loss(K, v, g=G):
    r"""Head loss in a fitting, :math:`h_m = K v^2/(2g)`.

    Parameters
    ----------
    K : float or array_like
        Loss coefficient of the fitting [-], nonnegative (about 0.5 for a
        sharp entrance, 1 for an exit, 0.3-1.5 for an elbow).
    v : float or array_like
        Mean velocity [m/s].
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Head loss [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    Examples
    --------
    >>> minor_head_loss(0.5, 2.0, g=10.0)
    0.1
    """
    _finite(K=K, v=v, g=g)
    _nonnegative(K=K)
    _positive(g=g)
    return K * v * v / (2 * g)


def pipe_resistance(f, L, D, g=G):
    r"""Resistance coefficient :math:`r = 8fL/(g\pi^2D^5)`, so that :math:`h_f = r\,Q|Q|`.

    Parameters
    ----------
    f : float or array_like
        Darcy friction factor [-], nonnegative.
    L : float or array_like
        Pipe length [m], nonnegative.
    D : float or array_like
        Pipe diameter [m], positive.
    g : float or array_like, optional
        Gravitational acceleration [m/s^2], positive.

    Returns
    -------
    float or ndarray
        Resistance coefficient [s^2/m^5].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The Darcy-Weisbach loss written in
    terms of volume flow, as pipe-network methods need it; ``f`` is held
    fixed.

    Examples
    --------
    >>> round(float(pipe_resistance(0.02, 100.0, 0.1)), 1)
    16531.0
    """
    _finite(f=f, L=L, D=D, g=g)
    _nonnegative(f=f, L=L)
    _positive(D=D, g=g)
    return 8 * f * L / (g * np.pi**2 * np.asarray(D, dtype=float) ** 5)


def smooth_pipe_friction_factor(Re):
    r"""Prandtl-Kármán friction factor of a hydraulically smooth pipe.

    Solves :math:`1/\sqrt f = 2\log_{10}(Re\sqrt f) - 0.8`, the fit to
    Nikuradse's smooth-pipe measurements.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number [-], at least 4000.

    Returns
    -------
    float or ndarray
        Darcy friction factor [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    References
    ----------
    J. Nikuradse, "Gesetzmäßigkeiten der turbulenten Strömung in glatten
    Rohren," *VDI-Forschungsheft* 356 (1932).

    Examples
    --------
    >>> round(smooth_pipe_friction_factor(1e5), 4)
    0.018
    """
    _finite(Re=Re)
    _interval("Re", Re, TURBULENT_LIMIT, np.inf)

    def solve(re):
        return optimize.brentq(lambda f: 1 / np.sqrt(f) - 2 * np.log10(re * np.sqrt(f)) + 0.8, 1e-4, 1.0, xtol=1e-15)

    return _elementwise(solve, Re)


def nikuradse_friction_factor(relative_roughness):
    r"""Fully rough friction factor from Nikuradse's sand-grain pipes.

    .. math::

        \frac{1}{\sqrt f} = 1.74 + 2\log_{10}\frac{R}{k_s}
        = 1.74 + 2\log_{10}\frac{1}{2\,k_s/D}

    Parameters
    ----------
    relative_roughness : float or array_like
        Sand-grain roughness over diameter, :math:`k_s/D` [-], in (0, 0.05].
        Nikuradse covered 0.001 to 0.033.

    Returns
    -------
    float or ndarray
        Darcy friction factor [-], independent of the Reynolds number.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. Valid once the roughness Reynolds
    number :math:`Re\,(k_s/D)\sqrt{f/8}` exceeds about 70.

    References
    ----------
    J. Nikuradse, "Strömungsgesetze in rauhen Rohren," *VDI-Forschungsheft*
    361 (1933).

    Examples
    --------
    >>> round(float(nikuradse_friction_factor(0.01)), 4)
    0.0379
    """
    _interval("relative_roughness", relative_roughness, 0, 0.05, closed="right")
    return (1.74 + 2 * np.log10(1 / (2 * np.asarray(relative_roughness, dtype=float)))) ** -2


def colebrook_friction_factor(Re, relative_roughness):
    r"""Colebrook's friction factor for commercial pipes.

    Solves :math:`1/\sqrt f = -2\log_{10}\left(\frac{\varepsilon/D}{3.7} + \frac{2.51}{Re\sqrt f}\right)`
    for turbulent flow, and returns :math:`64/Re` for laminar flow.

    Parameters
    ----------
    Re : float or array_like
        Reynolds number [-], positive; the transition range 2300-4000 is rejected.
    relative_roughness : float or array_like
        Roughness over diameter :math:`\varepsilon/D` [-], in [0, 0.05].

    Returns
    -------
    float or ndarray
        Darcy friction factor [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`. The implicit equation is solved by
    bracketed root finding (:func:`scipy.optimize.brentq`), element by
    element.

    References
    ----------
    C. F. Colebrook, "Turbulent flow in pipes, with particular reference to
    the transition region between the smooth and rough pipe laws," *J. Inst.
    Civil Eng.* 11 (1939) 133-156.

    Examples
    --------
    >>> round(colebrook_friction_factor(1e5, 1e-4), 6)
    0.018514
    """
    _finite(Re=Re, relative_roughness=relative_roughness)
    _positive(Re=Re)
    _interval("relative_roughness", relative_roughness, 0, 0.05)
    re = np.asarray(Re, dtype=float)
    if np.any((re >= LAMINAR_LIMIT) & (re < TURBULENT_LIMIT)):
        raise ValueError("transition regime requires a separate model")

    def solve(re, eps):
        if re < LAMINAR_LIMIT:
            return 64 / re

        def residual(f):
            return 1 / np.sqrt(f) + 2 * np.log10(eps / 3.7 + 2.51 / (re * np.sqrt(f)))

        return optimize.brentq(residual, 1e-4, 1.0, xtol=1e-15)

    return _elementwise(solve, Re, relative_roughness)


def swamee_jain_friction_factor(Re, relative_roughness):
    r"""Swamee and Jain's explicit approximation to Colebrook's equation.

    .. math::

        f = \frac{0.25}{\left[\log_{10}\left(\frac{\varepsilon/D}{3.7} + \frac{5.74}{Re^{0.9}}\right)\right]^2}

    Parameters
    ----------
    Re : float or array_like
        Reynolds number [-], in [5000, 1e8].
    relative_roughness : float or array_like
        Roughness over diameter :math:`\varepsilon/D` [-], in [0, 0.05].

    Returns
    -------
    float or ndarray
        Darcy friction factor [-], within 3% of Colebrook in that range
        (largest near Re = 5000) and usually within 1%.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/fluid_engineering`.

    References
    ----------
    P. K. Swamee and A. K. Jain, "Explicit equations for pipe-flow
    problems," *J. Hydraul. Div. ASCE* 102 (1976) 657-664.

    Examples
    --------
    >>> round(float(swamee_jain_friction_factor(1e5, 1e-4)), 4)
    0.0185
    """
    _finite(Re=Re, relative_roughness=relative_roughness)
    _interval("Re", Re, 5000, 1e8)
    _interval("relative_roughness", relative_roughness, 0, 0.05)
    return 0.25 / (np.log10(relative_roughness / 3.7 + 5.74 / Re**0.9) ** 2)

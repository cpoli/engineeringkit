r"""Transmission lines: the telegrapher's equations, reflections and the Smith chart.

Heaviside described a line by its series resistance and inductance and its
shunt conductance and capacitance per metre. A wave travels along it with
the propagation constant :math:`\gamma = \sqrt{(R + j\omega L)(G + j\omega C)}`
and sees the characteristic impedance :math:`Z_0`. A load that differs from
:math:`Z_0` reflects part of the wave, which Smith's chart lets an engineer
follow graphically.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "transmission_line_z0",
    "transmission_line_velocity",
    "propagation_constant",
    "distortionless_condition",
    "reflection_coefficient",
    "standing_wave_ratio",
    "input_impedance",
]


def transmission_line_z0(L_per_m, C_per_m):
    r"""Characteristic impedance :math:`Z_0 = \sqrt{L/C}` of a lossless line.

    Parameters
    ----------
    L_per_m : float or array_like
        Series inductance per unit length [H/m], positive.
    C_per_m : float or array_like
        Shunt capacitance per unit length [F/m], positive.

    Returns
    -------
    float or ndarray
        Characteristic impedance [ohm].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> round(float(transmission_line_z0(250e-9, 100e-12)), 1)
    50.0
    """
    _finite(L_per_m=L_per_m, C_per_m=C_per_m)
    _positive(L_per_m=L_per_m, C_per_m=C_per_m)
    return np.sqrt(L_per_m / C_per_m)


def transmission_line_velocity(L_per_m, C_per_m):
    r"""Wave velocity :math:`v = 1/\sqrt{LC}` on a lossless line.

    Parameters
    ----------
    L_per_m : float or array_like
        Series inductance per unit length [H/m], positive.
    C_per_m : float or array_like
        Shunt capacitance per unit length [F/m], positive.

    Returns
    -------
    float or ndarray
        Phase velocity [m/s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> round(float(transmission_line_velocity(250e-9, 100e-12)) / 1e8, 3)
    2.0
    """
    _finite(L_per_m=L_per_m, C_per_m=C_per_m)
    _positive(L_per_m=L_per_m, C_per_m=C_per_m)
    return 1 / np.sqrt(L_per_m * C_per_m)


def propagation_constant(R_per_m, L_per_m, G_per_m, C_per_m, omega):
    r"""Propagation constant :math:`\gamma = \alpha + j\beta = \sqrt{(R + j\omega L)(G + j\omega C)}`.

    A sinusoidal wave on the line varies as :math:`e^{-\gamma x}`: its
    amplitude decays by :math:`\alpha` nepers per metre and its phase
    advances by :math:`\beta` radians per metre, so it travels at
    :math:`\omega/\beta`.

    Parameters
    ----------
    R_per_m : float
        Series resistance per unit length [ohm/m], nonnegative.
    L_per_m : float
        Series inductance per unit length [H/m], positive.
    G_per_m : float
        Shunt conductance per unit length [S/m], nonnegative.
    C_per_m : float
        Shunt capacitance per unit length [F/m], positive.
    omega : float or array_like
        Angular frequency [rad/s], positive.

    Returns
    -------
    complex or ndarray of complex
        :math:`\gamma` [1/m], with :math:`\alpha, \beta \ge 0`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    O. Heaviside, "On the extra current," *Phil. Mag.* 2 (1876) 135-145.

    Examples
    --------
    >>> complex(np.round(propagation_constant(0, 1e-6, 0, 1e-12, 1e6), 12))
    0.001j
    """
    _nonnegative(R_per_m=R_per_m, G_per_m=G_per_m)
    _positive(L_per_m=L_per_m, C_per_m=C_per_m, omega=omega)
    w = np.asarray(omega, dtype=float)
    return np.sqrt((R_per_m + 1j * w * L_per_m) * (G_per_m + 1j * w * C_per_m))


def distortionless_condition(R_per_m, L_per_m, G_per_m, C_per_m, rtol=1e-9):
    r"""Whether a line meets Heaviside's distortionless condition :math:`R/L = G/C`.

    On such a line the attenuation :math:`\alpha = \sqrt{RG}` and the
    velocity :math:`1/\sqrt{LC}` are the same at every frequency, so a pulse
    shrinks without changing shape.

    Parameters
    ----------
    R_per_m : float
        Series resistance per unit length [ohm/m], nonnegative.
    L_per_m : float
        Series inductance per unit length [H/m], positive.
    G_per_m : float
        Shunt conductance per unit length [S/m], nonnegative.
    C_per_m : float
        Shunt capacitance per unit length [F/m], positive.
    rtol : float, optional
        Relative tolerance on :math:`RC = GL` [-], nonnegative.

    Returns
    -------
    bool
        ``True`` if the line is distortionless.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    O. Heaviside, "Electromagnetic induction and its propagation," *The
    Electrician* (1887); reprinted in *Electrical Papers*, vol. 2 (London:
    Macmillan, 1892).

    Examples
    --------
    >>> distortionless_condition(0.01, 1e-6, 1e-8, 1e-12)
    True
    """
    _nonnegative(R_per_m=R_per_m, G_per_m=G_per_m, rtol=rtol)
    _positive(L_per_m=L_per_m, C_per_m=C_per_m)
    return bool(np.isclose(R_per_m * C_per_m, G_per_m * L_per_m, rtol=rtol, atol=0))


def reflection_coefficient(z_load, z0):
    r"""Voltage reflection coefficient :math:`\Gamma = (Z_L - Z_0)/(Z_L + Z_0)`.

    Parameters
    ----------
    z_load : complex or array_like
        Load impedance [ohm], finite with nonnegative real part (a passive
        load); an open circuit is the limit :math:`\Gamma = 1`.
    z0 : float
        Characteristic impedance [ohm], positive.

    Returns
    -------
    complex or ndarray of complex
        :math:`\Gamma` [-], inside or on the unit circle.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> complex(reflection_coefficient(150.0, 50.0))
    (0.5+0j)
    """
    _finite(z_load=z_load)
    _positive(z0=z0)
    z = np.asarray(z_load, dtype=complex)
    if np.any(z.real < 0):
        raise ValueError("z_load must have a nonnegative real part")
    return (z - z0) / (z + z0)


def standing_wave_ratio(gamma):
    r"""Voltage standing-wave ratio :math:`(1 + |\Gamma|)/(1 - |\Gamma|)`.

    Parameters
    ----------
    gamma : complex or array_like
        Reflection coefficient [-], with :math:`|\Gamma| < 1`.

    Returns
    -------
    float or ndarray
        VSWR [-], one for a matched line.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> float(standing_wave_ratio(0.5))
    3.0
    """
    _finite(gamma=gamma)
    m = np.abs(np.asarray(gamma, dtype=complex))
    if np.any(m >= 1):
        raise ValueError("|gamma| must be less than one")
    return (1 + m) / (1 - m)


def input_impedance(z_load, z0, electrical_length):
    r"""Input impedance of a lossless line of electrical length :math:`\beta\ell` ending in :math:`Z_L`.

    The reflection coefficient rotates clockwise around the Smith chart,
    :math:`\Gamma_{in} = \Gamma_L e^{-2j\beta\ell}`, and the impedance is
    :math:`Z_{in} = Z_0(1 + \Gamma_{in})/(1 - \Gamma_{in})`, equivalent to
    :math:`Z_0 (Z_L + jZ_0\tan\beta\ell)/(Z_0 + jZ_L\tan\beta\ell)`.

    Parameters
    ----------
    z_load : complex or array_like
        Load impedance [ohm], finite with nonnegative real part.
    z0 : float
        Characteristic impedance [ohm], positive.
    electrical_length : float or array_like
        :math:`\beta\ell = 2\pi\ell/\lambda` [rad], nonnegative.

    Returns
    -------
    complex or ndarray of complex
        Input impedance [ohm]; very large where the line transforms the load
        into an open circuit.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    P. H. Smith, "Transmission line calculator," *Electronics* 12 (January
    1939) 29-31.

    Examples
    --------
    >>> round(float(abs(input_impedance(100.0, 50.0, np.pi / 2))), 9)
    25.0
    """
    gamma_load = reflection_coefficient(z_load, z0)
    _nonnegative(electrical_length=electrical_length)
    gamma_in = gamma_load * np.exp(-2j * np.asarray(electrical_length, dtype=float))
    with np.errstate(divide="ignore", invalid="ignore"):
        return z0 * (1 + gamma_in) / (1 - gamma_in)

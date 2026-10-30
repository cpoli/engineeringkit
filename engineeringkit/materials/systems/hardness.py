r"""Hardness: Brinell and Vickers indentation and the estimate of tensile strength.

Brinell (1900) pressed a hardened steel ball into the part and divided the
load by the curved area of the dent; Smith and Sandland's Vickers test
(1922) used a square diamond pyramid instead, so that the number no longer
depends on the load. Indentation is plastic flow under a constrained punch,
so the hardness of steels tracks their tensile strength,
:math:`S_{ut} \approx 3.45\,\mathrm{HB}` MPa.

Hardness numbers are conventional: the load in kilogram-force over the area
in mm^2. These functions take SI inputs and return the conventional number.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import positive as _positive
from engineeringkit.constants import STANDARD_GRAVITY

__all__ = ["brinell_hardness", "vickers_hardness", "brinell_tensile_strength"]


def brinell_hardness(load, ball_diameter, indent_diameter):
    r"""Brinell hardness number from the indentation diameter.

    .. math:: \mathrm{HB} = \frac{2P}{\pi D\left(D - \sqrt{D^2 - d^2}\right)}

    with :math:`P` in kgf and :math:`D`, :math:`d` in mm.

    Parameters
    ----------
    load : float or array_like
        Test load [N], positive (29.4 kN, i.e. 3000 kgf, for steel).
    ball_diameter : float or array_like
        Ball diameter :math:`D` [m], positive (10 mm standard).
    indent_diameter : float or array_like
        Diameter of the dent :math:`d` [m], positive and less than ``D``.

    Returns
    -------
    float or ndarray
        HB [kgf/mm^2, conventionally written without units].

    Notes
    -----
    Valid for :math:`0.24D < d < 0.6D`. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    J. A. Brinell, "Mémoire sur les épreuves à bille en acier," *Congrès
    International des Méthodes d'Essai des Matériaux de Construction*
    (Paris, 1900), vol. 2, 83-94.

    Examples
    --------
    >>> round(float(brinell_hardness(500 * 9.80665, 10e-3, 1.62e-3)), 1)
    241.0
    """
    _positive(load=load, ball_diameter=ball_diameter, indent_diameter=indent_diameter)
    D, d = 1e3 * np.asarray(ball_diameter), 1e3 * np.asarray(indent_diameter)
    if np.any(d >= D):
        raise ValueError("indent_diameter must be smaller than ball_diameter")
    P = np.asarray(load) / STANDARD_GRAVITY
    return 2 * P / (np.pi * D * (D - np.sqrt(D * D - d * d)))


def vickers_hardness(load, diagonal):
    r"""Vickers hardness number, :math:`\mathrm{HV} = 2P\sin(68^\circ)/d^2 \approx 1.8544\,P/d^2`.

    Parameters
    ----------
    load : float or array_like
        Test load [N], positive.
    diagonal : float or array_like
        Mean diagonal of the square indentation :math:`d` [m], positive.

    Returns
    -------
    float or ndarray
        HV [kgf/mm^2, conventionally written without units].

    Notes
    -----
    The 136° pyramid makes geometrically similar dents at every load.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> round(float(vickers_hardness(10 * 9.80665, 0.3e-3)), 1)
    206.0
    """
    _positive(load=load, diagonal=diagonal)
    P, d = np.asarray(load) / STANDARD_GRAVITY, 1e3 * np.asarray(diagonal)
    return 2 * P * np.sin(np.radians(68)) / (d * d)


def brinell_tensile_strength(HB, factor=3.45e6):
    r"""Tensile strength estimated from Brinell hardness, :math:`S_{ut} \approx 3.45\,\mathrm{HB}` MPa.

    Parameters
    ----------
    HB : float or array_like
        Brinell hardness number, positive.
    factor : float, optional
        Strength per hardness point [Pa], positive: 3.45 MPa for most
        steels (Callister, Eq. 6.20a); brass and cast iron need their own
        factors.

    Returns
    -------
    float or ndarray
        Ultimate tensile strength [Pa].

    Notes
    -----
    An empirical correlation, good to about 10% for plain carbon and
    low-alloy steels below 500 HB. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/materials`.

    Examples
    --------
    >>> float(brinell_tensile_strength(200.0)) / 1e6
    690.0
    """
    _positive(HB=HB, factor=factor)
    return factor * np.asarray(HB)

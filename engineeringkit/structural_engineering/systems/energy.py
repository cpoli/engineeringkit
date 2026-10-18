r"""Energy methods: bending strain energy, Castigliano's theorem and Maxwell-Betti reciprocity.

Castigliano's second theorem gives the deflection under a load :math:`P` as
the derivative of the strain energy, :math:`\delta = \partial U/\partial P`.
For a beam in bending this is

.. math::

    \delta = \int_0^L \frac{M}{EI}\,\frac{\partial M}{\partial P}\,dx.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import simpson

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector

__all__ = ["strain_energy_bending_constant_moment", "castigliano_deflection", "maxwell_betti_reciprocity"]


def strain_energy_bending_constant_moment(moment, length, E, I):
    r"""Bending strain energy under a constant moment, :math:`U = M^2L/(2EI)`.

    Parameters
    ----------
    moment : float or array_like
        Bending moment [N m].
    length : float or array_like
        Length [m], positive.
    E : float or array_like
        Young's modulus [Pa], positive.
    I : float or array_like
        Second moment of area [m^4], positive.

    Returns
    -------
    float or ndarray
        Strain energy [J].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(moment=moment, length=length, E=E, I=I)
    _positive(length=length, E=E, I=I)
    return np.asarray(moment) ** 2 * length / (2 * np.asarray(E) * I)


def castigliano_deflection(dU_dP=None, *, x=None, moment=None, dmoment_dP=None, E=None, I=None):
    r"""Deflection by Castigliano's second theorem.

    Called with a single number, it returns that strain-energy derivative
    :math:`\partial U/\partial P` unchanged (the original interface). Called
    with the keyword arguments, it integrates

    .. math::

        \delta = \int \frac{M(x)}{E I(x)}\,\frac{\partial M}{\partial P}(x)\,dx

    along a beam with Simpson's rule, where :math:`M` is the bending moment
    under the actual loads and :math:`\partial M/\partial P` the moment
    produced by a unit load at the point and in the direction of the
    wanted deflection (a dummy load if no real load acts there).

    Parameters
    ----------
    dU_dP : float, optional
        Strain-energy derivative [m], returned unchanged.
    x : array_like, optional
        Stations along the beam [m], strictly increasing. An odd number of
        equally spaced stations, with a station at every load point,
        integrates piecewise-cubic integrands exactly.
    moment : array_like, optional
        Bending moment at each station [N m].
    dmoment_dP : array_like, optional
        Derivative of the moment with respect to the load [m] (N m/N).
    E : float or array_like, optional
        Young's modulus [Pa], positive.
    I : float or array_like, optional
        Second moment of area [m^4], positive; may vary along the beam.

    Returns
    -------
    float
        Deflection [m] in the direction of the load ``P``.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata. Only bending energy is included; axial and shear
    energy are negligible for slender beams.

    References
    ----------
    A. Castigliano, *Théorie de l'équilibre des systèmes élastiques et ses
    applications* (Turin: Negro, 1879).
    """
    beam = (x, moment, dmoment_dP, E, I)
    if dU_dP is not None:
        if any(value is not None for value in beam):
            raise ValueError("pass either dU_dP or the beam arguments x, moment, dmoment_dP, E, I, not both")
        _finite(dU_dP=dU_dP)
        return dU_dP
    if any(value is None for value in beam):
        raise ValueError("x, moment, dmoment_dP, E and I are all required to integrate along a beam")
    xs = _vector("x", x, min_length=2)
    if np.any(np.diff(xs) <= 0):
        raise ValueError("x must be strictly increasing")
    _finite(moment=moment, dmoment_dP=dmoment_dP)
    _positive(E=E, I=I)
    integrand = np.broadcast_to(np.asarray(moment, float) * np.asarray(dmoment_dP, float) / (np.asarray(E) * I), xs.shape)
    return float(simpson(integrand, x=xs))


def maxwell_betti_reciprocity(delta_ij, delta_ji, rtol=1e-9, atol=1e-12):
    """Check Maxwell's reciprocal theorem, :math:`\\delta_{ij} = \\delta_{ji}`.

    In a linear elastic structure, the deflection at point i due to a unit
    load at point j equals the deflection at j due to a unit load at i.

    Parameters
    ----------
    delta_ij : float
        Deflection at i due to a unit load at j [m].
    delta_ji : float
        Deflection at j due to a unit load at i [m].
    rtol : float, optional
        Relative tolerance [-].
    atol : float, optional
        Absolute displacement tolerance [m].

    Returns
    -------
    bool
        Whether the two deflections agree within the tolerances.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.

    References
    ----------
    J. C. Maxwell, "On the calculation of the equilibrium and stiffness of
    frames," *Phil. Mag.* 27 (1864) 294-299; E. Betti, *Il Nuovo Cimento*
    7-8 (1872) 5-21.
    """
    _finite(delta_ij=delta_ij, delta_ji=delta_ji)
    return bool(np.isclose(delta_ij, delta_ji, rtol=rtol, atol=atol))

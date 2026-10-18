r"""Dead time: Padé approximation and Smith's predictor.

A transport delay :math:`e^{-s\theta}` adds phase lag
:math:`\omega\theta` without changing the gain, so it erodes the phase
margin of any loop closed around it. A Padé approximant replaces it with a
rational all-pass factor that the polynomial tools can handle. Smith (1957)
showed that a controller holding a delay-free model of the plant can be
tuned as if the delay were not there.
"""

from __future__ import annotations

from math import factorial

import numpy as np

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive
from engineeringkit.control_systems.systems.transfer import transfer_function
from engineeringkit.control_systems.utils.polynomials import num_den

__all__ = ["pade_delay", "smith_predictor"]


def pade_delay(delay, order=1):
    r"""Diagonal Padé approximant of a time delay :math:`e^{-s\theta}`.

    Parameters
    ----------
    delay : float
        Dead time :math:`\theta` [s], positive.
    order : int, optional
        Degree :math:`n \ge 1` of numerator and denominator.

    Returns
    -------
    scipy.signal.TransferFunction
        :math:`D(-s)/D(s)` with
        :math:`D(s) = \sum_{k=0}^{n} \frac{(2n-k)!\,n!}{(2n)!\,k!\,(n-k)!}(s\theta)^k`.

    Notes
    -----
    The approximant is all-pass, :math:`|D(-j\omega)/D(j\omega)| = 1`, and
    matches the Taylor series of :math:`e^{-s\theta}` to order :math:`2n`, so
    its phase is accurate for :math:`\omega\theta` up to about :math:`n`. It
    has :math:`n` right-half-plane zeros, so its step response starts with an
    undershoot that the true delay does not have. See
    :doc:`/guides/control_systems`.

    References
    ----------
    H. Padé, "Sur la représentation approchée d'une fonction par des fractions
    rationnelles," *Ann. Sci. École Norm. Sup.* 9 (1892) 3-93.

    Examples
    --------
    >>> D = pade_delay(2.0)
    >>> D.num.tolist(), D.den.tolist()
    ([-1.0, 1.0], [1.0, 1.0])
    """
    _positive(delay=delay)
    _integer("order", order, minimum=1)
    n = int(order)
    c = np.array([factorial(2 * n - k) * factorial(n) / (factorial(2 * n) * factorial(k) * factorial(n - k)) for k in range(n + 1)])
    c *= float(delay) ** np.arange(n + 1)
    den = c[::-1]
    num = den * (-1.0) ** np.arange(n, -1, -1)
    return transfer_function(num, den)


def smith_predictor(controller, plant, delay, order=3):
    r"""Equivalent controller of a Smith predictor around a delayed plant.

    The plant is :math:`G(s)e^{-s\theta}`. The predictor feeds back the
    output of a delay-free model :math:`G` and corrects it by the difference
    between the measured output and the delayed model output. The result is
    the single controller

    .. math::

        C_{\mathrm{SP}} = \frac{C}{1 + C G\,(1 - e^{-s\theta})},

    which makes the closed loop :math:`CGe^{-s\theta}/(1 + CG)`: the
    delay-free design, delayed.

    Parameters
    ----------
    controller : scipy.signal.lti
        Primary controller :math:`C(s)`, designed for the delay-free plant.
    plant : scipy.signal.lti
        Delay-free plant model :math:`G(s)`.
    delay : float
        Dead time :math:`\theta` [s], positive.
    order : int, optional
        Padé order used for the delay inside the predictor.

    Returns
    -------
    scipy.signal.TransferFunction
        :math:`C_{\mathrm{SP}}(s)` with :math:`e^{-s\theta}` replaced by
        :func:`pade_delay` of the given order.

    Notes
    -----
    The predictor relies on a perfect model: mismatch in :math:`G` or
    :math:`\theta` reappears in the characteristic equation. It cannot
    stabilise an unstable plant, whose mode stays in the open loop. In a
    computer the delay is a buffer of past commands, not a Padé factor. See
    :doc:`/guides/control_systems`.

    References
    ----------
    O. J. M. Smith, "Closer control of loops with dead time," *Chem. Eng.
    Prog.* 53 (1957) 217-219.

    Examples
    --------
    >>> from engineeringkit.control_systems import transfer_function
    >>> Csp = smith_predictor(transfer_function([1], [1]), transfer_function([1], [1, 1]), 1.0, order=1)
    >>> Csp.num.tolist(), Csp.den.tolist()
    ([1.0, 3.0, 2.0], [1.0, 5.0, 2.0])
    """
    cn, cd = num_den(controller)
    gn, gd = num_den(plant)
    dn, dd = num_den(pade_delay(delay, order))
    num = np.polymul(np.polymul(cn, gd), dd)
    den = np.polyadd(np.polymul(np.polymul(cd, gd), dd), np.polymul(np.polymul(cn, gn), np.polysub(dd, dn)))
    return transfer_function(num, den)

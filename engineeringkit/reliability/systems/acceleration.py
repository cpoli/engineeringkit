r"""Accelerated life testing: Arrhenius, Eyring and inverse-power-law acceleration factors.

A test at high temperature or voltage fails parts in weeks that would last
years in service, and an acceleration factor :math:`AF = L_{use}/L_{test}`
carries the result back to use conditions. Arrhenius's 1889 law of reaction
rates, :math:`k \propto e^{-E_a/(k_B T)}`, gives the thermal factor;
Eyring's 1935 transition-state theory adds a :math:`T` prefactor; and the
inverse power law :math:`L \propto V^{-n}` models voltage, pressure or
mechanical stress. Activation energies are in joules (multiply electronvolts
by ``scipy.constants.eV``) and temperatures in kelvin.
"""

from __future__ import annotations

import numpy as np
from scipy import constants

from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector

__all__ = [
    "arrhenius_acceleration_factor",
    "arrhenius_fit",
    "eyring_acceleration_factor",
    "inverse_power_law_acceleration_factor",
]


def arrhenius_acceleration_factor(activation_energy, use_temperature, test_temperature):
    r"""Arrhenius acceleration factor :math:`AF = \exp\left[\frac{E_a}{k_B}\left(\frac{1}{T_u} - \frac{1}{T_t}\right)\right]`.

    Parameters
    ----------
    activation_energy : float or array_like
        Activation energy :math:`E_a` of the failure mechanism [J],
        nonnegative; 0.7 eV is ``0.7 * scipy.constants.eV``.
    use_temperature : float or array_like
        Absolute temperature in service :math:`T_u` [K], positive.
    test_temperature : float or array_like
        Absolute test temperature :math:`T_t` [K], positive.

    Returns
    -------
    float or ndarray
        Ratio of life in service to life on test [-]; above one when the
        test is hotter.

    Notes
    -----
    Assumes the test does not switch on a failure mechanism absent in
    service. Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    S. Arrhenius, "Über die Reaktionsgeschwindigkeit bei der Inversion von
    Rohrzucker durch Säuren," *Z. Phys. Chem.* 4 (1889) 226-248.

    Examples
    --------
    >>> from scipy.constants import eV
    >>> round(float(arrhenius_acceleration_factor(0.7 * eV, 328.15, 398.15)), 1)
    77.6
    """
    _nonnegative(activation_energy=activation_energy)
    _positive(use_temperature=use_temperature, test_temperature=test_temperature)
    return np.exp(activation_energy / constants.k * (1 / np.asarray(use_temperature) - 1 / np.asarray(test_temperature)))


def arrhenius_fit(temperatures, lives):
    r"""Least-squares Arrhenius life model :math:`L = C\,e^{E_a/(k_B T)}` from lives at several temperatures.

    :math:`\ln L` is linear in :math:`1/T` with slope :math:`E_a/k_B`, so a
    straight line through the test points gives the activation energy and
    the constant; two temperatures fix it exactly.

    Parameters
    ----------
    temperatures : array_like
        Absolute test temperatures [K], positive, at least two different.
    lives : array_like
        Matching characteristic lives (for example Weibull :math:`\eta`
        or median lives) [s], positive.

    Returns
    -------
    tuple of float
        ``(activation_energy, constant)``: :math:`E_a` [J] and :math:`C` [s].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    W. Nelson, *Accelerated Testing: Statistical Models, Test Plans, and
    Data Analyses* (Wiley, New York, 1990), ch. 2.

    Examples
    --------
    >>> from scipy.constants import eV
    >>> energy, constant = arrhenius_fit([398.15, 423.15], [5000.0, 1800.0])
    >>> round(energy / eV, 3)
    0.593
    """
    temperature = _vector("temperatures", temperatures, min_length=2)
    life = _vector("lives", lives, min_length=2)
    if temperature.shape != life.shape:
        raise ValueError("temperatures and lives must have the same length")
    _positive(temperatures=temperature, lives=life)
    if np.ptp(temperature) == 0:
        raise ValueError("temperatures must not all be equal")
    slope, intercept = np.polyfit(1 / temperature, np.log(life), 1)
    return float(slope * constants.k), float(np.exp(intercept))


def eyring_acceleration_factor(activation_energy, use_temperature, test_temperature):
    r"""Eyring acceleration factor :math:`AF = \frac{T_t}{T_u}\exp\left[\frac{E_a}{k_B}\left(\frac{1}{T_u} - \frac{1}{T_t}\right)\right]`.

    Eyring's rate :math:`(k_B T/h)\,e^{-\Delta G^\ddagger/(k_B T)}` puts the
    temperature in front of the exponential, so life is
    :math:`L \propto T^{-1} e^{E_a/(k_B T)}` and the factor is the Arrhenius
    one times :math:`T_t/T_u`.

    Parameters
    ----------
    activation_energy : float or array_like
        Activation energy [J], nonnegative.
    use_temperature : float or array_like
        Absolute temperature in service [K], positive.
    test_temperature : float or array_like
        Absolute test temperature [K], positive.

    Returns
    -------
    float or ndarray
        Ratio of life in service to life on test [-].

    Notes
    -----
    Fitted to the same data, Eyring's model needs a slightly smaller
    activation energy than Arrhenius's, since the prefactor supplies part of
    the temperature dependence. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/reliability`.

    References
    ----------
    H. Eyring, "The activated complex in chemical reactions," *J. Chem.
    Phys.* 3 (1935) 107-115.

    Examples
    --------
    >>> from scipy.constants import eV
    >>> round(float(eyring_acceleration_factor(0.7 * eV, 328.15, 398.15)), 1)
    94.2
    """
    thermal = arrhenius_acceleration_factor(activation_energy, use_temperature, test_temperature)
    return thermal * np.asarray(test_temperature) / np.asarray(use_temperature)


def inverse_power_law_acceleration_factor(use_stress, test_stress, exponent):
    r"""Inverse-power-law acceleration factor :math:`AF = (V_t/V_u)^n` for life :math:`L \propto V^{-n}`.

    Parameters
    ----------
    use_stress : float or array_like
        Stress in service :math:`V_u` (voltage [V], pressure [Pa],
        mechanical stress [Pa], ...), positive.
    test_stress : float or array_like
        Test stress :math:`V_t`, positive, in the same unit.
    exponent : float or array_like
        Life-stress exponent :math:`n` [-], nonnegative.

    Returns
    -------
    float or ndarray
        Ratio of life in service to life on test [-].

    Notes
    -----
    Basquin's fatigue law is the same model with stress amplitude. Inputs
    are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/reliability`.

    References
    ----------
    W. Nelson, *Accelerated Testing: Statistical Models, Test Plans, and
    Data Analyses* (Wiley, New York, 1990), ch. 2.

    Examples
    --------
    >>> round(float(inverse_power_law_acceleration_factor(230.0, 460.0, 3.0)), 1)
    8.0
    """
    _positive(use_stress=use_stress, test_stress=test_stress)
    _nonnegative(exponent=exponent)
    return (np.asarray(test_stress) / np.asarray(use_stress)) ** exponent

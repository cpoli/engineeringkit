r"""Fatigue: stress-life, mean stress, strain-life and cumulative damage.

Basquin fitted Wöhler's endurance data with a straight line on log-log axes,
:math:`\sigma_a = \sigma_f'(2N)^b`. Goodman drew the line that trades
alternating stress against mean stress. Coffin and Manson found the same
power law for the plastic strain range of low-cycle fatigue, and Palmgren and
Miner added up the fractions of life used by blocks of different amplitude.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = ["basquin_cycles", "goodman_safety_factor", "coffin_manson_cycles", "strain_life_amplitude", "miner_damage"]


def _negative_exponent(name, value):
    _finite(**{name: value})
    if np.any(np.asarray(value) >= 0):
        raise ValueError(f"{name} must be negative")


def basquin_cycles(stress_amplitude, sigma_f, b):
    r"""Cycles to failure from Basquin's stress-life law.

    .. math:: \sigma_a = \sigma_f'(2N)^b \quad\Rightarrow\quad
              N = \tfrac{1}{2}\left(\sigma_a/\sigma_f'\right)^{1/b}

    Parameters
    ----------
    stress_amplitude : float or array_like
        Fully reversed stress amplitude :math:`\sigma_a` [Pa], positive.
    sigma_f : float or array_like
        Fatigue strength coefficient :math:`\sigma_f'` [Pa], positive.
    b : float or array_like
        Fatigue strength exponent [-], negative (typically -0.05 to -0.12).

    Returns
    -------
    float or ndarray
        Cycles to failure :math:`N` (not reversals :math:`2N`).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    O. H. Basquin, "The exponential law of endurance tests," *Proc. Am. Soc.
    Test. Mater.* 10 (1910) 625-630.

    Examples
    --------
    >>> float(basquin_cycles(500e6, 1000e6, -0.1))
    512.0
    """
    _positive(stress_amplitude=stress_amplitude, sigma_f=sigma_f)
    _negative_exponent("b", b)
    return 0.5 * (np.asarray(stress_amplitude) / sigma_f) ** (1 / np.asarray(b))


def goodman_safety_factor(sa, sm, Se, Sut):
    r"""Safety factor against fatigue on Goodman's mean-stress line.

    .. math:: \frac{1}{n} = \frac{\sigma_a}{S_e} + \frac{\sigma_m}{S_{ut}}

    Parameters
    ----------
    sa : float or array_like
        Alternating stress amplitude [Pa], nonnegative.
    sm : float or array_like
        Mean stress [Pa], nonnegative (tensile).
    Se : float or array_like
        Endurance limit or fully reversed fatigue strength [Pa], positive.
    Sut : float or array_like
        Ultimate tensile strength [Pa], positive.

    Returns
    -------
    float or ndarray
        Safety factor [-], the factor by which both stresses can grow
        proportionally before reaching the Goodman line.

    Raises
    ------
    ValueError
        If both stresses are zero.

    Notes
    -----
    Compressive mean stresses, where the line is usually replaced by
    :math:`\sigma_a = S_e`, are not accepted. Inputs are raw SI numbers with
    no dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    J. Goodman, *Mechanics Applied to Engineering* (Longmans, Green, London,
    1899).
    """
    _nonnegative(sa=sa, sm=sm)
    _positive(Se=Se, Sut=Sut)
    sa, sm = np.asarray(sa), np.asarray(sm)
    if np.any((sa == 0) & (sm == 0)):
        raise ValueError("zero demand has unbounded safety factor")
    return 1 / (sa / Se + sm / Sut)


def coffin_manson_cycles(plastic_strain_amplitude, epsilon_f, c):
    r"""Cycles to failure from the Coffin-Manson low-cycle fatigue law.

    .. math:: \frac{\Delta\varepsilon_p}{2} = \varepsilon_f'(2N)^c
              \quad\Rightarrow\quad
              N = \tfrac{1}{2}\left(\frac{\Delta\varepsilon_p/2}{\varepsilon_f'}\right)^{1/c}

    Parameters
    ----------
    plastic_strain_amplitude : float or array_like
        Plastic strain amplitude :math:`\Delta\varepsilon_p/2` [-], positive.
    epsilon_f : float or array_like
        Fatigue ductility coefficient :math:`\varepsilon_f'` [-], positive.
    c : float or array_like
        Fatigue ductility exponent [-], negative (typically -0.5 to -0.7).

    Returns
    -------
    float or ndarray
        Cycles to failure :math:`N`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    L. F. Coffin Jr., "A study of the effects of cyclic thermal stresses on a
    ductile metal," *Trans. ASME* 76 (1954) 931-950; S. S. Manson, "Behavior
    of materials under conditions of thermal stress," NACA Report 1170 (1954).

    Examples
    --------
    >>> float(coffin_manson_cycles(0.01, 0.64, -0.5))
    2048.0
    """
    _positive(plastic_strain_amplitude=plastic_strain_amplitude, epsilon_f=epsilon_f)
    _negative_exponent("c", c)
    return 0.5 * (np.asarray(plastic_strain_amplitude) / epsilon_f) ** (1 / np.asarray(c))


def strain_life_amplitude(cycles, E, sigma_f, b, epsilon_f, c):
    r"""Total strain amplitude of the Basquin-Coffin-Manson strain-life curve.

    .. math:: \frac{\Delta\varepsilon}{2} = \frac{\sigma_f'}{E}(2N)^b + \varepsilon_f'(2N)^c

    Parameters
    ----------
    cycles : float or array_like
        Cycles to failure :math:`N`, positive.
    E : float or array_like
        Young's modulus [Pa], positive.
    sigma_f : float or array_like
        Fatigue strength coefficient [Pa], positive.
    b : float or array_like
        Fatigue strength exponent [-], negative.
    epsilon_f : float or array_like
        Fatigue ductility coefficient [-], positive.
    c : float or array_like
        Fatigue ductility exponent [-], negative.

    Returns
    -------
    tuple of float or ndarray
        ``(total, elastic, plastic)`` strain amplitudes [-].

    Notes
    -----
    The elastic and plastic lines cross at the transition life
    :math:`2N_t = (\varepsilon_f' E/\sigma_f')^{1/(b - c)}`. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/materials`.
    """
    _positive(cycles=cycles, E=E, sigma_f=sigma_f, epsilon_f=epsilon_f)
    _negative_exponent("b", b)
    _negative_exponent("c", c)
    reversals = 2 * np.asarray(cycles, float)
    elastic = sigma_f / np.asarray(E) * reversals**b
    plastic = epsilon_f * reversals**c
    return elastic + plastic, elastic, plastic


def miner_damage(cycles, cycles_to_failure):
    r"""Palmgren-Miner cumulative damage :math:`D = \sum_i n_i / N_i`.

    Parameters
    ----------
    cycles : array_like
        Applied cycles :math:`n_i` at each amplitude, nonnegative.
    cycles_to_failure : array_like
        Matching lives :math:`N_i` at those amplitudes, positive.

    Returns
    -------
    float
        Damage [-]; failure is predicted at :math:`D = 1`.

    Notes
    -----
    The sum ignores the order of the blocks. Inputs are raw SI numbers with
    no dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    A. Palmgren, "Die Lebensdauer von Kugellagern," *Z. Ver. Dtsch. Ing.* 68
    (1924) 339-341; M. A. Miner, "Cumulative damage in fatigue," *J. Appl.
    Mech.* 12 (1945) A159-A164.
    """
    cycles, cycles_to_failure = np.asarray(cycles, float), np.asarray(cycles_to_failure, float)
    _nonnegative(cycles=cycles)
    _positive(cycles_to_failure=cycles_to_failure)
    if cycles.ndim != 1 or cycles.shape != cycles_to_failure.shape:
        raise ValueError("cycle counts and lives must be matching vectors")
    return float(np.sum(cycles / cycles_to_failure))

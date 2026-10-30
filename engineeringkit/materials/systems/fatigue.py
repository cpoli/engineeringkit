r"""Fatigue: stress-life, mean stress, endurance limit, strain-life and cumulative damage.

Basquin fitted Wöhler's endurance data with a straight line on log-log axes,
:math:`\sigma_a = \sigma_f'(2N)^b`. Gerber's parabola, Goodman's line and
Soderberg's line trade alternating stress against mean stress, and Marin's
factors carry a polished test bar's endurance limit over to a real part. Coffin and Manson found the same
power law for the plastic strain range of low-cycle fatigue, and Palmgren and
Miner added up the fractions of life used by blocks of different amplitude.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "basquin_cycles",
    "goodman_safety_factor",
    "gerber_safety_factor",
    "soderberg_safety_factor",
    "endurance_limit_estimate",
    "marin_surface_factor",
    "marin_size_factor",
    "marin_load_factor",
    "marin_reliability_factor",
    "marin_endurance_limit",
    "coffin_manson_cycles",
    "strain_life_amplitude",
    "miner_damage",
]

# Shigley's surface-finish constants a [MPa^-b] and b for ka = a Sut^b, Sut in MPa.
_SURFACE = {"ground": (1.58, -0.085), "machined": (4.51, -0.265), "hot-rolled": (57.7, -0.718), "forged": (272.0, -0.995)}
_LOAD = {"bending": 1.0, "axial": 0.85, "torsion": 0.59}


def _negative_exponent(name, value):
    _finite(**{name: value})
    if np.any(np.asarray(value) >= 0):
        raise ValueError(f"{name} must be negative")


def _mean_stress_demand(sa, sm):
    _nonnegative(sa=sa, sm=sm)
    sa, sm = np.asarray(sa), np.asarray(sm)
    if np.any((sa == 0) & (sm == 0)):
        raise ValueError("zero demand has unbounded safety factor")
    return sa, sm


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
    sa, sm = _mean_stress_demand(sa, sm)
    _positive(Se=Se, Sut=Sut)
    return 1 / (sa / Se + sm / Sut)


def gerber_safety_factor(sa, sm, Se, Sut):
    r"""Safety factor against fatigue on Gerber's mean-stress parabola.

    .. math:: \frac{n\sigma_a}{S_e} + \left(\frac{n\sigma_m}{S_{ut}}\right)^2 = 1
              \quad\Rightarrow\quad
              n = \frac{2}{\sigma_a/S_e + \sqrt{(\sigma_a/S_e)^2 + 4(\sigma_m/S_{ut})^2}}

    Parameters
    ----------
    sa : float or array_like
        Alternating stress amplitude [Pa], nonnegative.
    sm : float or array_like
        Mean stress [Pa], nonnegative (tensile).
    Se : float or array_like
        Endurance limit [Pa], positive.
    Sut : float or array_like
        Ultimate tensile strength [Pa], positive.

    Returns
    -------
    float or ndarray
        Safety factor [-] along a load line through the origin; never below
        Goodman's.

    Raises
    ------
    ValueError
        If both stresses are zero.

    Notes
    -----
    Gerber's parabola runs through the middle of ductile-steel test data, so
    it is a best estimate rather than a design line. Inputs are raw SI
    numbers with no dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    H. Gerber, "Bestimmung der zulässigen Spannungen in Eisen-Constructionen,"
    *Z. Bayer. Archit. Ing.-Ver.* 6 (1874) 101-110.

    Examples
    --------
    >>> round(float(gerber_safety_factor(100e6, 0.0, 200e6, 600e6)), 6)
    2.0
    """
    sa, sm = _mean_stress_demand(sa, sm)
    _positive(Se=Se, Sut=Sut)
    a, m = sa / Se, sm / Sut
    return 2 / (a + np.sqrt(a * a + 4 * m * m))


def soderberg_safety_factor(sa, sm, Se, Sy):
    r"""Safety factor on Soderberg's line, :math:`1/n = \sigma_a/S_e + \sigma_m/S_y`.

    Parameters
    ----------
    sa : float or array_like
        Alternating stress amplitude [Pa], nonnegative.
    sm : float or array_like
        Mean stress [Pa], nonnegative (tensile).
    Se : float or array_like
        Endurance limit [Pa], positive.
    Sy : float or array_like
        Yield strength [Pa], positive.

    Returns
    -------
    float or ndarray
        Safety factor [-]; the line to the yield strength also guards
        against first-cycle yielding, and is the most conservative of the
        three.

    Raises
    ------
    ValueError
        If both stresses are zero.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    C. R. Soderberg, "Factor of safety and working stress," *Trans. ASME* 52
    (1930) APM-52-2, 13-28.

    Examples
    --------
    >>> float(soderberg_safety_factor(50e6, 150e6, 200e6, 300e6))
    1.3333333333333333
    """
    sa, sm = _mean_stress_demand(sa, sm)
    _positive(Se=Se, Sy=Sy)
    return 1 / (sa / Se + sm / Sy)


def endurance_limit_estimate(Sut):
    r"""Rotating-beam endurance limit of a steel estimated from its strength.

    .. math:: S_e' = 0.5\,S_{ut} \ (S_{ut} \le 1400\ \mathrm{MPa}), \qquad S_e' = 700\ \mathrm{MPa}\ \text{above}

    Parameters
    ----------
    Sut : float or array_like
        Ultimate tensile strength [Pa], positive.

    Returns
    -------
    float or ndarray
        Endurance limit of a polished specimen [Pa].

    Notes
    -----
    For steels only; aluminium and copper alloys have no endurance limit.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> float(endurance_limit_estimate(600e6)) / 1e6
    300.0
    """
    _positive(Sut=Sut)
    return np.minimum(0.5 * np.asarray(Sut), 700e6)


def marin_surface_factor(Sut, finish):
    r"""Marin surface factor :math:`k_a = a\,S_{ut}^b` (Sut in MPa).

    Parameters
    ----------
    Sut : float or array_like
        Ultimate tensile strength [Pa], positive.
    finish : {"ground", "machined", "hot-rolled", "forged"}
        Surface finish; ``"machined"`` also covers cold-drawn.

    Returns
    -------
    float or ndarray
        :math:`k_a` [-], below one and falling with strength.

    Notes
    -----
    Constants from Shigley's *Mechanical Engineering Design*, Table 6-2.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    J. Marin, *Mechanical Behavior of Engineering Materials* (Prentice-Hall,
    Englewood Cliffs, 1962).

    Examples
    --------
    >>> round(float(marin_surface_factor(520e6, "machined")), 3)
    0.86
    """
    _positive(Sut=Sut)
    if finish not in _SURFACE:
        raise ValueError(f"finish must be one of {sorted(_SURFACE)}")
    a, b = _SURFACE[finish]
    return a * (np.asarray(Sut) / 1e6) ** b


def marin_size_factor(diameter):
    r"""Marin size factor of a round bar in rotating bending or torsion.

    .. math:: k_b = 1.24\,d^{-0.107}\ (2.79 \le d \le 51\ \mathrm{mm}), \qquad
              k_b = 1.51\,d^{-0.157}\ (51 < d \le 254\ \mathrm{mm})

    Parameters
    ----------
    diameter : float or array_like
        Bar diameter [m], between 2.79 mm and 254 mm.

    Returns
    -------
    float or ndarray
        :math:`k_b` [-]; about one at the 7.6 mm test-bar size.

    Notes
    -----
    Use :math:`k_b = 1` for axial loading. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/materials`.

    Examples
    --------
    >>> round(float(marin_size_factor(32e-3)), 3)
    0.856
    """
    _interval("diameter", diameter, 2.79e-3, 254e-3)
    d = 1e3 * np.asarray(diameter)
    return np.where(d <= 51, 1.24 * d**-0.107, 1.51 * d**-0.157)


def marin_load_factor(loading):
    r"""Marin load factor :math:`k_c`: 1 in bending, 0.85 axial, 0.59 in torsion.

    Parameters
    ----------
    loading : {"bending", "axial", "torsion"}
        Type of loading.

    Returns
    -------
    float
        :math:`k_c` [-]. The torsion value is the von Mises shear ratio
        :math:`1/\sqrt3`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    if loading not in _LOAD:
        raise ValueError(f"loading must be one of {sorted(_LOAD)}")
    return _LOAD[loading]


def marin_reliability_factor(reliability):
    r"""Marin reliability factor :math:`k_e = 1 - 0.08\,z_R`.

    Parameters
    ----------
    reliability : float or array_like
        Required probability of survival [-], in [0.5, 1).

    Returns
    -------
    float or ndarray
        :math:`k_e` [-]: 1 at 50%, 0.814 at 99%, assuming endurance limits
        scatter normally with an 8% coefficient of variation.

    Notes
    -----
    :math:`z_R` is the standard normal quantile of the reliability. Inputs
    are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> round(float(marin_reliability_factor(0.99)), 3)
    0.814
    """
    _interval("reliability", reliability, 0.5, 1, closed="left")
    return 1 - 0.08 * norm.ppf(reliability)


def marin_endurance_limit(Se_prime, ka=1.0, kb=1.0, kc=1.0, kd=1.0, ke=1.0, kf=1.0):
    r"""Endurance limit of a part, :math:`S_e = k_a k_b k_c k_d k_e k_f\,S_e'`.

    Parameters
    ----------
    Se_prime : float or array_like
        Rotating-beam endurance limit of a polished specimen [Pa], positive.
    ka, kb, kc, kd, ke, kf : float or array_like, optional
        Surface, size, load, temperature, reliability and miscellaneous
        factors [-], positive.

    Returns
    -------
    float or ndarray
        Endurance limit of the part [Pa].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    J. Marin, *Mechanical Behavior of Engineering Materials* (Prentice-Hall,
    Englewood Cliffs, 1962).

    Examples
    --------
    >>> float(marin_endurance_limit(300e6, ka=0.8, kb=0.9)) / 1e6
    216.0
    """
    _positive(Se_prime=Se_prime, ka=ka, kb=kb, kc=kc, kd=kd, ke=ke, kf=kf)
    return np.asarray(Se_prime) * ka * kb * kc * kd * ke * kf


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

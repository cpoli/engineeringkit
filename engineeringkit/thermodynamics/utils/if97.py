r"""IAPWS-IF97 basic equations for regions 1 to 5, in reduced variables.

The coefficients are those of the IAPWS release (Lucerne 2007 revision of
the 1997 formulation): region 1 is compressed liquid, region 2 superheated
vapour, region 3 the near-critical fluid, region 4 the saturation line and
region 5 the high-temperature steam of gas turbines. Regions 1, 2 and 5 are dimensionless Gibbs free energies
:math:`\gamma(\pi, \tau)` and region 3 a dimensionless Helmholtz free
energy :math:`\phi(\delta, \tau)`; every property follows from the
function and its derivatives. Pressures are in MPa and temperatures in
K inside this module, as in the release; the public wrappers in
:mod:`engineeringkit.thermodynamics.systems.steam` convert to Pa and J.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

__all__ = [
    "R_WATER",
    "T_CRIT",
    "P_CRIT",
    "region1",
    "region2",
    "region3",
    "region3_density",
    "region5",
    "saturation_pressure_mpa",
    "saturation_temperature_k",
    "b23_pressure_mpa",
    "b23_temperature_k",
]

R_WATER = 461.526
"""Specific gas constant of water used by IF97 [J/(kg K)]."""

T_CRIT = 647.096
"""Critical temperature of water [K]."""

P_CRIT = 22.064e6
"""Critical pressure of water [Pa]."""

RHO_CRIT = 322.0
"""Critical density of water [kg/m^3]."""

# Region 1 (Table 2 of the release): I, J, n.
_R1 = np.array(
    [
        (0, -2, 0.14632971213167),
        (0, -1, -0.84548187169114),
        (0, 0, -0.37563603672040e1),
        (0, 1, 0.33855169168385e1),
        (0, 2, -0.95791963387872),
        (0, 3, 0.15772038513228),
        (0, 4, -0.16616417199501e-1),
        (0, 5, 0.81214629983568e-3),
        (1, -9, 0.28319080123804e-3),
        (1, -7, -0.60706301565874e-3),
        (1, -1, -0.18990068218419e-1),
        (1, 0, -0.32529748770505e-1),
        (1, 1, -0.21841717175414e-1),
        (1, 3, -0.52838357969930e-4),
        (2, -3, -0.47184321073267e-3),
        (2, 0, -0.30001780793026e-3),
        (2, 1, 0.47661393906987e-4),
        (2, 3, -0.44141845330846e-5),
        (2, 17, -0.72694996297594e-15),
        (3, -4, -0.31679644845054e-4),
        (3, 0, -0.28270797985312e-5),
        (3, 6, -0.85205128120103e-9),
        (4, -5, -0.22425281908000e-5),
        (4, -2, -0.65171222895601e-6),
        (4, 10, -0.14341729937924e-12),
        (5, -8, -0.40516996860117e-6),
        (8, -11, -0.12734301741641e-8),
        (8, -6, -0.17424871230634e-9),
        (21, -29, -0.68762131295531e-18),
        (23, -31, 0.14478307828521e-19),
        (29, -38, 0.26335781662795e-22),
        (30, -39, -0.11947622640071e-22),
        (31, -40, 0.18228094581404e-23),
        (32, -41, -0.93537087292458e-25),
    ]
)

# Region 2 ideal-gas part (Table 10): J0, n0.
_R2_IDEAL = np.array(
    [
        (0, -0.96927686500217e1),
        (1, 0.10086655968018e2),
        (-5, -0.56087911283020e-2),
        (-4, 0.71452738081455e-1),
        (-3, -0.40710498223928),
        (-2, 0.14240819171444e1),
        (-1, -0.43839511319450e1),
        (2, -0.28408632460772),
        (3, 0.21268463753307e-1),
    ]
)

# Region 2 residual part (Table 11): I, J, n.
_R2_RES = np.array(
    [
        (1, 0, -0.17731742473213e-2),
        (1, 1, -0.17834862292358e-1),
        (1, 2, -0.45996013696365e-1),
        (1, 3, -0.57581259083432e-1),
        (1, 6, -0.50325278727930e-1),
        (2, 1, -0.33032641670203e-4),
        (2, 2, -0.18948987516315e-3),
        (2, 4, -0.39392777243355e-2),
        (2, 7, -0.43797295650573e-1),
        (2, 36, -0.26674547914087e-4),
        (3, 0, 0.20481737692309e-7),
        (3, 1, 0.43870667284435e-6),
        (3, 3, -0.32277677238570e-4),
        (3, 6, -0.15033924542148e-2),
        (3, 35, -0.40668253562649e-1),
        (4, 1, -0.78847309559367e-9),
        (4, 2, 0.12790717852285e-7),
        (4, 3, 0.48225372718507e-6),
        (5, 7, 0.22922076337661e-5),
        (6, 3, -0.16714766451061e-10),
        (6, 16, -0.21171472321355e-2),
        (6, 35, -0.23895741934104e2),
        (7, 0, -0.59059564324270e-17),
        (7, 11, -0.12621808899101e-5),
        (7, 25, -0.38946842435739e-1),
        (8, 8, 0.11256211360459e-10),
        (8, 36, -0.82311340897998e1),
        (9, 13, 0.19809712802088e-7),
        (10, 4, 0.10406965210174e-18),
        (10, 10, -0.10234747095929e-12),
        (10, 14, -0.10018179379511e-8),
        (16, 29, -0.80882908646985e-10),
        (16, 50, 0.10693031879409),
        (18, 57, -0.33662250574171),
        (20, 20, 0.89185845355421e-24),
        (20, 35, 0.30629316876232e-12),
        (20, 48, -0.42002467698208e-5),
        (21, 21, -0.59056029685639e-25),
        (22, 53, 0.37826947613457e-5),
        (23, 39, -0.12768608934681e-14),
        (24, 26, 0.73087610595061e-28),
        (24, 40, 0.55414715350778e-16),
        (24, 58, -0.94369707241210e-6),
    ]
)

# Region 3 (Table 30): n1 multiplies ln(delta); then I, J, n.
_R3_N1 = 0.10658070028513e1
_R3 = np.array(
    [
        (0, 0, -0.15732845290239e2),
        (0, 1, 0.20944396974307e2),
        (0, 2, -0.76867707878716e1),
        (0, 7, 0.26185947787954e1),
        (0, 10, -0.28080781148620e1),
        (0, 12, 0.12053369696517e1),
        (0, 23, -0.84566812812502e-2),
        (1, 2, -0.12654315477714e1),
        (1, 6, -0.11524407806681e1),
        (1, 15, 0.88521043984318),
        (1, 17, -0.64207765181607),
        (2, 0, 0.38493460186671),
        (2, 2, -0.85214708824206),
        (2, 6, 0.48972281541877e1),
        (2, 7, -0.30502617256965e1),
        (2, 22, 0.39420536879154e-1),
        (2, 26, 0.12558408424308),
        (3, 0, -0.27999329698710),
        (3, 2, 0.13899799569460e1),
        (3, 4, -0.20189915023570e1),
        (3, 16, -0.82147637173963e-2),
        (3, 26, -0.47596035734923),
        (4, 0, 0.43984074473500e-1),
        (4, 2, -0.44476435428739),
        (4, 4, 0.90572070719733),
        (4, 26, 0.70522450087967),
        (5, 1, 0.10770512626332),
        (5, 3, -0.32913623258954),
        (5, 26, -0.50871062041158),
        (6, 0, -0.22175400873096e-1),
        (6, 2, 0.94260751665092e-1),
        (6, 26, 0.16436278447961),
        (7, 2, -0.13503372241348e-1),
        (8, 26, -0.14834345352472e-1),
        (9, 2, 0.57922953628084e-3),
        (9, 26, 0.32308904703711e-2),
        (10, 0, 0.80964802996215e-4),
        (10, 1, -0.16557679795037e-3),
        (11, 26, -0.44923899061815e-4),
    ]
)

# Region 5 (Tables 37 and 38 of the 2007 revision, valid to 50 MPa): ideal part J0, n0; residual I, J, n.
_R5_IDEAL = np.array(
    [
        (0, -0.13179983674201e2),
        (1, 0.68540841634434e1),
        (-3, -0.24805148933466e-1),
        (-2, 0.36901534980333),
        (-1, -0.31161318213925e1),
        (2, -0.32961626538917),
    ]
)
_R5_RES = np.array(
    [
        (1, 1, 0.15736404855259e-2),
        (1, 2, 0.90153761673944e-3),
        (1, 3, -0.50270077677648e-2),
        (2, 3, 0.22440037409485e-5),
        (2, 9, -0.41163275453471e-5),
        (3, 7, 0.37919454822955e-7),
    ]
)

# Region 4 saturation line (Table 34): n1 ... n10.
_R4 = np.array(
    [
        0.11670521452767e4,
        -0.72421316703206e6,
        -0.17073846940092e2,
        0.12020824702470e5,
        -0.32325550322333e7,
        0.14915108613530e2,
        -0.48232657361591e4,
        0.40511340542057e6,
        -0.23855557567849,
        0.65017534844798e3,
    ]
)

# Boundary between regions 2 and 3 (Table 1): n1 ... n5.
_B23 = (0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2)


def _props(T, p_mpa, g, g_p, g_t, g_tt, pi, tau):
    """Specific properties from a reduced Gibbs function and its derivatives (SI: m^3/kg, J/kg, J/(kg K))."""
    RT = R_WATER * T
    return {
        "v": pi * g_p * RT / (p_mpa * 1e6),
        "h": RT * tau * g_t,
        "u": RT * (tau * g_t - pi * g_p),
        "s": R_WATER * (tau * g_t - g),
        "cp": -R_WATER * tau**2 * g_tt,
    }


def region1(T, p_mpa):
    """Region 1 (compressed liquid) properties from the IF97 Gibbs function.

    Parameters
    ----------
    T : float or array_like
        Temperature [K], 273.15 K to 623.15 K.
    p_mpa : float or array_like
        Pressure [MPa], above saturation and at most 100 MPa.

    Returns
    -------
    dict of ndarray
        Specific volume ``v`` [m^3/kg], enthalpy ``h`` [J/kg], internal
        energy ``u`` [J/kg], entropy ``s`` [J/(kg K)] and isobaric heat
        capacity ``cp`` [J/(kg K)], broadcast over ``T`` and ``p_mpa``.
    """
    T, p_mpa = np.broadcast_arrays(np.asarray(T, dtype=float), np.asarray(p_mpa, dtype=float))
    pi, tau = p_mpa / 16.53, 1386.0 / T
    I, J, n = _R1[:, 0], _R1[:, 1], _R1[:, 2]
    a, b = (7.1 - pi)[..., None], (tau - 1.222)[..., None]
    g = np.sum(n * a**I * b**J, axis=-1)
    g_p = np.sum(-n * I * a ** (I - 1) * b**J, axis=-1)
    g_t = np.sum(n * a**I * J * b ** (J - 1), axis=-1)
    g_tt = np.sum(n * a**I * J * (J - 1) * b ** (J - 2), axis=-1)
    return _props(T, p_mpa, g, g_p, g_t, g_tt, pi, tau)


def region2(T, p_mpa):
    """Region 2 (vapour) properties from the IF97 Gibbs function.

    Parameters
    ----------
    T : float or array_like
        Temperature [K], 273.15 K to 1073.15 K.
    p_mpa : float or array_like
        Pressure [MPa], positive and below saturation (and below the
        region 2-3 boundary above 623.15 K).

    Returns
    -------
    dict of ndarray
        Specific volume ``v`` [m^3/kg], enthalpy ``h`` [J/kg], internal
        energy ``u`` [J/kg], entropy ``s`` [J/(kg K)] and isobaric heat
        capacity ``cp`` [J/(kg K)], broadcast over ``T`` and ``p_mpa``.
    """
    T, p_mpa = np.broadcast_arrays(np.asarray(T, dtype=float), np.asarray(p_mpa, dtype=float))
    pi, tau = p_mpa / 1.0, 540.0 / T
    J0, n0 = _R2_IDEAL[:, 0], _R2_IDEAL[:, 1]
    I, J, n = _R2_RES[:, 0], _R2_RES[:, 1], _R2_RES[:, 2]
    t0, P, t = tau[..., None], pi[..., None], (tau - 0.5)[..., None]
    g = np.log(pi) + np.sum(n0 * t0**J0, axis=-1) + np.sum(n * P**I * t**J, axis=-1)
    g_p = 1.0 / pi + np.sum(n * I * P ** (I - 1) * t**J, axis=-1)
    g_t = np.sum(n0 * J0 * t0 ** (J0 - 1), axis=-1) + np.sum(n * P**I * J * t ** (J - 1), axis=-1)
    g_tt = np.sum(n0 * J0 * (J0 - 1) * t0 ** (J0 - 2), axis=-1) + np.sum(n * P**I * J * (J - 1) * t ** (J - 2), axis=-1)
    return _props(T, p_mpa, g, g_p, g_t, g_tt, pi, tau)


def region5(T, p_mpa):
    """Region 5 (high-temperature steam) properties from the IF97 Gibbs function.

    Parameters
    ----------
    T : float or array_like
        Temperature [K], 1073.15 K to 2273.15 K.
    p_mpa : float or array_like
        Pressure [MPa], positive and at most 50 MPa.

    Returns
    -------
    dict of ndarray
        Specific volume ``v`` [m^3/kg], enthalpy ``h`` [J/kg], internal
        energy ``u`` [J/kg], entropy ``s`` [J/(kg K)] and isobaric heat
        capacity ``cp`` [J/(kg K)], broadcast over ``T`` and ``p_mpa``.
    """
    T, p_mpa = np.broadcast_arrays(np.asarray(T, dtype=float), np.asarray(p_mpa, dtype=float))
    pi, tau = p_mpa / 1.0, 1000.0 / T
    J0, n0 = _R5_IDEAL[:, 0], _R5_IDEAL[:, 1]
    I, J, n = _R5_RES[:, 0], _R5_RES[:, 1], _R5_RES[:, 2]
    t, P = tau[..., None], pi[..., None]
    g = np.log(pi) + np.sum(n0 * t**J0, axis=-1) + np.sum(n * P**I * t**J, axis=-1)
    g_p = 1.0 / pi + np.sum(n * I * P ** (I - 1) * t**J, axis=-1)
    g_t = np.sum(n0 * J0 * t ** (J0 - 1), axis=-1) + np.sum(n * P**I * J * t ** (J - 1), axis=-1)
    g_tt = np.sum(n0 * J0 * (J0 - 1) * t ** (J0 - 2), axis=-1) + np.sum(n * P**I * J * (J - 1) * t ** (J - 2), axis=-1)
    return _props(T, p_mpa, g, g_p, g_t, g_tt, pi, tau)


def region3(rho, T):
    """Region 3 (near-critical) properties from the IF97 Helmholtz function.

    Parameters
    ----------
    rho : float or array_like
        Density [kg/m^3], positive.
    T : float or array_like
        Temperature [K], 623.15 K to 863.15 K.

    Returns
    -------
    dict of ndarray
        Pressure ``p`` [Pa] and, as for :func:`region1`, ``v``, ``h``,
        ``u``, ``s`` and ``cp``, broadcast over ``rho`` and ``T``.
    """
    rho, T = np.broadcast_arrays(np.asarray(rho, dtype=float), np.asarray(T, dtype=float))
    delta, tau = rho / RHO_CRIT, T_CRIT / T
    I, J, n = _R3[:, 0], _R3[:, 1], _R3[:, 2]
    d, t = delta[..., None], tau[..., None]
    f = _R3_N1 * np.log(delta) + np.sum(n * d**I * t**J, axis=-1)
    f_d = _R3_N1 / delta + np.sum(n * I * d ** (I - 1) * t**J, axis=-1)
    f_dd = -_R3_N1 / delta**2 + np.sum(n * I * (I - 1) * d ** (I - 2) * t**J, axis=-1)
    f_t = np.sum(n * d**I * J * t ** (J - 1), axis=-1)
    f_tt = np.sum(n * d**I * J * (J - 1) * t ** (J - 2), axis=-1)
    f_dt = np.sum(n * I * d ** (I - 1) * J * t ** (J - 1), axis=-1)
    RT = R_WATER * T
    coupling = delta * f_d - delta * tau * f_dt
    return {
        "p": rho * RT * delta * f_d,
        "v": 1 / rho,
        "h": RT * (tau * f_t + delta * f_d),
        "u": RT * tau * f_t,
        "s": R_WATER * (tau * f_t - f),
        "cp": R_WATER * (-(tau**2) * f_tt + coupling**2 / (2 * delta * f_d + delta**2 * f_dd)),
    }


def region3_density(T, p_mpa, phase):
    """Density at which the region 3 pressure equals ``p_mpa`` on the isotherm ``T``.

    Below the critical temperature an isotherm of the Helmholtz function has
    a van der Waals loop, so up to three densities share a pressure; the
    largest is the liquid and the smallest the vapour.

    Parameters
    ----------
    T : float
        Temperature [K], 623.15 K to 863.15 K.
    p_mpa : float
        Pressure [MPa].
    phase : {"liquid", "vapour"}
        Which root to return when there are several.

    Returns
    -------
    float
        Density [kg/m^3].

    Raises
    ------
    ValueError
        If no density between 100 and 800 kg/m^3 gives that pressure.

    Notes
    -----
    Region 3 spans about 113 to 765 kg/m^3. The search stops at 800 kg/m^3
    because the fit, extrapolated to denser states, turns back down and
    would give a spurious liquid root.
    """
    grid = np.linspace(100.0, 800.0, 2801)
    gap = region3(grid, T)["p"] / 1e6 - p_mpa
    roots = np.flatnonzero(np.sign(gap[:-1]) != np.sign(gap[1:]))
    if not roots.size:
        raise ValueError("no region 3 density gives this pressure and temperature")
    i = roots[-1] if phase == "liquid" else roots[0]
    return float(brentq(lambda r: float(region3(r, T)["p"]) / 1e6 - p_mpa, grid[i], grid[i + 1], xtol=1e-12))


def saturation_pressure_mpa(T):
    """Region 4 saturation pressure.

    Parameters
    ----------
    T : float or array_like
        Temperature [K], 273.15 K to 647.096 K.

    Returns
    -------
    float or ndarray
        Saturation pressure [MPa].
    """
    n = _R4
    T = np.asarray(T, dtype=float)
    theta = T + n[8] / (T - n[9])
    A = theta**2 + n[0] * theta + n[1]
    B = n[2] * theta**2 + n[3] * theta + n[4]
    C = n[5] * theta**2 + n[6] * theta + n[7]
    return (2 * C / (-B + np.sqrt(B**2 - 4 * A * C))) ** 4


def saturation_temperature_k(p_mpa):
    """Region 4 saturation temperature, the inverse of :func:`saturation_pressure_mpa`.

    Parameters
    ----------
    p_mpa : float or array_like
        Pressure [MPa], 611.213 Pa to 22.064 MPa.

    Returns
    -------
    float or ndarray
        Saturation temperature [K].
    """
    n = _R4
    beta = np.asarray(p_mpa, dtype=float) ** 0.25
    E = beta**2 + n[2] * beta + n[5]
    F = n[0] * beta**2 + n[3] * beta + n[6]
    G = n[1] * beta**2 + n[4] * beta + n[7]
    D = 2 * G / (-F - np.sqrt(F**2 - 4 * E * G))
    return (n[9] + D - np.sqrt((n[9] + D) ** 2 - 4 * (n[8] + n[9] * D))) / 2


def b23_pressure_mpa(T):
    """Pressure on the boundary between regions 2 and 3.

    Parameters
    ----------
    T : float or array_like
        Temperature [K], 623.15 K to 863.15 K.

    Returns
    -------
    float or ndarray
        Boundary pressure [MPa].
    """
    n1, n2, n3 = _B23
    T = np.asarray(T, dtype=float)
    return n1 + n2 * T + n3 * T**2


def b23_temperature_k(p_mpa):
    """Temperature on the boundary between regions 2 and 3, the inverse of :func:`b23_pressure_mpa`.

    Parameters
    ----------
    p_mpa : float or array_like
        Pressure [MPa], 16.529 MPa to 100 MPa.

    Returns
    -------
    float or ndarray
        Boundary temperature [K].
    """
    n4, n5 = 0.57254459862746e3, 0.13918839778870e2
    return n4 + np.sqrt((np.asarray(p_mpa, dtype=float) - n5) / _B23[2])

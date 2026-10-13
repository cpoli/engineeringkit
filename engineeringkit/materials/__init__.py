"""Engineering materials: elasticity, yielding, fracture, fatigue, creep and composites."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive


def normal_stress(force, area):
    """Normal stress.

    Parameters
    ----------
    force : N
        Input in the units stated above.
    area : m^2
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(force=force, area=area)
    _positive(area=area)
    return force / area


def normal_strain(delta_length, length):
    """Normal strain.

    Parameters
    ----------
    delta_length : m
        Input in the units stated above.
    length : m
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(delta_length=delta_length, length=length)
    _positive(length=length)
    return delta_length / length


def hooke_stress(E, strain):
    """Hooke stress.

    Parameters
    ----------
    E : Pa
        Input in the units stated above.
    strain : 1
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(E=E, strain=strain)
    _positive(E=E)
    return E * np.asarray(strain)


def shear_modulus(E, nu):
    """Shear modulus.

    Parameters
    ----------
    E : Pa
        Input in the units stated above.
    nu : 1
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(E=E, nu=nu)
    _positive(E=E)
    _interval("nu", nu, -1, 0.5, closed=False)
    return E / (2 * (1 + nu))


def bulk_modulus(E, nu):
    """Bulk modulus.

    Parameters
    ----------
    E : Pa
        Input in the units stated above.
    nu : 1
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(E=E, nu=nu)
    _positive(E=E)
    _interval("nu", nu, -1, 0.5, closed=False)
    return E / (3 * (1 - 2 * nu))


def lame_parameters(E, nu):
    """Lame parameters.

    Parameters
    ----------
    E : Pa
        Input in the units stated above.
    nu : 1
        Input in the units stated above.

    Returns
    -------
    result
        (lambda [Pa], shear modulus [Pa]).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(E=E, nu=nu)
    _positive(E=E)
    _interval("nu", nu, -1, 0.5, closed=False)
    return E * nu / ((1 + nu) * (1 - 2 * nu)), shear_modulus(E, nu)


def thermal_strain(alpha, dT):
    """Thermal strain.

    Parameters
    ----------
    alpha : 1/K
        Input in the units stated above.
    dT : K (difference)
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return alpha * dT


def constrained_thermal_stress(E, alpha, dT):
    """Constrained thermal stress.

    Parameters
    ----------
    E : Pa
        Input in the units stated above.
    alpha : 1/K
        Input in the units stated above.
    dT : K (difference)
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(E=E, alpha=alpha, dT=dT)
    _positive(E=E)
    return -E * alpha * dT


def von_mises_plane_stress(sx, sy, txy):
    """Von mises plane stress.

    Parameters
    ----------
    sx : Pa
        Input in the units stated above.
    sy : Pa
        Input in the units stated above.
    txy : Pa
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return np.sqrt(sx * sx - sx * sy + sy * sy + 3 * txy * txy)


def tresca_equivalent(principal):
    """Tresca equivalent.

    Parameters
    ----------
    principal : Pa
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    p = np.sort(np.asarray(principal, float))
    return p[-1] - p[0]


def factor_of_safety(strength, stress):
    """Factor of safety.

    Parameters
    ----------
    strength : Pa
        Input in the units stated above.
    stress : Pa
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(strength=strength, stress=stress)
    _positive(strength=strength, stress=stress)
    return strength / stress


def resilience_modulus(sigma_y, E):
    """Resilience modulus.

    Parameters
    ----------
    sigma_y : Pa
        Input in the units stated above.
    E : Pa
        Input in the units stated above.

    Returns
    -------
    result
        J/m^3.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(sigma_y=sigma_y, E=E)
    _positive(sigma_y=sigma_y, E=E)
    return sigma_y**2 / (2 * E)


def toughness(strain, stress):
    """Toughness.

    Parameters
    ----------
    strain : 1
        Input in the units stated above.
    stress : Pa
        Input in the units stated above.

    Returns
    -------
    result
        J/m^3.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    strain, stress = np.asarray(strain, float), np.asarray(stress, float)
    _finite(strain=strain, stress=stress)
    if strain.ndim != 1 or strain.shape != stress.shape or len(strain) < 2 or np.any(np.diff(strain) <= 0):
        raise ValueError("stress and increasing strain must be matching vectors")
    return float(trapezoid(stress, strain))


def stress_intensity(K_geometry, sigma, a):
    """Stress intensity.

    Parameters
    ----------
    K_geometry : 1
        Input in the units stated above.
    sigma : Pa
        Input in the units stated above.
    a : m
        Input in the units stated above.

    Returns
    -------
    result
        Pa sqrt(m).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(K_geometry=K_geometry, sigma=sigma, a=a)
    _positive(K_geometry=K_geometry, a=a)
    _nonnegative(sigma=sigma)
    return K_geometry * sigma * np.sqrt(np.pi * a)


def critical_crack_length(Kic, sigma, Y=1.0):
    """Critical crack length.

    Parameters
    ----------
    Kic : Pa sqrt(m)
        Input in the units stated above.
    sigma : Pa
        Input in the units stated above.
    Y : 1
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Kic=Kic, sigma=sigma, Y=Y)
    _positive(Kic=Kic, sigma=sigma, Y=Y)
    return (Kic / (Y * sigma)) ** 2 / np.pi


def griffith_fracture_stress(E, gamma, a):
    """Griffith fracture stress.

    Parameters
    ----------
    E : Pa
        Input in the units stated above.
    gamma : J/m^2
        Input in the units stated above.
    a : m
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(E=E, gamma=gamma, a=a)
    _positive(E=E, gamma=gamma, a=a)
    return np.sqrt(2 * E * gamma / (np.pi * a))


def paris_crack_growth(C, m, delta_K):
    """Paris crack growth.

    Parameters
    ----------
    C : m/(cycle (Pa sqrt(m))^m)
        Input in the units stated above.
    m : 1
        Input in the units stated above.
    delta_K : Pa sqrt(m)
        Input in the units stated above.

    Returns
    -------
    result
        m/cycle.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(C=C, m=m, delta_K=delta_K)
    _nonnegative(C=C, delta_K=delta_K)
    _positive(m=m)
    return C * np.asarray(delta_K) ** m


def basquin_cycles(stress_amplitude, sigma_f, b):
    """Estimate fatigue life in cycles from Basquin stress-life parameters.

    Parameters
    ----------
    stress_amplitude : Pa
        Input in the units stated above.
    sigma_f : Pa
        Input in the units stated above.
    b : 1
        Input in the units stated above.

    Returns
    -------
    result
        cycles (not reversals).

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(stress_amplitude=stress_amplitude, sigma_f=sigma_f, b=b)
    _positive(stress_amplitude=stress_amplitude, sigma_f=sigma_f)
    if b >= 0:
        raise ValueError("Basquin exponent b must be negative")
    return 0.5 * (stress_amplitude / sigma_f) ** (1 / b)


def goodman_safety_factor(sa, sm, Se, Sut):
    """Goodman safety factor.

    Parameters
    ----------
    sa : Pa
        Input in the units stated above.
    sm : Pa
        Input in the units stated above.
    Se : Pa
        Input in the units stated above.
    Sut : Pa
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(sa=sa, sm=sm, Se=Se, Sut=Sut)
    _nonnegative(sa=sa, sm=sm)
    _positive(Se=Se, Sut=Sut)
    if sa == 0 and sm == 0:
        raise ValueError("zero demand has unbounded safety factor")
    return 1 / (sa / Se + sm / Sut)


def miner_damage(cycles, cycles_to_failure):
    """Miner damage.

    Parameters
    ----------
    cycles : cycles
        Input in the units stated above.
    cycles_to_failure : cycles
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    cycles, cycles_to_failure = np.asarray(cycles, float), np.asarray(cycles_to_failure, float)
    _nonnegative(cycles=cycles)
    _positive(cycles_to_failure=cycles_to_failure)
    if cycles.ndim != 1 or cycles.shape != cycles_to_failure.shape:
        raise ValueError("cycle counts and lives must be matching vectors")
    return float(np.sum(cycles / cycles_to_failure))


def arrhenius_creep_rate(A, Q, R, T, sigma, n=1):
    """Arrhenius creep rate.

    Parameters
    ----------
    A : 1/(s Pa^n)
        Input in the units stated above.
    Q : J/mol
        Input in the units stated above.
    R : J/(mol K)
        Input in the units stated above.
    T : K (absolute)
        Input in the units stated above.
    sigma : Pa
        Input in the units stated above.
    n : 1
        Input in the units stated above.

    Returns
    -------
    result
        1/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(A=A, Q=Q, R=R, T=T, sigma=sigma, n=n)
    _nonnegative(A=A, Q=Q, sigma=sigma)
    _positive(R=R, T=T, n=n)
    return A * sigma**n * np.exp(-Q / (R * T))


def hall_petch(sigma0, k, d):
    """Hall petch.

    Parameters
    ----------
    sigma0 : Pa
        Input in the units stated above.
    k : Pa sqrt(m)
        Input in the units stated above.
    d : m
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(sigma0=sigma0, k=k, d=d)
    _nonnegative(sigma0=sigma0, k=k)
    _positive(d=d)
    return sigma0 + k / np.sqrt(d)


def rule_of_mixtures_E_longitudinal(Ef, Em, Vf):
    """Rule of mixtures e longitudinal.

    Parameters
    ----------
    Ef : Pa
        Input in the units stated above.
    Em : Pa
        Input in the units stated above.
    Vf : 1
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Ef=Ef, Em=Em, Vf=Vf)
    _positive(Ef=Ef, Em=Em)
    _interval("Vf", Vf, 0, 1)
    return Ef * Vf + Em * (1 - Vf)


def rule_of_mixtures_E_transverse(Ef, Em, Vf):
    """Rule of mixtures e transverse.

    Parameters
    ----------
    Ef : Pa
        Input in the units stated above.
    Em : Pa
        Input in the units stated above.
    Vf : 1
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Ef=Ef, Em=Em, Vf=Vf)
    _positive(Ef=Ef, Em=Em)
    _interval("Vf", Vf, 0, 1)
    return 1 / (Vf / Ef + (1 - Vf) / Em)


def density_rule_of_mixtures(rhof, rhom, Vf):
    """Density rule of mixtures.

    Parameters
    ----------
    rhof : kg/m^3
        Input in the units stated above.
    rhom : kg/m^3
        Input in the units stated above.
    Vf : 1
        Input in the units stated above.

    Returns
    -------
    result
        kg/m^3.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/materials`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(rhof=rhof, rhom=rhom, Vf=Vf)
    _positive(rhof=rhof, rhom=rhom)
    _interval("Vf", Vf, 0, 1)
    return rhof * Vf + rhom * (1 - Vf)


@dataclass(frozen=True)
class StressStrainCurve:
    strain: np.ndarray
    stress: np.ndarray

    def __post_init__(self):
        strain, stress = np.asarray(self.strain, float), np.asarray(self.stress, float)
        toughness(strain, stress)
        object.__setattr__(self, "strain", strain.copy())
        object.__setattr__(self, "stress", stress.copy())

    def young_modulus(self, n=5):
        _integer("n", n, 2)
        if n > len(self.strain):
            raise ValueError("n cannot exceed the number of samples")
        return float(np.polyfit(self.strain[:n], self.stress[:n], 1)[0])

    def toughness(self):
        return toughness(self.strain, self.stress)

    def plot(self):
        fig, ax = plt.subplots()
        ax.plot(self.strain, self.stress)
        ax.set(xlabel="strain", ylabel="stress [Pa]")
        return fig, ax

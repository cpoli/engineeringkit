"""Mechanical engineering fundamentals: motion, mechanisms, machines and vibration."""

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.integrators import integrate_ode


def kinetic_energy(m, v):
    """Kinetic energy.

    Parameters
    ----------
    m : kg
        Input in the units stated above.
    v : m/s
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(m=m, v=v)
    _positive(m=m)
    return 0.5 * m * v * v


def rotational_energy(I, omega):
    """Rotational energy.

    Parameters
    ----------
    I : kg m^2
        Input in the units stated above.
    omega : rad/s
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(I=I, omega=omega)
    _positive(I=I)
    return 0.5 * I * omega * omega


def work_constant_force(F, d, theta=0):
    """Work constant force.

    Parameters
    ----------
    F : N
        Input in the units stated above.
    d : m
        Input in the units stated above.
    theta : rad
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return F * d * np.cos(theta)


def power_torque(torque, omega):
    """Power torque.

    Parameters
    ----------
    torque : N m
        Input in the units stated above.
    omega : rad/s
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return torque * omega


def angular_momentum(I, omega):
    """Angular momentum.

    Parameters
    ----------
    I : kg m^2
        Input in the units stated above.
    omega : rad/s
        Input in the units stated above.

    Returns
    -------
    result
        kg m^2/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(I=I, omega=omega)
    _positive(I=I)
    return I * omega


def centripetal_force(m, v, r):
    """Centripetal force.

    Parameters
    ----------
    m : kg
        Input in the units stated above.
    v : m/s
        Input in the units stated above.
    r : m
        Input in the units stated above.

    Returns
    -------
    result
        N.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(m=m, v=v, r=r)
    _positive(m=m, r=r)
    return m * v * v / r


def gear_ratio(teeth_driven, teeth_driver):
    """Gear ratio.

    Parameters
    ----------
    teeth_driven : 1
        Input in the units stated above.
    teeth_driver : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(teeth_driven=teeth_driven, teeth_driver=teeth_driver)
    _integer("teeth_driven", teeth_driven, 1)
    _integer("teeth_driver", teeth_driver, 1)
    return teeth_driven / teeth_driver


def gear_output_speed(input_speed, teeth_driver, teeth_driven):
    """Gear output speed.

    Parameters
    ----------
    input_speed : rad/s
        Input in the units stated above.
    teeth_driver : 1
        Input in the units stated above.
    teeth_driven : 1
        Input in the units stated above.

    Returns
    -------
    result
        rad/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(input_speed=input_speed, teeth_driver=teeth_driver, teeth_driven=teeth_driven)
    _integer("teeth_driven", teeth_driven, 1)
    _integer("teeth_driver", teeth_driver, 1)
    return input_speed * teeth_driver / teeth_driven


def gear_output_torque(input_torque, teeth_driver, teeth_driven, efficiency=1):
    """Gear output torque.

    Parameters
    ----------
    input_torque : N m
        Input in the units stated above.
    teeth_driver : 1
        Input in the units stated above.
    teeth_driven : 1
        Input in the units stated above.
    efficiency : 1
        Input in the units stated above.

    Returns
    -------
    result
        N m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(input_torque=input_torque, teeth_driver=teeth_driver, teeth_driven=teeth_driven, efficiency=efficiency)
    _integer("teeth_driven", teeth_driven, 1)
    _integer("teeth_driver", teeth_driver, 1)
    _interval("efficiency", efficiency, 0, 1)
    return input_torque * teeth_driven / teeth_driver * efficiency


def belt_speed(diameter, rpm):
    """Belt speed.

    Parameters
    ----------
    diameter : m
        Input in the units stated above.
    rpm : rev/min
        Input in the units stated above.

    Returns
    -------
    result
        m/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(diameter=diameter, rpm=rpm)
    _positive(diameter=diameter)
    return np.pi * diameter * rpm / 60


def bearing_l10_life(C, P, p=3):
    """Bearing l10 life.

    Parameters
    ----------
    C : N
        Input in the units stated above.
    P : N
        Input in the units stated above.
    p : 1
        Input in the units stated above.

    Returns
    -------
    result
        revolutions at 90% survival.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(C=C, P=P, p=p)
    _positive(C=C, P=P, p=p)
    return (C / P) ** p * 1e6


def screw_torque_raise(load, mean_diameter, lead, mu):
    """Screw torque raise.

    Parameters
    ----------
    load : N
        Input in the units stated above.
    mean_diameter : m
        Input in the units stated above.
    lead : m
        Input in the units stated above.
    mu : 1
        Input in the units stated above.

    Returns
    -------
    result
        N m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(load=load, mean_diameter=mean_diameter, lead=lead, mu=mu)
    _nonnegative(load=load, mu=mu)
    _positive(mean_diameter=mean_diameter, lead=lead)
    lam = np.arctan(lead / (np.pi * mean_diameter))
    phi = np.arctan(mu)
    if phi + lam >= np.pi / 2:
        raise ValueError("friction and lead angles must sum to less than pi/2")
    return load * mean_diameter / 2 * np.tan(phi + lam)


def flywheel_energy(I, w1, w2):
    """Flywheel energy.

    Parameters
    ----------
    I : kg m^2
        Input in the units stated above.
    w1 : rad/s
        Input in the units stated above.
    w2 : rad/s
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(I=I, w1=w1, w2=w2)
    _positive(I=I)
    return 0.5 * I * (w1 * w1 - w2 * w2)


def natural_frequency(k, m):
    """Natural frequency.

    Parameters
    ----------
    k : N/m
        Input in the units stated above.
    m : kg
        Input in the units stated above.

    Returns
    -------
    result
        rad/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(k=k, m=m)
    _positive(k=k, m=m)
    return np.sqrt(k / m)


def natural_frequency_hz(k, m):
    """Natural frequency hz.

    Parameters
    ----------
    k : N/m
        Input in the units stated above.
    m : kg
        Input in the units stated above.

    Returns
    -------
    result
        Hz.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return natural_frequency(k, m) / (2 * np.pi)


def damping_ratio(c, k, m):
    """Damping ratio.

    Parameters
    ----------
    c : N s/m
        Input in the units stated above.
    k : N/m
        Input in the units stated above.
    m : kg
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(c=c, k=k, m=m)
    _nonnegative(c=c)
    _positive(k=k, m=m)
    return c / (2 * np.sqrt(k * m))


def damped_frequency(k, m, c):
    """Return underdamped frequency, or zero when no oscillation exists.

    Parameters
    ----------
    k : N/m
        Input in the units stated above.
    m : kg
        Input in the units stated above.
    c : N s/m
        Input in the units stated above.

    Returns
    -------
    result
        rad/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    wn = natural_frequency(k, m)
    z = damping_ratio(c, k, m)
    return wn * np.sqrt(max(0, 1 - z * z))


def transmissibility(r, zeta):
    """Transmissibility.

    Parameters
    ----------
    r : 1 (excitation/natural frequency)
        Input in the units stated above.
    zeta : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _nonnegative(r=r, zeta=zeta)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sqrt(np.divide(1 + (2 * zeta * r) ** 2, (1 - r * r) ** 2 + (2 * zeta * r) ** 2))


def rotating_unbalance_force(m_e, e, omega):
    """Rotating unbalance force.

    Parameters
    ----------
    m_e : kg
        Input in the units stated above.
    e : m
        Input in the units stated above.
    omega : rad/s
        Input in the units stated above.

    Returns
    -------
    result
        N.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(m_e=m_e, e=e, omega=omega)
    _nonnegative(m_e=m_e, e=e)
    return m_e * e * omega**2


def slider_crank_displacement(r, l, theta):
    """Slider crank displacement.

    Parameters
    ----------
    r : m
        Input in the units stated above.
    l : m
        Input in the units stated above.
    theta : rad
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(r=r, l=l, theta=theta)
    _positive(r=r, l=l)
    if l < r:
        raise ValueError("connecting rod l must be >= crank radius r")
    return r * np.cos(theta) + np.sqrt(l * l - r * r * np.sin(theta) ** 2)


def fourbar_grashof(lengths):
    """Fourbar grashof.

    Parameters
    ----------
    lengths : m
        Input in the units stated above.

    Returns
    -------
    result
        bool.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _positive(lengths=lengths)
    if np.asarray(lengths).ndim != 1 or len(lengths) != 4:
        raise ValueError("four positive link lengths are required")
    a = sorted(lengths)
    return a[0] + a[-1] <= a[1] + a[2]


def disk_inertia(m, r):
    """Disk inertia.

    Parameters
    ----------
    m : kg
        Input in the units stated above.
    r : m
        Input in the units stated above.

    Returns
    -------
    result
        kg m^2.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(m=m, r=r)
    _positive(m=m, r=r)
    return 0.5 * m * r * r


def solid_sphere_inertia(m, r):
    """Solid sphere inertia.

    Parameters
    ----------
    m : kg
        Input in the units stated above.
    r : m
        Input in the units stated above.

    Returns
    -------
    result
        kg m^2.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(m=m, r=r)
    _positive(m=m, r=r)
    return 0.4 * m * r * r


def parallel_axis(I_centroid, m, d):
    """Parallel axis.

    Parameters
    ----------
    I_centroid : kg m^2
        Input in the units stated above.
    m : kg
        Input in the units stated above.
    d : m
        Input in the units stated above.

    Returns
    -------
    result
        kg m^2.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(I_centroid=I_centroid, m=m, d=d)
    _nonnegative(I_centroid=I_centroid)
    _positive(m=m)
    return I_centroid + m * d * d


@dataclass(frozen=True)
class SDOFResult:
    t: np.ndarray
    displacement: np.ndarray
    velocity: np.ndarray

    def summary(self):
        return {"peak_displacement_m": float(np.max(np.abs(self.displacement))), "final_velocity_m_s": float(self.velocity[-1])}

    @property
    def history(self):
        return np.column_stack((self.displacement, self.velocity))

    def plot(self):
        from engineeringkit.visualizers import plot_history

        return plot_history(self.t, self.displacement, labels=["displacement [m]"])

    def animate(self, interval=30):
        from engineeringkit.visualizers import animate_displacement

        return animate_displacement(self.t, self.displacement, interval=interval)


def simulate_sdof(m, c, k, x0=0, v0=0, t=None, force=None):
    """Integrate m*x'' + c*x' + k*x = force(t) for a passive linear oscillator.

        SI inputs: mass [kg], damping [N s/m], stiffness [N/m], displacement [m],
        velocity [m/s], time [s], and external force [N]. Constant coefficients
        and small linear displacements are assumed.

    Parameters
    ----------
    m : kg
        Input in the units stated above.
    c : N s/m
        Input in the units stated above.
    k : N/m
        Input in the units stated above.
    x0 : m
        Input in the units stated above.
    v0 : m/s
        Input in the units stated above.
    t : s
        Input in the units stated above.
    force : callable(time [s]) -> force [N] or None
        Input in the units stated above.

    Returns
    -------
    result
        SDOFResult: time [s], displacement [m], velocity [m/s].

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/mechanical_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    if not np.all(np.isfinite([m, c, k, x0, v0])) or m <= 0 or c < 0 or k < 0:
        raise ValueError("mass must be positive; damping/stiffness nonnegative; inputs finite")
    t = np.linspace(0, 10, 1001) if t is None else np.asarray(t)

    def rhs(tt, y):
        f = 0 if force is None else force(tt)
        return [y[1], (f - c * y[1] - k * y[0]) / m]

    result = integrate_ode(rhs, [x0, v0], t)
    return SDOFResult(result.t, result.history[:, 0], result.history[:, 1])

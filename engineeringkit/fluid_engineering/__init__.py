"""Fluid statics, internal flows, dimensional analysis, pumps and measurement."""

import numpy as np
from scipy import optimize

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

G = 9.80665


def hydrostatic_pressure(rho, h, g=G):
    """Hydrostatic pressure.

    Parameters
    ----------
    rho : kg/m^3
        Input in the units stated above.
    h : m
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(rho=rho, h=h, g=g)
    _positive(rho=rho, g=g)
    _nonnegative(h=h)
    return rho * g * h


def buoyancy_force(rho, displaced_volume, g=G):
    """Buoyancy force.

    Parameters
    ----------
    rho : kg/m^3
        Input in the units stated above.
    displaced_volume : m^3
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        N.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(rho=rho, displaced_volume=displaced_volume, g=g)
    _positive(rho=rho, g=g)
    _nonnegative(displaced_volume=displaced_volume)
    return rho * g * displaced_volume


def bernoulli_pressure(p1, rho, v1, v2, z1=0, z2=0, g=G):
    """Bernoulli pressure.

    Parameters
    ----------
    p1 : Pa
        Input in the units stated above.
    rho : kg/m^3
        Input in the units stated above.
    v1 : m/s
        Input in the units stated above.
    v2 : m/s
        Input in the units stated above.
    z1 : m
        Input in the units stated above.
    z2 : m
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(p1=p1, rho=rho, v1=v1, v2=v2, z1=z1, z2=z2, g=g)
    _positive(rho=rho, g=g)
    return p1 + 0.5 * rho * (v1 * v1 - v2 * v2) + rho * g * (z1 - z2)


def reynolds_number(rho, v, D, mu):
    """Reynolds number.

    Parameters
    ----------
    rho : kg/m^3
        Input in the units stated above.
    v : m/s
        Input in the units stated above.
    D : m
        Input in the units stated above.
    mu : Pa s
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(rho=rho, v=v, D=D, mu=mu)
    _positive(rho=rho, D=D, mu=mu)
    return rho * abs(v) * D / mu


def poiseuille_flow(delta_p, radius, length, mu):
    """Poiseuille flow.

    Parameters
    ----------
    delta_p : Pa
        Input in the units stated above.
    radius : m
        Input in the units stated above.
    length : m
        Input in the units stated above.
    mu : Pa s
        Input in the units stated above.

    Returns
    -------
    result
        m^3/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(delta_p=delta_p, radius=radius, length=length, mu=mu)
    _positive(radius=radius, length=length, mu=mu)
    return np.pi * radius**4 * delta_p / (8 * mu * length)


def darcy_head_loss(f, L, D, v, g=G):
    """Darcy head loss.

    Parameters
    ----------
    f : 1
        Input in the units stated above.
    L : m
        Input in the units stated above.
    D : m
        Input in the units stated above.
    v : m/s
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(f=f, L=L, D=D, v=v, g=g)
    _nonnegative(f=f, L=L)
    _positive(D=D, g=g)
    return f * (L / D) * v * v / (2 * g)


def laminar_friction_factor(Re):
    """Laminar friction factor.

    Parameters
    ----------
    Re : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Re=Re)
    _positive(Re=Re)
    if Re >= 2300:
        raise ValueError("laminar model requires Re < 2300")
    return 64 / Re


def swamee_jain_friction_factor(Re, relative_roughness):
    """Swamee jain friction factor.

    Parameters
    ----------
    Re : 1
        Input in the units stated above.
    relative_roughness : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Re=Re, relative_roughness=relative_roughness)
    _interval("Re", Re, 5000, 1e8)
    _interval("relative_roughness", relative_roughness, 0, 0.05)
    return 0.25 / (np.log10(relative_roughness / 3.7 + 5.74 / Re**0.9) ** 2)


def colebrook_friction_factor(Re, relative_roughness):
    """Colebrook friction factor.

    Parameters
    ----------
    Re : 1
        Input in the units stated above.
    relative_roughness : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Re=Re, relative_roughness=relative_roughness)
    _positive(Re=Re)
    _interval("relative_roughness", relative_roughness, 0, 0.05)
    if 2300 <= Re < 4000:
        raise ValueError("transition regime requires a separate model")
    if Re < 2300:
        return 64 / Re

    def fun(f):
        return 1 / np.sqrt(f) + 2 * np.log10(relative_roughness / 3.7 + 2.51 / (Re * np.sqrt(f)))

    return optimize.brentq(fun, 0.005, 0.2)


def dynamic_pressure(rho, v):
    """Dynamic pressure.

    Parameters
    ----------
    rho : kg/m^3
        Input in the units stated above.
    v : m/s
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(rho=rho, v=v)
    _positive(rho=rho)
    return 0.5 * rho * v * v


def pitot_velocity(delta_p, rho):
    """Pitot velocity.

    Parameters
    ----------
    delta_p : Pa
        Input in the units stated above.
    rho : kg/m^3
        Input in the units stated above.

    Returns
    -------
    result
        m/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(delta_p=delta_p, rho=rho)
    _nonnegative(delta_p=delta_p)
    _positive(rho=rho)
    return np.sqrt(2 * delta_p / rho)


def venturi_flow(delta_p, rho, A1, A2):
    """Venturi flow.

    Parameters
    ----------
    delta_p : Pa
        Input in the units stated above.
    rho : kg/m^3
        Input in the units stated above.
    A1 : m^2
        Input in the units stated above.
    A2 : m^2
        Input in the units stated above.

    Returns
    -------
    result
        m^3/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(delta_p=delta_p, rho=rho, A1=A1, A2=A2)
    _nonnegative(delta_p=delta_p)
    _positive(rho=rho, A1=A1, A2=A2)
    if A2 >= A1:
        raise ValueError("throat A2 must be smaller than A1")
    return A1 * A2 * np.sqrt(2 * delta_p / (rho * (A1 * A1 - A2 * A2)))


def cavitation_number(p, p_vapor, rho, v):
    """Cavitation number.

    Parameters
    ----------
    p : Pa
        Input in the units stated above.
    p_vapor : Pa
        Input in the units stated above.
    rho : kg/m^3
        Input in the units stated above.
    v : m/s
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(p=p, p_vapor=p_vapor, rho=rho, v=v)
    _positive(rho=rho)
    if v == 0:
        raise ValueError("velocity must be nonzero")
    return (p - p_vapor) / (0.5 * rho * v * v)


def pump_hydraulic_power(rho, Q, head, g=G):
    """Pump hydraulic power.

    Parameters
    ----------
    rho : kg/m^3
        Input in the units stated above.
    Q : m^3/s
        Input in the units stated above.
    head : m
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(rho=rho, Q=Q, head=head, g=g)
    _positive(rho=rho, g=g)
    _nonnegative(Q=Q, head=head)
    return rho * g * Q * head


def pump_input_power(rho, Q, head, efficiency, g=G):
    """Pump input power.

    Parameters
    ----------
    rho : kg/m^3
        Input in the units stated above.
    Q : m^3/s
        Input in the units stated above.
    head : m
        Input in the units stated above.
    efficiency : 1
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        W.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(rho=rho, Q=Q, head=head, efficiency=efficiency, g=g)
    _interval("efficiency", efficiency, 0, 1)
    if efficiency == 0:
        raise ValueError("efficiency must be positive")
    return pump_hydraulic_power(rho, Q, head, g) / efficiency


def mach_number(v, speed_of_sound):
    """Mach number.

    Parameters
    ----------
    v : m/s
        Input in the units stated above.
    speed_of_sound : m/s
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(v=v, speed_of_sound=speed_of_sound)
    _positive(speed_of_sound=speed_of_sound)
    return v / speed_of_sound


def continuity_velocity(A1, v1, A2):
    """Continuity velocity.

    Parameters
    ----------
    A1 : m^2
        Input in the units stated above.
    v1 : m/s
        Input in the units stated above.
    A2 : m^2
        Input in the units stated above.

    Returns
    -------
    result
        m/s.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(A1=A1, v1=v1, A2=A2)
    _positive(A1=A1, A2=A2)
    return A1 * v1 / A2


def hydraulic_diameter(area, wetted_perimeter):
    """Hydraulic diameter.

    Parameters
    ----------
    area : m^2
        Input in the units stated above.
    wetted_perimeter : m
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(area=area, wetted_perimeter=wetted_perimeter)
    _positive(area=area, wetted_perimeter=wetted_perimeter)
    return 4 * area / wetted_perimeter


def pressure_head(p, rho, g=G):
    """Pressure head.

    Parameters
    ----------
    p : Pa
        Input in the units stated above.
    rho : kg/m^3
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(p=p, rho=rho, g=g)
    _positive(rho=rho, g=g)
    return p / (rho * g)


def velocity_head(v, g=G):
    """Velocity head.

    Parameters
    ----------
    v : m/s
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(v=v, g=g)
    _positive(g=g)
    return v * v / (2 * g)


def froude_number(v, length, g=G):
    """Froude number.

    Parameters
    ----------
    v : m/s
        Input in the units stated above.
    length : m
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(v=v, length=length, g=g)
    _positive(length=length, g=g)
    return v / np.sqrt(g * length)


def weber_number(rho, v, length, surface_tension):
    """Weber number.

    Parameters
    ----------
    rho : kg/m^3
        Input in the units stated above.
    v : m/s
        Input in the units stated above.
    length : m
        Input in the units stated above.
    surface_tension : N/m
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(rho=rho, v=v, length=length, surface_tension=surface_tension)
    _positive(rho=rho, length=length, surface_tension=surface_tension)
    return rho * v * v * length / surface_tension


def euler_number(delta_p, rho, v):
    """Euler number.

    Parameters
    ----------
    delta_p : Pa
        Input in the units stated above.
    rho : kg/m^3
        Input in the units stated above.
    v : m/s
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(delta_p=delta_p, rho=rho, v=v)
    _positive(rho=rho)
    if v == 0:
        raise ValueError("velocity must be nonzero")
    return delta_p / (rho * v * v)


def pump_specific_speed(omega, Q, head, g=G):
    """Pump specific speed.

    Parameters
    ----------
    omega : rad/s
        Input in the units stated above.
    Q : m^3/s
        Input in the units stated above.
    head : m
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(omega=omega, Q=Q, head=head, g=g)
    _positive(head=head, g=g)
    _nonnegative(Q=Q, omega=omega)
    return omega * np.sqrt(Q) / (g * head) ** 0.75


def minor_head_loss(K, v, g=G):
    """Minor head loss.

    Parameters
    ----------
    K : 1
        Input in the units stated above.
    v : m/s
        Input in the units stated above.
    g : m/s^2
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(K=K, v=v, g=g)
    _nonnegative(K=K)
    _positive(g=g)
    return K * v * v / (2 * g)


def water_hammer_pressure(rho, wave_speed, delta_v):
    """Water hammer pressure.

    Parameters
    ----------
    rho : kg/m^3
        Input in the units stated above.
    wave_speed : m/s
        Input in the units stated above.
    delta_v : m/s
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/fluid_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(rho=rho, wave_speed=wave_speed, delta_v=delta_v)
    _positive(rho=rho, wave_speed=wave_speed)
    return rho * wave_speed * delta_v

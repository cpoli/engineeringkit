r"""Mechanisms: the slider-crank, the four-bar linkage and serial robot arms.

Reuleaux reduced machines to kinematic chains of rigid links and pairs.
Grashof's condition says which four-bar chains have a link that turns fully.
Freudenstein wrote the four-bar's loop closure as an equation linear in
three coefficients, so a linkage can be synthesised from three input-output
angle pairs by solving a 3 by 3 system. Differentiating the loop-closure
equation once and twice gives linear equations for the angular velocities
and accelerations; Aronhold and Kennedy read the velocity ratios off the
instant centres, which lie three by three on straight lines. Denavit and Hartenberg described any
lower-pair chain with four parameters per joint and a 4 by 4 matrix.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive
from engineeringkit._validation import vector as _vector
from engineeringkit.mechanical_engineering.core.base import FourBarMotion, SliderCrankMotion

__all__ = [
    "slider_crank_displacement",
    "slider_crank_kinematics",
    "fourbar_grashof",
    "fourbar_output_angle",
    "fourbar_kinematics",
    "freudenstein_coefficients",
    "freudenstein_link_lengths",
    "dh_transform",
    "forward_kinematics",
]


def slider_crank_displacement(r, l, theta):
    r"""Slider position of an in-line slider-crank, measured from the crank axis.

    .. math:: x = r\cos\theta + \sqrt{l^2 - r^2\sin^2\theta}

    Parameters
    ----------
    r : float
        Crank radius [m], positive.
    l : float
        Connecting-rod length [m], at least ``r``.
    theta : float or array_like
        Crank angle from top dead centre [rad].

    Returns
    -------
    float or ndarray
        Slider position [m], from :math:`l + r` at :math:`\theta = 0` to
        :math:`l - r` at :math:`\theta = \pi`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    F. Reuleaux, *Theoretische Kinematik* (Vieweg, Braunschweig, 1875).

    Examples
    --------
    >>> float(slider_crank_displacement(1.0, 3.0, 0.0))
    4.0
    """
    _finite(r=r, l=l, theta=theta)
    _positive(r=r, l=l)
    if np.any(np.asarray(l) < r):
        raise ValueError("connecting rod l must be >= crank radius r")
    return r * np.cos(theta) + np.sqrt(l * l - r * r * np.sin(theta) ** 2)


def slider_crank_kinematics(r, l, theta, omega, alpha=0.0) -> SliderCrankMotion:
    r"""Exact slider velocity and acceleration of an in-line slider-crank.

    With the rod at angle :math:`\beta` to the line of stroke,
    :math:`l\sin\beta = r\sin\theta` and :math:`x = r\cos\theta + l\cos\beta`.
    Differentiating,

    .. math::

        \dot\beta = \frac{r\omega\cos\theta}{l\cos\beta}, \qquad
        \dot x = -r\omega\sin\theta - l\dot\beta\sin\beta,

    and once more for :math:`\ddot\beta` and :math:`\ddot x`. For a rod long
    compared with the crank,
    :math:`\ddot x \approx -r\omega^2(\cos\theta + (r/l)\cos 2\theta)`.

    Parameters
    ----------
    r : float
        Crank radius [m], positive.
    l : float
        Connecting-rod length [m], greater than ``r``.
    theta : float or array_like
        Crank angle from top dead centre [rad].
    omega : float
        Crank angular velocity [rad/s].
    alpha : float, optional
        Crank angular acceleration [rad/s^2].

    Returns
    -------
    SliderCrankMotion
        Slider position [m], velocity [m/s] and acceleration [m/s^2], and the
        rod's angle [rad], angular velocity [rad/s] and acceleration
        [rad/s^2].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> motion = slider_crank_kinematics(0.05, 0.2, 0.0, 100.0)
    >>> round(float(motion.acceleration), 6)
    -625.0
    """
    _finite(theta=theta, omega=omega, alpha=alpha)
    _positive(r=r, l=l)
    if l <= r:
        raise ValueError("connecting rod l must be longer than the crank radius r")
    t = np.asarray(theta, dtype=float)
    beta = np.arcsin(r * np.sin(t) / l)
    cb, sb = np.cos(beta), np.sin(beta)
    beta_dot = r * omega * np.cos(t) / (l * cb)
    beta_ddot = (r * alpha * np.cos(t) - r * omega**2 * np.sin(t) + l * beta_dot**2 * sb) / (l * cb)
    return SliderCrankMotion(
        theta=t,
        position=r * np.cos(t) + l * cb,
        velocity=-r * omega * np.sin(t) - l * beta_dot * sb,
        acceleration=-r * alpha * np.sin(t) - r * omega**2 * np.cos(t) - l * beta_ddot * sb - l * beta_dot**2 * cb,
        rod_angle=beta,
        rod_angular_velocity=beta_dot,
        rod_angular_acceleration=beta_ddot,
    )


def fourbar_grashof(lengths):
    r"""Grashof's condition :math:`s + l \le p + q` for a four-bar chain.

    Parameters
    ----------
    lengths : array_like
        The four link lengths [m], positive, in any order.

    Returns
    -------
    bool
        True if the shortest link :math:`s` can rotate fully relative to the
        others. Which link is grounded then decides crank-rocker, double
        crank or double rocker. If False, no link rotates fully.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    F. Grashof, *Theoretische Maschinenlehre*, vol. 2 (Voss, Leipzig, 1883).

    Examples
    --------
    >>> fourbar_grashof([1.0, 2.0, 2.0, 3.0])
    True
    """
    _positive(lengths=lengths)
    if np.asarray(lengths).ndim != 1 or len(lengths) != 4:
        raise ValueError("four positive link lengths are required")
    a = np.sort(np.asarray(lengths, dtype=float))
    return bool(a[0] + a[3] <= a[1] + a[2])


def fourbar_output_angle(a, b, c, d, theta2, crossed=False):
    r"""Output-rocker angle of a four-bar linkage from its loop-closure equation.

    The ground link :math:`d` lies along the x axis from the input pivot to
    the output pivot. The input crank :math:`a` at angle :math:`\theta_2` and
    the output link :math:`c` at angle :math:`\theta_4` are joined by the
    coupler :math:`b`:

    .. math:: a e^{i\theta_2} + b e^{i\theta_3} = d + c e^{i\theta_4}

    Parameters
    ----------
    a, b, c, d : float
        Input, coupler, output and ground link lengths [m], positive.
    theta2 : float or array_like
        Input crank angle [rad].
    crossed : bool, optional
        False for the open circuit, True for the crossed circuit.

    Returns
    -------
    float or ndarray
        :math:`\theta_4` [rad] in (-pi, pi].

    Raises
    ------
    ValueError
        If the linkage cannot be assembled at some ``theta2``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> round(float(np.degrees(fourbar_output_angle(1.0, 3.0, 2.0, 3.0, np.pi / 2))), 4)
    94.8488
    """
    _finite(theta2=theta2)
    _positive(a=a, b=b, c=c, d=d)
    t = np.asarray(theta2, dtype=float)
    A = 2 * c * (d - a * np.cos(t))
    B = -2 * a * c * np.sin(t)
    C = b * b - a * a - c * c - d * d + 2 * a * d * np.cos(t)
    ratio = C / np.hypot(A, B)
    if np.any(np.abs(ratio) > 1 + 1e-12):
        raise ValueError("the linkage cannot be assembled at every input angle")
    theta4 = np.arctan2(B, A) + (-1 if crossed else 1) * np.arccos(np.clip(ratio, -1, 1))
    theta4 = np.angle(np.exp(1j * theta4))
    return float(theta4) if t.ndim == 0 else theta4


def _solve_pair(u, v, w):
    # Real x, y with x u + y v = w for complex u, v, w (two real equations).
    return np.imag(np.conj(v) * w) / np.imag(np.conj(v) * u), np.imag(np.conj(u) * w) / np.imag(np.conj(u) * v)


def fourbar_kinematics(a, b, c, d, theta2, omega2, alpha2=0.0, crossed=False) -> FourBarMotion:
    r"""Angles, angular velocities and angular accelerations of a four-bar linkage.

    Differentiating the loop closure
    :math:`a e^{i\theta_2} + b e^{i\theta_3} = d + c e^{i\theta_4}` gives

    .. math::

        i b\omega_3 e^{i\theta_3} - i c\omega_4 e^{i\theta_4} = -i a\omega_2 e^{i\theta_2}

    and, once more,

    .. math::

        i b\alpha_3 e^{i\theta_3} - i c\alpha_4 e^{i\theta_4} =
        -(i\alpha_2 - \omega_2^2) a e^{i\theta_2} + b\omega_3^2 e^{i\theta_3} - c\omega_4^2 e^{i\theta_4},

    each a pair of real linear equations. The velocity solution is
    :math:`\omega_4/\omega_2 = a\sin(\theta_3 - \theta_2)/(c\sin(\theta_3 - \theta_4))`.

    Parameters
    ----------
    a, b, c, d : float
        Input, coupler, output and ground link lengths [m], positive, as for
        :func:`fourbar_output_angle`.
    theta2 : float or array_like
        Input crank angle [rad].
    omega2 : float
        Input angular velocity [rad/s].
    alpha2 : float, optional
        Input angular acceleration [rad/s^2].
    crossed : bool, optional
        False for the open circuit, True for the crossed circuit.

    Returns
    -------
    FourBarMotion
        :math:`\theta_{2,3,4}` [rad], :math:`\omega_{3,4}` [rad/s] and
        :math:`\alpha_{3,4}` [rad/s^2]; infinite at a change point, where
        coupler and output are in line.

    Raises
    ------
    ValueError
        If the linkage cannot be assembled at some ``theta2``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> motion = fourbar_kinematics(1.0, 3.0, 2.0, 3.0, np.pi / 2, 1.0)
    >>> round(float(motion.omega4), 6)
    0.487298
    """
    _finite(omega2=omega2, alpha2=alpha2)
    theta4 = np.asarray(fourbar_output_angle(a, b, c, d, theta2, crossed=crossed), dtype=float)
    t2 = np.asarray(theta2, dtype=float)
    e2, e4 = np.exp(1j * t2), np.exp(1j * theta4)
    e3 = (d + c * e4 - a * e2) / b
    theta3 = np.angle(e3)
    e3 = np.exp(1j * theta3)
    with np.errstate(divide="ignore", invalid="ignore"):
        omega3, omega4 = _solve_pair(1j * b * e3, -1j * c * e4, -1j * a * omega2 * e2)
        rhs = -(1j * alpha2 - omega2**2) * a * e2 + b * omega3**2 * e3 - c * omega4**2 * e4
        alpha3, alpha4 = _solve_pair(1j * b * e3, -1j * c * e4, rhs)
    return FourBarMotion(theta2=t2, theta3=theta3, theta4=theta4, omega3=omega3, omega4=omega4, alpha3=alpha3, alpha4=alpha4)


def freudenstein_coefficients(input_angles, output_angles):
    r"""Freudenstein's coefficients from prescribed input and output angles.

    Freudenstein's equation

    .. math:: K_1\cos\theta_4 - K_2\cos\theta_2 + K_3 = \cos(\theta_2 - \theta_4)

    with :math:`K_1 = d/a`, :math:`K_2 = d/c` and
    :math:`K_3 = (a^2 - b^2 + c^2 + d^2)/(2ac)` is linear in the :math:`K`.
    Three precision points determine them exactly; more are fitted by least
    squares.

    Parameters
    ----------
    input_angles, output_angles : array_like
        Matching crank angles :math:`\theta_2` and rocker angles
        :math:`\theta_4` [rad], at least three pairs.

    Returns
    -------
    ndarray
        ``[K1, K2, K3]`` [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    References
    ----------
    F. Freudenstein, "An analytical approach to the design of four-link
    mechanisms," *Trans. ASME* 76 (1954) 483-492.

    Examples
    --------
    >>> t2 = np.radians([40.0, 70.0, 100.0])
    >>> t4 = fourbar_output_angle(1.0, 3.0, 2.0, 3.0, t2)
    >>> np.round(freudenstein_coefficients(t2, t4), 6).tolist()
    [3.0, 1.5, 1.25]
    """
    t2 = _vector("input_angles", input_angles, min_length=3)
    t4 = _vector("output_angles", output_angles, min_length=3)
    if t2.shape != t4.shape:
        raise ValueError("one output angle is required per input angle")
    matrix = np.column_stack((np.cos(t4), -np.cos(t2), np.ones_like(t2)))
    if np.linalg.matrix_rank(matrix) < 3:
        raise ValueError("the precision points do not determine the coefficients")
    return np.linalg.lstsq(matrix, np.cos(t2 - t4), rcond=None)[0]


def freudenstein_link_lengths(K, ground=1.0):
    r"""Four-bar link lengths ``(a, b, c, d)`` from Freudenstein's coefficients.

    Parameters
    ----------
    K : array_like
        ``[K1, K2, K3]`` from :func:`freudenstein_coefficients`.
    ground : float, optional
        Ground link length :math:`d` [m], positive; it sets the scale.

    Returns
    -------
    ndarray
        Input, coupler, output and ground lengths [m].

    Raises
    ------
    ValueError
        If ``K1`` or ``K2`` is not positive (a link of negative length) or
        no real coupler exists.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> freudenstein_link_lengths([3.0, 1.5, 1.25], ground=3.0).round(6).tolist()
    [1.0, 3.0, 2.0, 3.0]
    """
    k = _vector("K", K, min_length=3)
    if k.size != 3:
        raise ValueError("K must hold exactly three coefficients")
    _positive(ground=ground)
    if k[0] <= 0 or k[1] <= 0:
        raise ValueError("K1 and K2 must be positive for links of positive length")
    a, c = ground / k[0], ground / k[1]
    b2 = a * a + c * c + ground * ground - 2 * a * c * k[2]
    if b2 <= 0:
        raise ValueError("these coefficients have no real coupler length")
    return np.array([a, np.sqrt(b2), c, float(ground)])


def dh_transform(a, alpha, d, theta):
    r"""Denavit-Hartenberg homogeneous transform from one link frame to the next.

    .. math::

        T = R_z(\theta)\,T_z(d)\,T_x(a)\,R_x(\alpha) =
        \begin{bmatrix}
        \cos\theta & -\sin\theta\cos\alpha & \sin\theta\sin\alpha & a\cos\theta \\
        \sin\theta & \cos\theta\cos\alpha & -\cos\theta\sin\alpha & a\sin\theta \\
        0 & \sin\alpha & \cos\alpha & d \\
        0 & 0 & 0 & 1
        \end{bmatrix}

    Parameters
    ----------
    a : float
        Link length along the new x axis [m].
    alpha : float
        Link twist about the new x axis [rad].
    d : float
        Joint offset along the old z axis [m].
    theta : float
        Joint angle about the old z axis [rad].

    Returns
    -------
    ndarray
        4 by 4 transform; its upper-left 3 by 3 block is a rotation matrix.

    Notes
    -----
    Uses the original (distal) convention. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/mechanical_engineering`.

    References
    ----------
    J. Denavit and R. S. Hartenberg, "A kinematic notation for lower-pair
    mechanisms based on matrices," *J. Appl. Mech.* 22 (1955) 215-221.

    Examples
    --------
    >>> dh_transform(1.0, 0.0, 0.0, np.pi / 2).round(12)[:2, 3].tolist()
    [0.0, 1.0]
    """
    _finite(a=a, alpha=alpha, d=d, theta=theta)
    ct, st, ca, sa = np.cos(theta), np.sin(theta), np.cos(alpha), np.sin(alpha)
    return np.array(
        [
            [ct, -st * ca, st * sa, a * ct],
            [st, ct * ca, -ct * sa, a * st],
            [0.0, sa, ca, d],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )


def forward_kinematics(dh_table):
    r"""Frames of a serial chain from its Denavit-Hartenberg table.

    Parameters
    ----------
    dh_table : array_like
        Shape ``(n, 4)``: one row ``[a, alpha, d, theta]`` per joint
        (lengths [m], angles [rad]).

    Returns
    -------
    ndarray
        Shape ``(n + 1, 4, 4)``: the base frame (identity) followed by the
        cumulative transform :math:`T_1 T_2 \cdots T_k` of each link. The
        last frame's column ``[:3, 3]`` is the end-effector position.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/mechanical_engineering`.

    Examples
    --------
    >>> frames = forward_kinematics([[1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, np.pi / 2]])
    >>> frames[-1][:2, 3].round(12).tolist()
    [1.0, 1.0]
    """
    table = np.asarray(dh_table, dtype=float)
    if table.ndim != 2 or table.shape[1] != 4 or not len(table):
        raise ValueError("dh_table must have shape (n, 4) with n >= 1")
    frames = [np.eye(4)]
    for a, alpha, d, theta in table:
        frames.append(frames[-1] @ dh_transform(a, alpha, d, theta))
    return np.array(frames)

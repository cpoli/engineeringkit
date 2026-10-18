r"""Power systems: AC power, three-phase supply, symmetrical components and load flow.

Tesla's polyphase system delivers constant power through three wires;
Budeanu split AC power into a real part that does work and a reactive part
that only sloshes back and forth; Fortescue decomposed any unbalanced
three-phase set into three balanced ones; and Ward and Hale solved the flow
of power through a whole network by iteration on a digital computer.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.electrical_engineering.core.base import LoadFlowResult

__all__ = [
    "real_power",
    "reactive_power",
    "power_factor",
    "three_phase_power",
    "symmetrical_components",
    "phase_components",
    "admittance_matrix",
    "gauss_seidel_load_flow",
]

_A = np.exp(2j * np.pi / 3)
_FORTESCUE = np.array([[1, 1, 1], [1, _A**2, _A], [1, _A, _A**2]])


def real_power(v_rms, i_rms, pf=1):
    r"""Real (active) power :math:`P = V_{rms} I_{rms} \cos\phi`.

    Parameters
    ----------
    v_rms : float or array_like
        RMS voltage [V], nonnegative.
    i_rms : float or array_like
        RMS current [A], nonnegative.
    pf : float or array_like, optional
        Power factor :math:`\cos\phi` [-], in :math:`[-1, 1]`.

    Returns
    -------
    float or ndarray
        Real power [W].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> real_power(230.0, 10.0, 0.8)
    1840.0
    """
    _finite(v_rms=v_rms, i_rms=i_rms, pf=pf)
    _nonnegative(v_rms=v_rms, i_rms=i_rms)
    _interval("pf", pf, -1, 1)
    return v_rms * i_rms * pf


def reactive_power(v_rms, i_rms, phase_rad):
    r"""Reactive power :math:`Q = V_{rms} I_{rms} \sin\phi`.

    Parameters
    ----------
    v_rms : float or array_like
        RMS voltage [V], nonnegative.
    i_rms : float or array_like
        RMS current [A], nonnegative.
    phase_rad : float or array_like
        Angle :math:`\phi` by which the current lags the voltage [rad].

    Returns
    -------
    float or ndarray
        Reactive power [var]; positive for an inductive (lagging) load.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    C. I. Budeanu, *Puissances réactives et fictives* (Bucharest: Institut
    National Roumain pour l'Étude de l'Aménagement et de l'Utilisation des
    Sources d'Énergie, 1927).

    Examples
    --------
    >>> round(float(reactive_power(230.0, 10.0, np.pi / 6)), 6)
    1150.0
    """
    _finite(v_rms=v_rms, i_rms=i_rms, phase_rad=phase_rad)
    _nonnegative(v_rms=v_rms, i_rms=i_rms)
    return v_rms * i_rms * np.sin(phase_rad)


def power_factor(real_power_value, apparent_power):
    r"""Power factor :math:`P/S`.

    Parameters
    ----------
    real_power_value : float
        Real power [W], with :math:`|P| \le S`.
    apparent_power : float
        Apparent power :math:`S = V_{rms} I_{rms}` [VA], positive.

    Returns
    -------
    float
        Power factor [-], in :math:`[-1, 1]`.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> power_factor(8.0, 10.0)
    0.8
    """
    _finite(real_power_value=real_power_value, apparent_power=apparent_power)
    _positive(apparent_power=apparent_power)
    if abs(real_power_value) > apparent_power:
        raise ValueError("real power magnitude cannot exceed apparent power")
    return real_power_value / apparent_power


def three_phase_power(v_line, i_line, pf=1):
    r"""Balanced three-phase real power :math:`P = \sqrt3 V_L I_L \cos\phi`.

    Parameters
    ----------
    v_line : float or array_like
        RMS line-to-line voltage [V], nonnegative.
    i_line : float or array_like
        RMS line current [A], nonnegative.
    pf : float or array_like, optional
        Power factor [-], in :math:`[-1, 1]`.

    Returns
    -------
    float or ndarray
        Total real power [W] of the three phases.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    N. Tesla, "A new system of alternate current motors and transformers,"
    *Trans. AIEE* 5 (1888) 308-327.

    Examples
    --------
    >>> round(float(three_phase_power(400.0, 10.0)), 1)
    6928.2
    """
    _finite(v_line=v_line, i_line=i_line, pf=pf)
    _nonnegative(v_line=v_line, i_line=i_line)
    _interval("pf", pf, -1, 1)
    return np.sqrt(3) * v_line * i_line * pf


def symmetrical_components(va, vb, vc):
    r"""Fortescue's zero-, positive- and negative-sequence components of three phasors.

    With :math:`a = e^{j2\pi/3}`,

    .. math::

        \begin{bmatrix}V_0\\V_1\\V_2\end{bmatrix} = \frac13
        \begin{bmatrix}1&1&1\\1&a&a^2\\1&a^2&a\end{bmatrix}
        \begin{bmatrix}V_a\\V_b\\V_c\end{bmatrix}.

    Parameters
    ----------
    va, vb, vc : complex or array_like
        Phase phasors [V or A].

    Returns
    -------
    tuple of complex or ndarray
        ``(V0, V1, V2)``, the phase-a members of the zero-, positive- and
        negative-sequence sets, same units as the input.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    C. L. Fortescue, "Method of symmetrical co-ordinates applied to the
    solution of polyphase networks," *Trans. AIEE* 37 (1918) 1027-1140.

    Examples
    --------
    >>> a = np.exp(2j * np.pi / 3)
    >>> [round(abs(complex(x)), 12) for x in symmetrical_components(1, a**2, a)]
    [0.0, 1.0, 0.0]
    """
    _finite(va=va, vb=vb, vc=vc)
    abc = np.array(np.broadcast_arrays(*(np.asarray(v, dtype=complex) for v in (va, vb, vc))))
    v0, v1, v2 = np.tensordot(np.linalg.inv(_FORTESCUE), abc, axes=1)
    if v0.ndim == 0:
        return complex(v0), complex(v1), complex(v2)
    return v0, v1, v2


def phase_components(v0, v1, v2):
    r"""Phase phasors :math:`(V_a, V_b, V_c)` from symmetrical components; the inverse of :func:`symmetrical_components`.

    Parameters
    ----------
    v0, v1, v2 : complex or array_like
        Zero-, positive- and negative-sequence components [V or A].

    Returns
    -------
    tuple of complex or ndarray
        ``(Va, Vb, Vc)``, same units as the input.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> [round(abs(x), 12) for x in phase_components(0, 1, 0)]
    [1.0, 1.0, 1.0]
    """
    _finite(v0=v0, v1=v1, v2=v2)
    seq = np.array(np.broadcast_arrays(*(np.asarray(v, dtype=complex) for v in (v0, v1, v2))))
    va, vb, vc = np.tensordot(_FORTESCUE, seq, axes=1)
    if va.ndim == 0:
        return complex(va), complex(vb), complex(vc)
    return va, vb, vc


def admittance_matrix(n_buses, lines):
    r"""Bus admittance matrix :math:`Y_{bus}` of a network of :math:`\pi`-model lines.

    Each line adds its series admittance :math:`y = 1/z` to the two diagonal
    entries and subtracts it from the two off-diagonal ones, and puts half
    its total shunt admittance on each end.

    Parameters
    ----------
    n_buses : int
        Number of buses, at least two.
    lines : sequence of (int, int, complex) or (int, int, complex, complex)
        ``(i, j, z)`` or ``(i, j, z, y_shunt)``: buses ``i`` and ``j``
        (``0..n_buses-1``), series impedance :math:`z` (nonzero) and total
        line-charging admittance :math:`y_{shunt}`, all per unit.

    Returns
    -------
    ndarray of complex
        :math:`Y_{bus}`, shape ``(n_buses, n_buses)``, per unit.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> admittance_matrix(2, [(0, 1, 0.1j)]).tolist()
    [[-10j, 10j], [10j, -10j]]
    """
    _integer("n_buses", n_buses, minimum=2)
    lines = list(lines)
    if not lines:
        raise ValueError("lines must contain at least one line")
    Y = np.zeros((n_buses, n_buses), dtype=complex)
    for line in lines:
        if len(line) not in (3, 4):
            raise ValueError("each line must be (i, j, z) or (i, j, z, y_shunt)")
        i, j, z = line[:3]
        y_shunt = line[3] if len(line) == 4 else 0.0
        for name, bus in (("i", i), ("j", j)):
            _integer(f"line bus {name}", bus)
            if bus >= n_buses:
                raise ValueError(f"line bus {name} must be in 0..{n_buses - 1}")
        if i == j:
            raise ValueError("a line must join two different buses")
        _finite(z=z, y_shunt=y_shunt)
        if z == 0:
            raise ValueError("line impedance must be nonzero")
        y = 1 / complex(z)
        Y[i, i] += y + y_shunt / 2
        Y[j, j] += y + y_shunt / 2
        Y[i, j] -= y
        Y[j, i] -= y
    return Y


def gauss_seidel_load_flow(Y, p_load, q_load, v_slack=1.0, tol=1e-8, max_iter=500, acceleration=1.0):
    r"""Gauss-Seidel load flow for a slack bus and PQ buses.

    Bus 0 is the slack bus with fixed voltage. At every other bus the net
    injected complex power :math:`S_i = P_i + jQ_i` is specified, and each
    sweep updates the bus voltages in turn,

    .. math::

        V_i \leftarrow \frac{1}{Y_{ii}}\left(\frac{S_i^*}{V_i^*}
        - \sum_{k \ne i} Y_{ik} V_k\right),

    using the newest value of every :math:`V_k`.

    Parameters
    ----------
    Y : array_like of complex
        Bus admittance matrix [pu], square, at least 2x2, nonzero diagonal;
        see :func:`admittance_matrix`.
    p_load, q_load : array_like
        Real [pu] and reactive [pu] power *consumed* at buses
        ``1..n-1`` (length ``n - 1``); negative values are generation.
    v_slack : complex, optional
        Slack-bus voltage [pu], nonzero.
    tol : float, optional
        Stop when the largest voltage update is below this [pu], positive.
    max_iter : int, optional
        Maximum number of sweeps, at least one.
    acceleration : float, optional
        Over-relaxation factor [-] in :math:`(0, 2)`; 1.4 to 1.6 is typical,
        1 is plain Gauss-Seidel.

    Returns
    -------
    LoadFlowResult
        Bus voltages, injected powers (bus 0 gives the slack power), sweep
        count, a convergence flag and the update history. No exception is
        raised on non-convergence; check ``converged``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. PV (voltage-controlled) buses are
    not modelled.

    References
    ----------
    J. B. Ward and H. W. Hale, "Digital computer solution of power-flow
    problems," *Trans. AIEE Part III* 75 (1956) 398-404.

    Examples
    --------
    >>> Y = admittance_matrix(2, [(0, 1, 0.1j)])
    >>> r = gauss_seidel_load_flow(Y, [0.5], [0.0])
    >>> r.converged, round(float(abs(r.voltages[1])), 4)
    (True, 0.9987)
    """
    Y = np.asarray(Y, dtype=complex)
    if Y.ndim != 2 or Y.shape[0] != Y.shape[1] or Y.shape[0] < 2:
        raise ValueError("Y must be a square matrix with at least two buses")
    _finite(Y=Y)
    n = Y.shape[0]
    if np.any(np.diag(Y) == 0):
        raise ValueError("Y must have a nonzero diagonal")
    p, q = np.asarray(p_load, dtype=float), np.asarray(q_load, dtype=float)
    if p.shape != (n - 1,) or q.shape != (n - 1,):
        raise ValueError("p_load and q_load must each have one entry per non-slack bus")
    _finite(p_load=p, q_load=q, v_slack=v_slack)
    if v_slack == 0:
        raise ValueError("v_slack must be nonzero")
    _positive(tol=tol)
    _integer("max_iter", max_iter, minimum=1)
    _interval("acceleration", acceleration, 0, 2, closed=False)

    s_spec = -(p + 1j * q)
    V = np.ones(n, dtype=complex)
    V[0] = v_slack
    history = []
    converged = False
    for _ in range(max_iter):
        largest = 0.0
        for i in range(1, n):
            new = ((s_spec[i - 1] / V[i]).conjugate() - Y[i] @ V + Y[i, i] * V[i]) / Y[i, i]
            new = V[i] + acceleration * (new - V[i])
            largest = max(largest, abs(new - V[i]))
            V[i] = new
        history.append(largest)
        if largest < tol:
            converged = True
            break
    power = V * (Y @ V).conjugate()
    return LoadFlowResult(V, power, len(history), converged, np.array(history))

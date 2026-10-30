r"""Direct-current circuits: Ohm's law, equivalent resistance, bridges and nodal analysis.

Ohm made current proportional to voltage, Wheatstone balanced two ratios of
resistances to measure an unknown one, and Kirchhoff's current and voltage
laws turned any resistive network into a linear system of equations, written
either for the node voltages or for the mesh currents of a planar circuit.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.electrical_engineering.core.base import MeshResult, NodalResult

__all__ = [
    "ohms_law",
    "series_resistance",
    "parallel_resistance",
    "capacitors_parallel",
    "capacitors_series",
    "wheatstone_balanced",
    "nodal_analysis",
    "mesh_analysis",
]


def ohms_law(voltage=None, current=None, resistance=None):
    r"""Solve :math:`V = IR` for whichever of the three quantities is omitted.

    Parameters
    ----------
    voltage : float, optional
        Voltage across the resistor [V].
    current : float, optional
        Current through the resistor [A].
    resistance : float, optional
        Resistance [ohm], positive when given.

    Returns
    -------
    float
        The missing voltage [V], current [A] or resistance [ohm].

    Raises
    ------
    ValueError
        Unless exactly two quantities are given, or if a resistance would be
        inferred from zero current or comes out negative.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    G. S. Ohm, *Die galvanische Kette, mathematisch bearbeitet* (Berlin:
    T. H. Riemann, 1827).

    Examples
    --------
    >>> ohms_law(voltage=12.0, resistance=4.0)
    3.0
    """
    vals = [voltage, current, resistance]
    if sum(v is None for v in vals) != 1:
        raise ValueError("provide exactly two quantities")
    for name, value in zip(("voltage", "current", "resistance"), vals):
        if value is not None:
            _finite(**{name: value})
    if resistance is not None:
        _positive(resistance=resistance)
    if resistance is None and current == 0:
        raise ValueError("nonzero current is needed to infer resistance")
    if voltage is None:
        return current * resistance
    if current is None:
        return voltage / resistance
    result = voltage / current
    _nonnegative(resistance=result)
    return result


def series_resistance(*r):
    r"""Equivalent resistance of resistors in series, :math:`R = \sum R_k`.

    Parameters
    ----------
    *r : float
        Resistances [ohm], nonnegative.

    Returns
    -------
    float
        Equivalent resistance [ohm].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> series_resistance(100.0, 220.0)
    320.0
    """
    _nonnegative(r=r)
    return sum(r)


def parallel_resistance(*r):
    r"""Equivalent resistance of resistors in parallel, :math:`1/R = \sum 1/R_k`.

    Parameters
    ----------
    *r : float
        Resistances [ohm], positive (model a short circuit separately).

    Returns
    -------
    float
        Equivalent resistance [ohm].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> parallel_resistance(100.0, 100.0)
    50.0
    """
    _positive(r=r)
    return 1 / sum(1 / x for x in r)


def capacitors_parallel(*c):
    r"""Equivalent capacitance of capacitors in parallel, :math:`C = \sum C_k`.

    Parameters
    ----------
    *c : float
        Capacitances [F], positive.

    Returns
    -------
    float
        Equivalent capacitance [F].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.
    """
    _positive(c=c)
    return sum(c)


def capacitors_series(*c):
    r"""Equivalent capacitance of capacitors in series, :math:`1/C = \sum 1/C_k`.

    Parameters
    ----------
    *c : float
        Capacitances [F], positive.

    Returns
    -------
    float
        Equivalent capacitance [F].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.
    """
    _positive(c=c)
    return 1 / sum(1 / x for x in c)


def wheatstone_balanced(r1, r2, r3, r4, rtol=1e-9):
    r"""Whether a Wheatstone bridge is balanced, :math:`R_1/R_2 = R_3/R_4`.

    The bridge has :math:`R_1, R_2` in one arm and :math:`R_3, R_4` in the
    other, both arms across the same source. At balance the two midpoints are
    at the same potential and the galvanometer between them reads zero.

    Parameters
    ----------
    r1, r2, r3, r4 : float
        Bridge resistances [ohm], positive.
    rtol : float, optional
        Relative tolerance on the ratio comparison [-], nonnegative.

    Returns
    -------
    bool
        ``True`` at balance.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    C. Wheatstone, "An account of several new instruments and processes for
    determining the constants of a voltaic circuit," *Phil. Trans. R. Soc.*
    133 (1843) 303-327.

    Examples
    --------
    >>> wheatstone_balanced(100.0, 200.0, 50.0, 100.0)
    True
    """
    _finite(r1=r1, r2=r2, r3=r3, r4=r4, rtol=rtol)
    _positive(r1=r1, r2=r2, r3=r3, r4=r4)
    _nonnegative(rtol=rtol)
    balanced = np.isclose(np.divide(r1, r2), np.divide(r3, r4), rtol=rtol)
    return bool(balanced) if balanced.ndim == 0 else balanced


def _check_node(name, node, n_nodes):
    _integer(name, node)
    if node > n_nodes:
        raise ValueError(f"{name} must be a node index in 0..{n_nodes}")


def nodal_analysis(n_nodes, resistors, current_sources=(), voltage_sources=()):
    r"""Solve a resistive DC network by modified nodal analysis.

    Kirchhoff's current law at each non-ground node, with Ohm's law for
    every resistor, gives the conductance system :math:`G\mathbf{v} =
    \mathbf{i}`. Each ideal voltage source adds one unknown (its current)
    and one equation (its terminal voltage), the *modified* nodal system.

    Parameters
    ----------
    n_nodes : int
        Number of non-ground nodes, at least one. Nodes are numbered
        ``1..n_nodes``; node ``0`` is ground.
    resistors : sequence of (int, int, float)
        ``(a, b, R)`` for each resistor between nodes ``a`` and ``b``, with
        :math:`R` [ohm] positive.
    current_sources : sequence of (int, int, float), optional
        ``(a, b, I)``: a source driving current :math:`I` [A] from node
        ``a`` through itself into node ``b``.
    voltage_sources : sequence of (int, int, float), optional
        ``(plus, minus, V)``: an ideal source holding
        :math:`v_{plus} - v_{minus} = V` [V].

    Returns
    -------
    NodalResult
        Node voltages [V] (ground first), resistor currents from ``a`` to
        ``b`` [A], voltage-source currents out of the ``+`` terminal [A] and
        the KCL residual [A].

    Raises
    ------
    ValueError
        For bad node indices or values, or if the network has no unique
        solution (a floating node or a loop of voltage sources).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    G. Kirchhoff, "Ueber den Durchgang eines elektrischen Stromes durch eine
    Ebene, insbesondere durch eine kreisförmige," *Ann. Phys.* 64 (1845)
    497-514; C.-W. Ho, A. E. Ruehli and P. A. Brennan, "The modified nodal
    approach to network analysis," *IEEE Trans. Circuits Syst.* 22 (1975)
    504-509.

    Examples
    --------
    >>> r = nodal_analysis(2, [(1, 2, 1e3), (2, 0, 1e3)], voltage_sources=[(1, 0, 10.0)])
    >>> [round(float(v), 6) for v in r.voltages]
    [0.0, 10.0, 5.0]
    """
    _integer("n_nodes", n_nodes, minimum=1)
    resistors, current_sources, voltage_sources = list(resistors), list(current_sources), list(voltage_sources)
    if not resistors:
        raise ValueError("resistors must contain at least one resistor")
    for group, kind in ((resistors, "resistor"), (current_sources, "current source"), (voltage_sources, "voltage source")):
        for element in group:
            if len(element) != 3:
                raise ValueError(f"each {kind} must be a (node, node, value) triple")
            _check_node(f"{kind} node", element[0], n_nodes)
            _check_node(f"{kind} node", element[1], n_nodes)
            _finite(**{f"{kind} value": element[2]})
            if element[0] == element[1]:
                raise ValueError(f"a {kind} must join two different nodes")
    _positive(resistance=[r for *_, r in resistors])

    m = len(voltage_sources)
    size = n_nodes + m
    A = np.zeros((size, size))
    z = np.zeros(size)

    def stamp(row, col, value):
        if row and col:
            A[row - 1, col - 1] += value

    for a, b, R in resistors:
        g = 1.0 / R
        stamp(a, a, g)
        stamp(b, b, g)
        stamp(a, b, -g)
        stamp(b, a, -g)
    for a, b, current in current_sources:
        if a:
            z[a - 1] -= current
        if b:
            z[b - 1] += current
    for k, (plus, minus, V) in enumerate(voltage_sources):
        row = n_nodes + k
        for node, sign in ((plus, 1.0), (minus, -1.0)):
            if node:
                A[node - 1, row] += sign
                A[row, node - 1] += sign
        z[row] = V
    if np.linalg.cond(A) > 1e12:
        raise ValueError("the network has no unique solution (floating node or voltage-source loop)")
    x = np.linalg.solve(A, z)

    voltages = np.concatenate(([0.0], x[:n_nodes]))
    source_currents = -x[n_nodes:]
    resistor_currents = np.array([(voltages[a] - voltages[b]) / R for a, b, R in resistors])
    leaving = np.zeros(n_nodes + 1)
    for (a, b, _), i in zip(resistors, resistor_currents):
        leaving[a] += i
        leaving[b] -= i
    for a, b, current in current_sources:
        leaving[a] += current
        leaving[b] -= current
    for (plus, minus, _), i in zip(voltage_sources, source_currents):
        leaving[plus] -= i
        leaving[minus] += i
    return NodalResult(voltages, resistor_currents, source_currents, float(np.max(np.abs(leaving[1:]))))


def mesh_analysis(n_meshes, resistors, voltage_sources=()):
    r"""Solve a planar resistive DC network by mesh analysis.

    Each mesh (window) of the planar circuit carries a clockwise mesh
    current. Kirchhoff's voltage law around every mesh, with Ohm's law for
    every resistor, gives the resistance system :math:`R\mathbf{i} =
    \mathbf{v}`: :math:`R_{kk}` is the total resistance around mesh
    :math:`k`, :math:`R_{kl}` minus the resistance mesh :math:`k` shares with
    mesh :math:`l`, and :math:`v_k` the sum of the source EMFs driving
    clockwise around mesh :math:`k`. It is the dual of
    :func:`nodal_analysis`.

    Parameters
    ----------
    n_meshes : int
        Number of meshes, at least one. Meshes are numbered
        ``1..n_meshes``; ``0`` stands for the outside of the circuit.
    resistors : sequence of (int, int, float)
        ``(a, b, R)`` for each resistor on the boundary between meshes
        ``a`` and ``b`` (``b = 0`` for an outer branch), with :math:`R`
        [ohm] positive.
    voltage_sources : sequence of (int, int, float), optional
        ``(a, b, V)``: an ideal source on the boundary between meshes ``a``
        and ``b`` whose EMF :math:`V` [V] drives current clockwise around
        mesh ``a`` (and so anticlockwise around mesh ``b``).

    Returns
    -------
    MeshResult
        Mesh currents [A], resistor currents in the clockwise sense of their
        first mesh [A] and the KVL residual [V].

    Raises
    ------
    ValueError
        For bad mesh indices or values, or if the network has no unique
        solution (a mesh with no resistance).

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`. Current sources are not
    accepted: replace one by its Thévenin equivalent, or use
    :func:`nodal_analysis`.

    References
    ----------
    W. H. Hayt, J. E. Kemmerly and S. M. Durbin, *Engineering Circuit
    Analysis*, 8th ed. (New York: McGraw-Hill, 2012), ch. 4.

    Examples
    --------
    >>> r = mesh_analysis(2, [(1, 0, 2.0), (1, 2, 4.0), (2, 0, 4.0)], voltage_sources=[(1, 0, 12.0)])
    >>> [round(float(i), 6) for i in r.mesh_currents]
    [3.0, 1.5]
    """
    _integer("n_meshes", n_meshes, minimum=1)
    resistors, voltage_sources = list(resistors), list(voltage_sources)
    if not resistors:
        raise ValueError("resistors must contain at least one resistor")
    for group, kind in ((resistors, "resistor"), (voltage_sources, "voltage source")):
        for element in group:
            if len(element) != 3:
                raise ValueError(f"each {kind} must be a (mesh, mesh, value) triple")
            _check_node(f"{kind} mesh", element[0], n_meshes)
            _check_node(f"{kind} mesh", element[1], n_meshes)
            _finite(**{f"{kind} value": element[2]})
            if element[0] == element[1]:
                raise ValueError(f"a {kind} must lie between two different meshes")
    _positive(resistance=[r for *_, r in resistors])

    A = np.zeros((n_meshes, n_meshes))
    v = np.zeros(n_meshes)
    for a, b, R in resistors:
        for row, col, value in ((a, a, R), (b, b, R), (a, b, -R), (b, a, -R)):
            if row and col:
                A[row - 1, col - 1] += value
    for a, b, V in voltage_sources:
        if a:
            v[a - 1] += V
        if b:
            v[b - 1] -= V
    if np.linalg.cond(A) > 1e12:
        raise ValueError("the network has no unique solution (a mesh with no resistance)")
    currents = np.concatenate(([0.0], np.linalg.solve(A, v)))

    resistor_currents = np.array([currents[a] - currents[b] for a, b, _ in resistors])
    loop = np.zeros(n_meshes + 1)
    for (a, b, R), i in zip(resistors, resistor_currents):
        loop[a] -= R * i
        loop[b] += R * i
    for a, b, V in voltage_sources:
        loop[a] += V
        loop[b] -= V
    return MeshResult(currents[1:], resistor_currents, float(np.max(np.abs(loop[1:]))))

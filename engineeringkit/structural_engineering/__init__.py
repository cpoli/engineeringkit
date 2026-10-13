"""Structural engineering: stress, beams, columns, trusses, and energy methods."""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import integer as _integer
from engineeringkit._validation import interval as _interval
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive


def hooke_stress(E: float, strain):
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
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(E=E, strain=strain)
    _positive(E=E)
    return E * np.asarray(strain)


def axial_deformation(force, length, area, E):
    """Axial deformation.

    Parameters
    ----------
    force : N
        Input in the units stated above.
    length : m
        Input in the units stated above.
    area : m^2
        Input in the units stated above.
    E : Pa
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(force=force, length=length, area=area, E=E)
    _positive(length=length, area=area, E=E)
    return force * length / (area * E)


def thermal_strain(alpha, delta_T):
    """Thermal strain.

    Parameters
    ----------
    alpha : 1/K
        Input in the units stated above.
    delta_T : K (temperature difference)
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return alpha * delta_T


def poisson_lateral_strain(axial_strain, nu):
    """Poisson lateral strain.

    Parameters
    ----------
    axial_strain : 1
        Input in the units stated above.
    nu : 1
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(axial_strain=axial_strain, nu=nu)
    _interval("nu", nu, -1, 0.5, closed=False)
    return -nu * axial_strain


def euler_buckling_load(E, I, L, K=1.0):
    """Euler buckling load.

    Parameters
    ----------
    E : Pa
        Input in the units stated above.
    I : m^4
        Input in the units stated above.
    L : m
        Input in the units stated above.
    K : 1
        Input in the units stated above.

    Returns
    -------
    result
        N.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(E=E, I=I, L=L, K=K)
    _positive(E=E, I=I, L=L, K=K)
    return np.pi**2 * E * I / (K * L) ** 2


def radius_of_gyration(I, A):
    """Radius of gyration.

    Parameters
    ----------
    I : m^4
        Input in the units stated above.
    A : m^2
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(I=I, A=A)
    _positive(I=I, A=A)
    return np.sqrt(I / A)


def slenderness_ratio(L, K, r):
    """Slenderness ratio.

    Parameters
    ----------
    L : m
        Input in the units stated above.
    K : 1
        Input in the units stated above.
    r : m
        Input in the units stated above.

    Returns
    -------
    result
        1.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(L=L, K=K, r=r)
    _positive(L=L, K=K, r=r)
    return K * L / r


def torsion_shear_stress(T, r, J):
    """Torsion shear stress.

    Parameters
    ----------
    T : N m
        Input in the units stated above.
    r : m
        Input in the units stated above.
    J : m^4
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T=T, r=r, J=J)
    _nonnegative(r=r)
    _positive(J=J)
    return T * r / J


def torsion_angle(T, L, G, J):
    """Torsion angle.

    Parameters
    ----------
    T : N m
        Input in the units stated above.
    L : m
        Input in the units stated above.
    G : Pa
        Input in the units stated above.
    J : m^4
        Input in the units stated above.

    Returns
    -------
    result
        rad.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(T=T, L=L, G=G, J=J)
    _positive(L=L, G=G, J=J)
    return T * L / (G * J)


def rectangular_second_moment(b, h):
    """Rectangular second moment.

    Parameters
    ----------
    b : m
        Input in the units stated above.
    h : m
        Input in the units stated above.

    Returns
    -------
    result
        m^4.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(b=b, h=h)
    _positive(b=b, h=h)
    return b * h**3 / 12


def parallel_axis(Ic, A, d):
    """Parallel axis.

    Parameters
    ----------
    Ic : m^4
        Input in the units stated above.
    A : m^2
        Input in the units stated above.
    d : m
        Input in the units stated above.

    Returns
    -------
    result
        m^4.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(Ic=Ic, A=A, d=d)
    _nonnegative(Ic=Ic)
    _positive(A=A)
    return Ic + A * d * d


def bending_stress(M, y, I):
    """Bending stress.

    Parameters
    ----------
    M : N m
        Input in the units stated above.
    y : m
        Input in the units stated above.
    I : m^4
        Input in the units stated above.

    Returns
    -------
    result
        Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(M=M, y=y, I=I)
    _positive(I=I)
    return M * y / I


def mohr_circle(sx, sy, txy):
    """Mohr circle.

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
        (center, radius, maximum principal stress, minimum principal stress), all Pa.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    c = (sx + sy) / 2
    r = np.hypot((sx - sy) / 2, txy)
    return c, r, c + r, c - r


def castigliano_deflection(dU_dP):
    """Return a caller-computed strain-energy derivative (not an energy solver).

    Parameters
    ----------
    dU_dP : m
        Input in the units stated above.

    Returns
    -------
    result
        m.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return dU_dP


@dataclass(frozen=True)
class BeamResult:
    x: np.ndarray
    shear: np.ndarray
    moment: np.ndarray
    deflection: np.ndarray
    reactions: tuple[float, float]

    def summary(self):
        return {
            "max_abs_shear": float(np.max(abs(self.shear))),
            "max_abs_moment": float(np.max(abs(self.moment))),
            "max_abs_deflection": float(np.max(abs(self.deflection))),
            "reactions": self.reactions,
        }

    def plot(self):
        fig, axs = plt.subplots(3, 1, sharex=True)
        axs[0].plot(self.x, self.shear)
        axs[1].plot(self.x, self.moment)
        axs[2].plot(self.x, self.deflection)
        axs[2].set_xlabel("x [m]")
        return fig, axs


@dataclass(frozen=True)
class SimplySupportedBeam:
    length: float
    young_modulus: float
    second_moment: float

    def __post_init__(self):
        _positive(length=self.length, young_modulus=self.young_modulus, second_moment=self.second_moment)

    def solve_point_load(self, load: float, position: float | None = None, n: int = 401):
        _finite(load=load)
        _integer("n", n, 2)
        L = self.length
        a = L / 2 if position is None else position
        if not 0 <= a <= L:
            raise ValueError("position must lie on beam")
        b = L - a
        ra = load * b / L
        rb = load * a / L
        x = np.linspace(0, L, n)
        V = np.where(x < a, ra, ra - load)
        M = np.where(x <= a, ra * x, rb * (L - x))
        # exact Euler-Bernoulli deflection, downward negative
        E, I = self.young_modulus, self.second_moment
        y = np.empty_like(x)
        left = x <= a
        y[left] = -load * b * x[left] * (L**2 - b**2 - x[left] ** 2) / (6 * L * E * I)
        xr = L - x[~left]
        y[~left] = -load * a * xr * (L**2 - a**2 - xr**2) / (6 * L * E * I)
        return BeamResult(x, V, M, y, (ra, rb))

    def solve_uniform_load(self, w: float, n: int = 401):
        _finite(w=w)
        _integer("n", n, 2)
        L = self.length
        x = np.linspace(0, L, n)
        R = w * L / 2
        V = R - w * x
        M = R * x - w * x * x / 2
        y = -w * x * (L**3 - 2 * L * x * x + x**3) / (24 * self.young_modulus * self.second_moment)
        return BeamResult(x, V, M, y, (R, R))


def rectangular_section_modulus(b, h):
    """Elastic section modulus ``Z = I/c`` for a rectangle [m^3].

    Parameters
    ----------
    b : m
        Input in the units stated above.
    h : m
        Input in the units stated above.

    Returns
    -------
    result
        m^3.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(b=b, h=h)
    _positive(b=b, h=h)
    return b * h**2 / 6


def strain_energy_axial(force, length, area, E):
    """Linear-elastic axial strain energy ``P^2 L/(2AE)`` [J].

    Parameters
    ----------
    force : N
        Input in the units stated above.
    length : m
        Input in the units stated above.
    area : m^2
        Input in the units stated above.
    E : Pa
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(force=force, length=length, area=area, E=E)
    _positive(length=length, area=area, E=E)
    return force**2 * length / (2 * area * E)


def strain_energy_bending_constant_moment(moment, length, E, I):
    """Bending strain energy for constant moment ``M^2 L/(2EI)`` [J].

    Parameters
    ----------
    moment : N m
        Input in the units stated above.
    length : m
        Input in the units stated above.
    E : Pa
        Input in the units stated above.
    I : m^4
        Input in the units stated above.

    Returns
    -------
    result
        J.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    _finite(moment=moment, length=length, E=E, I=I)
    _positive(length=length, E=E, I=I)
    return moment**2 * length / (2 * E * I)


def maxwell_betti_reciprocity(delta_ij, delta_ji, rtol=1e-9, atol=1e-12):
    """Check Maxwell-Betti reciprocal displacement equality.

    Parameters
    ----------
    delta_ij : m
        Input in the units stated above.
    delta_ji : m
        Input in the units stated above.
    rtol : 1
        Input in the units stated above.
    atol : m (absolute displacement tolerance)
        Input in the units stated above.

    Returns
    -------
    result
        bool.

    Notes
    -----
    Assumptions, governing equations and validity limits are discussed in
    :doc:`/guides/structural_engineering`. SI inputs are required; raw numbers carry no
    dimensional metadata.
    """
    return bool(np.isclose(delta_ij, delta_ji, rtol=rtol, atol=atol))


@dataclass(frozen=True)
class TrussResult:
    displacements: np.ndarray
    reactions: np.ndarray
    member_forces: np.ndarray
    member_stresses: np.ndarray

    def summary(self):
        return {
            "max_abs_displacement": float(np.max(np.abs(self.displacements))),
            "max_abs_member_force": float(np.max(np.abs(self.member_forces))),
        }

    def plot(self, nodes, members, scale=1.0):
        fig, ax = plt.subplots()
        nodes = np.asarray(nodes, float)
        d = self.displacements.reshape(-1, 2)
        deformed = nodes + scale * d
        for i, j in members:
            ax.plot(nodes[[i, j], 0], nodes[[i, j], 1], "--", alpha=0.5)
            ax.plot(deformed[[i, j], 0], deformed[[i, j], 1])
        ax.set_aspect("equal")
        ax.set(xlabel="x [m]", ylabel="y [m]")
        return fig, ax


@dataclass(frozen=True)
class PlanarTruss:
    nodes: np.ndarray
    members: tuple
    young_modulus: float
    area: float

    def solve(self, loads, fixed_dofs):
        """Solve a pin-jointed 2-D truss by the direct stiffness method."""
        xy = np.asarray(self.nodes, float)
        loads = np.asarray(loads, float).reshape(-1)
        if xy.ndim != 2 or xy.shape[1] != 2:
            raise ValueError("nodes must have shape (n,2)")
        if self.young_modulus <= 0 or self.area <= 0:
            raise ValueError("E and area must be positive")
        _finite(nodes=xy, loads=loads)
        _positive(E=self.young_modulus, area=self.area)
        if len(xy) < 2 or not len(self.members):
            raise ValueError("truss requires nodes and members")
        ndof = 2 * len(xy)
        if loads.size != ndof:
            raise ValueError("loads must contain two components per node")
        K = np.zeros((ndof, ndof))
        geom = []
        for i, j in self.members:
            _integer("member node", i, 0)
            _integer("member node", j, 0)
            if i >= len(xy) or j >= len(xy):
                raise ValueError("member node index out of range")
            dx, dy = xy[j] - xy[i]
            L = np.hypot(dx, dy)
            if L <= 0:
                raise ValueError("members must have positive length")
            c, s = dx / L, dy / L
            k = self.young_modulus * self.area / L
            ke = k * np.array(
                [
                    [c * c, c * s, -c * c, -c * s],
                    [c * s, s * s, -c * s, -s * s],
                    [-c * c, -c * s, c * c, c * s],
                    [-c * s, -s * s, c * s, s * s],
                ]
            )
            dofs = [2 * i, 2 * i + 1, 2 * j, 2 * j + 1]
            K[np.ix_(dofs, dofs)] += ke
            geom.append((i, j, L, c, s, dofs))
        raw_fixed = np.asarray(fixed_dofs)
        if raw_fixed.ndim != 1 or not np.issubdtype(raw_fixed.dtype, np.integer):
            raise ValueError("fixed_dofs must be an integer vector")
        if np.any(raw_fixed < 0) or np.any(raw_fixed >= ndof):
            raise ValueError("fixed degree of freedom out of range")
        fixed = np.unique(raw_fixed.astype(int))
        free = np.setdiff1d(np.arange(ndof), fixed)
        u = np.zeros(ndof)
        try:
            u[free] = np.linalg.solve(K[np.ix_(free, free)], loads[free])
        except np.linalg.LinAlgError as exc:
            raise ValueError("truss is unstable or insufficiently constrained") from exc
        reactions = K @ u - loads
        forces = []
        stresses = []
        for i, j, L, c, s, dofs in geom:
            extension = np.dot(np.array([-c, -s, c, s]), u[dofs])
            force = self.young_modulus * self.area / L * extension
            forces.append(force)
            stresses.append(force / self.area)
        return TrussResult(u, reactions, np.asarray(forces), np.asarray(stresses))

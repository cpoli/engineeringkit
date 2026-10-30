r"""Crystal structures: packing, theoretical density and Bragg diffraction.

Metals crystallise mostly as simple cubic, body-centred cubic, face-centred
cubic or hexagonal close-packed arrays of touching spheres. The geometry of
the unit cell fixes the packing factor and, with the molar mass, the density.
X-rays reflect from a family of lattice planes when the path difference
between neighbouring planes is a whole number of wavelengths,
:math:`n\lambda = 2d\sin\theta` (W. L. Bragg, 1913), so a diffraction pattern
measures the lattice parameter.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive

__all__ = [
    "cubic_lattice_parameter",
    "atomic_packing_factor",
    "theoretical_density",
    "interplanar_spacing",
    "bragg_angle",
    "diffraction_peaks",
]

AVOGADRO = 6.02214076e23
"""float: Avogadro constant [1/mol], exact in the 2019 SI."""

# Atoms per cell and edge length per atomic radius of the cubic cells.
_CUBIC = {"sc": (1, 2.0), "bcc": (2, 4 / np.sqrt(3)), "fcc": (4, 2 * np.sqrt(2))}


def _structure(structure, allowed):
    if structure not in allowed:
        raise ValueError(f"structure must be one of {sorted(allowed)}")


def cubic_lattice_parameter(atomic_radius, structure):
    r"""Edge length of a cubic unit cell of touching hard spheres.

    .. math:: a_{sc} = 2R, \quad a_{bcc} = \frac{4R}{\sqrt3}, \quad a_{fcc} = 2\sqrt2\,R

    Parameters
    ----------
    atomic_radius : float or array_like
        Atomic radius :math:`R` [m], positive.
    structure : {"sc", "bcc", "fcc"}
        Cubic crystal structure. Atoms touch along the cube edge, the body
        diagonal and the face diagonal respectively.

    Returns
    -------
    float or ndarray
        Lattice parameter :math:`a` [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> round(float(cubic_lattice_parameter(0.128e-9, "fcc")) * 1e9, 4)
    0.362
    """
    _positive(atomic_radius=atomic_radius)
    _structure(structure, _CUBIC)
    return _CUBIC[structure][1] * np.asarray(atomic_radius)


def atomic_packing_factor(structure):
    r"""Fraction of a unit cell's volume filled by its hard-sphere atoms.

    .. math:: \mathrm{APF} = \frac{n\,\tfrac43\pi R^3}{V_c}

    Parameters
    ----------
    structure : {"sc", "bcc", "fcc", "hcp"}
        Crystal structure; ``"hcp"`` uses the ideal ratio
        :math:`c/a = \sqrt{8/3}`.

    Returns
    -------
    float
        Packing factor [-]: 0.52, 0.68, 0.74 and 0.74.

    Notes
    -----
    FCC and ideal HCP reach Kepler's close-packing limit
    :math:`\pi/\sqrt{18}`. Inputs are raw SI numbers with no dimensional
    metadata; see :doc:`/guides/materials`.

    Examples
    --------
    >>> round(atomic_packing_factor("bcc"), 4)
    0.6802
    """
    _structure(structure, {*_CUBIC, "hcp"})
    sphere = 4 / 3 * np.pi
    if structure == "hcp":
        # Six atoms in a hexagonal prism of side a = 2R and height c = sqrt(8/3) a.
        a = 2.0
        return float(6 * sphere / (3 * np.sqrt(3) / 2 * a * a * np.sqrt(8 / 3) * a))
    n, a = _CUBIC[structure]
    return float(n * sphere / a**3)


def theoretical_density(atoms_per_cell, molar_mass, cell_volume):
    r"""Density of a perfect crystal from its unit cell, :math:`\rho = nA/(V_c N_A)`.

    Parameters
    ----------
    atoms_per_cell : float or array_like
        Atoms per unit cell :math:`n` [-], positive (1, 2, 4 or 6).
    molar_mass : float or array_like
        Molar mass :math:`A` [kg/mol], positive.
    cell_volume : float or array_like
        Unit cell volume :math:`V_c` [m^3], positive (:math:`a^3` for a
        cubic cell).

    Returns
    -------
    float or ndarray
        Density [kg/m^3].

    Notes
    -----
    Vacancies and other defects make the measured density slightly lower.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> a = cubic_lattice_parameter(0.128e-9, "fcc")
    >>> round(float(theoretical_density(4, 63.5e-3, a**3)))
    8888
    """
    _positive(atoms_per_cell=atoms_per_cell, molar_mass=molar_mass, cell_volume=cell_volume)
    return np.asarray(atoms_per_cell) * molar_mass / (np.asarray(cell_volume) * AVOGADRO)


def interplanar_spacing(a, h, k, l):  # noqa: E741
    r"""Spacing of the :math:`(hkl)` planes of a cubic crystal, :math:`d = a/\sqrt{h^2 + k^2 + l^2}`.

    Parameters
    ----------
    a : float or array_like
        Lattice parameter [m], positive.
    h, k, l : int or array_like
        Miller indices, integers not all zero.

    Returns
    -------
    float or ndarray
        Interplanar spacing :math:`d_{hkl}` [m].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> round(float(interplanar_spacing(0.2866e-9, 2, 2, 0)) * 1e9, 4)
    0.1013
    """
    _positive(a=a)
    h, k, l = (np.asarray(i) for i in (h, k, l))  # noqa: E741
    if any(not np.issubdtype(i.dtype, np.integer) for i in (h, k, l)):
        raise ValueError("Miller indices must be integers")
    s = h * h + k * k + l * l
    if np.any(s == 0):
        raise ValueError("Miller indices must not all be zero")
    return np.asarray(a) / np.sqrt(s)


def bragg_angle(d, wavelength, order=1):
    r"""Bragg angle :math:`\theta` from :math:`n\lambda = 2d\sin\theta`.

    Parameters
    ----------
    d : float or array_like
        Interplanar spacing [m], positive.
    wavelength : float or array_like
        X-ray wavelength :math:`\lambda` [m], positive.
    order : int, optional
        Order of reflection :math:`n`, at least 1.

    Returns
    -------
    float or ndarray
        Angle between the beam and the planes [rad]; the diffractometer
        reads :math:`2\theta`.

    Raises
    ------
    ValueError
        If :math:`n\lambda > 2d`, when the planes cannot reflect.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    W. L. Bragg, "The diffraction of short electromagnetic waves by a
    crystal," *Proc. Camb. Phil. Soc.* 17 (1913) 43-57.

    Examples
    --------
    >>> round(float(np.degrees(2 * bragg_angle(0.1013e-9, 0.1790e-9))), 1)
    124.1
    """
    _positive(d=d, wavelength=wavelength)
    _integer("order", order, 1)
    s = order * np.asarray(wavelength) / (2 * np.asarray(d))
    if np.any(s > 1):
        raise ValueError("n * wavelength exceeds 2d, so the planes cannot reflect")
    return np.arcsin(s)


def diffraction_peaks(a, wavelength, structure, max_index=4):
    r"""First-order powder diffraction peaks of a cubic crystal, in order of angle.

    Planes whose atoms scatter out of phase give no peak: in BCC
    :math:`h + k + l` must be even, and in FCC :math:`h, k, l` must be all odd
    or all even.

    Parameters
    ----------
    a : float
        Lattice parameter [m], positive.
    wavelength : float
        X-ray wavelength [m], positive.
    structure : {"sc", "bcc", "fcc"}
        Cubic crystal structure.
    max_index : int, optional
        Largest Miller index searched, at least 1.

    Returns
    -------
    tuple of ndarray
        ``(planes, two_theta)``: one :math:`(hkl)` row with
        :math:`h \ge k \ge l \ge 0` per peak, shape ``(n, 3)``, and the
        diffraction angles :math:`2\theta` [rad], increasing. Families with
        the same :math:`h^2 + k^2 + l^2`, such as (300) and (221), share a
        peak and are listed once.

    Notes
    -----
    Families with :math:`\lambda > 2d` are dropped. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/materials`.

    Examples
    --------
    >>> planes, _ = diffraction_peaks(0.2866e-9, 0.1790e-9, "bcc")
    >>> planes[:3].tolist()
    [[1, 1, 0], [2, 0, 0], [2, 1, 1]]
    """
    _positive(a=a, wavelength=wavelength)
    _structure(structure, _CUBIC)
    _integer("max_index", max_index, 1)
    planes, seen = [], set()
    for h in range(max_index + 1):
        for k in range(h + 1):
            for l in range(k + 1):  # noqa: E741
                s = h * h + k * k + l * l
                if s == 0 or s in seen:
                    continue
                if structure == "bcc" and (h + k + l) % 2:
                    continue
                if structure == "fcc" and len({h % 2, k % 2, l % 2}) > 1:
                    continue
                if wavelength * np.sqrt(s) > 2 * a:
                    continue
                seen.add(s)
                planes.append((h, k, l))
    planes_array = np.array(sorted(planes, key=lambda p: p[0] ** 2 + p[1] ** 2 + p[2] ** 2), dtype=int).reshape(-1, 3)
    d = a / np.sqrt(np.sum(planes_array**2, axis=1))
    return planes_array, 2 * np.arcsin(wavelength / (2 * d))

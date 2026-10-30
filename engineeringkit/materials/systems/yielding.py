r"""Yield criteria: Tresca's maximum shear and von Mises's distortion energy.

Tresca found that metals forced through dies flow when the largest shear
stress reaches a constant, :math:`\sigma_1 - \sigma_3 = \sigma_y`. Von Mises
replaced the hexagonal prism this gives in principal-stress space with the
smooth cylinder through its corners, the second deviatoric invariant
:math:`J_2 = \sigma_y^2/3`. Both are written in principal stresses, the
eigenvalues of the symmetric stress tensor.
"""

from __future__ import annotations

import numpy as np

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import positive as _positive

__all__ = ["principal_stresses", "von_mises_stress", "von_mises_plane_stress", "tresca_equivalent", "factor_of_safety"]


def _stress_components(sx, sy, sz, txy, tyz, tzx):
    _finite(sx=sx, sy=sy, sz=sz, txy=txy, tyz=tyz, tzx=tzx)
    return np.broadcast_arrays(*(np.asarray(c, float) for c in (sx, sy, sz, txy, tyz, tzx)))


def principal_stresses(sx, sy, sz, txy=0.0, tyz=0.0, tzx=0.0):
    r"""Principal stresses of a three-dimensional stress state, largest first.

    The eigenvalues of

    .. math:: \boldsymbol\sigma = \begin{pmatrix} \sigma_x & \tau_{xy} & \tau_{zx} \\
              \tau_{xy} & \sigma_y & \tau_{yz} \\ \tau_{zx} & \tau_{yz} & \sigma_z \end{pmatrix}

    Parameters
    ----------
    sx, sy, sz : float or array_like
        Normal stresses [Pa].
    txy, tyz, tzx : float or array_like, optional
        Shear stresses [Pa].

    Returns
    -------
    ndarray
        :math:`(\sigma_1, \sigma_2, \sigma_3)` [Pa] along the last axis,
        :math:`\sigma_1 \ge \sigma_2 \ge \sigma_3`; pass them to
        :func:`tresca_equivalent`.

    Notes
    -----
    ``numpy.linalg.eigvalsh`` finds the eigenvalues of each symmetric tensor,
    the roots of :math:`\sigma^3 - I_1\sigma^2 + I_2\sigma - I_3 = 0`.
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    Examples
    --------
    >>> principal_stresses(0.0, 0.0, 0.0, txy=50.0).round(9).tolist()
    [50.0, 0.0, -50.0]
    """
    sx, sy, sz, txy, tyz, tzx = _stress_components(sx, sy, sz, txy, tyz, tzx)
    tensor = np.stack([np.stack([sx, txy, tzx], -1), np.stack([txy, sy, tyz], -1), np.stack([tzx, tyz, sz], -1)], -2)
    return np.linalg.eigvalsh(tensor)[..., ::-1]


def von_mises_stress(sx, sy, sz, txy=0.0, tyz=0.0, tzx=0.0):
    r"""Von Mises equivalent stress of a three-dimensional stress state.

    .. math:: \sigma_{vM} = \sqrt{\tfrac12\left[(\sigma_x - \sigma_y)^2 + (\sigma_y - \sigma_z)^2
              + (\sigma_z - \sigma_x)^2\right] + 3(\tau_{xy}^2 + \tau_{yz}^2 + \tau_{zx}^2)} = \sqrt{3J_2}

    Parameters
    ----------
    sx, sy, sz : float or array_like
        Normal stresses [Pa].
    txy, tyz, tzx : float or array_like, optional
        Shear stresses [Pa].

    Returns
    -------
    float or ndarray
        Equivalent stress [Pa], unchanged by a hydrostatic stress.

    Notes
    -----
    With :math:`\sigma_z = \tau_{yz} = \tau_{zx} = 0` it reduces to
    :func:`von_mises_plane_stress`. Inputs are raw SI numbers with no
    dimensional metadata; see :doc:`/guides/materials`.

    Examples
    --------
    >>> float(von_mises_stress(100.0, 100.0, 100.0))
    0.0
    """
    sx, sy, sz, txy, tyz, tzx = _stress_components(sx, sy, sz, txy, tyz, tzx)
    normal = 0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2)
    return np.sqrt(normal + 3 * (txy * txy + tyz * tyz + tzx * tzx))


def von_mises_plane_stress(sx, sy, txy):
    r"""Von Mises equivalent stress for plane stress.

    .. math:: \sigma_{vM} = \sqrt{\sigma_x^2 - \sigma_x\sigma_y + \sigma_y^2 + 3\tau_{xy}^2}

    Parameters
    ----------
    sx, sy : float or array_like
        Normal stresses [Pa].
    txy : float or array_like
        In-plane shear stress [Pa].

    Returns
    -------
    float or ndarray
        Equivalent stress [Pa]; yield is predicted when it reaches the
        uniaxial yield strength.

    Notes
    -----
    Assumes an isotropic, pressure-insensitive ductile metal with
    :math:`\sigma_z = \tau_{xz} = \tau_{yz} = 0`. Inputs are raw SI numbers
    with no dimensional metadata; see :doc:`/guides/materials`.

    References
    ----------
    R. von Mises, "Mechanik der festen Körper im plastisch-deformablen
    Zustand," *Nachr. Ges. Wiss. Göttingen, Math.-Phys. Kl.* (1913) 582-592.

    Examples
    --------
    >>> round(float(von_mises_plane_stress(0.0, 0.0, 1.0)), 6)
    1.732051
    """
    _finite(sx=sx, sy=sy, txy=txy)
    sx, sy, txy = np.asarray(sx), np.asarray(sy), np.asarray(txy)
    return np.sqrt(sx * sx - sx * sy + sy * sy + 3 * txy * txy)


def tresca_equivalent(principal):
    r"""Tresca equivalent stress, the largest minus the smallest principal stress.

    Parameters
    ----------
    principal : array_like
        Principal stresses [Pa] along the last axis, in any order. Include
        the zero out-of-plane stress for plane stress.

    Returns
    -------
    float or ndarray
        :math:`\sigma_{max} - \sigma_{min} = 2\tau_{max}` [Pa], one value per
        stress state.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.

    References
    ----------
    H. Tresca, "Mémoire sur l'écoulement des corps solides soumis à de
    fortes pressions," *C. R. Acad. Sci. Paris* 59 (1864) 754-758.

    Examples
    --------
    >>> float(tresca_equivalent([100.0, 0.0, -50.0]))
    150.0
    """
    p = np.asarray(principal, float)
    _finite(principal=p)
    if p.ndim == 0:
        raise ValueError("principal must list the principal stresses along its last axis")
    return np.ptp(p, axis=-1)


def factor_of_safety(strength, stress):
    r"""Factor of safety :math:`n = S / \sigma`.

    Parameters
    ----------
    strength : float or array_like
        Material strength [Pa], positive.
    stress : float or array_like
        Equivalent working stress [Pa], positive.

    Returns
    -------
    float or ndarray
        Safety factor [-].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/materials`.
    """
    _positive(strength=strength, stress=stress)
    return np.asarray(strength) / stress

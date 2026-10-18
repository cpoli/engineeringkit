"""Result containers for engineeringkit.materials.

Most materials models are closed-form relations that return plain floats or
arrays. A measured stress-strain curve and an integrated crack-growth history
are sampled curves, so they are frozen dataclasses that can be summarised and
plotted consistently.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import integer as _integer
from engineeringkit._validation import positive as _positive

__all__ = ["StressStrainCurve", "CrackGrowthResult"]


@dataclass(frozen=True)
class StressStrainCurve:
    """A sampled engineering stress-strain curve from a tensile test.

    The arrays are validated (matching vectors, strictly increasing strain)
    and copied on construction.
    """

    strain: np.ndarray
    """ndarray: Engineering strain, strictly increasing [-]."""

    stress: np.ndarray
    """ndarray: Engineering stress at each strain [Pa]."""

    def __post_init__(self):
        from engineeringkit.materials.systems.elasticity import toughness

        strain, stress = np.asarray(self.strain, float), np.asarray(self.stress, float)
        toughness(strain, stress)
        object.__setattr__(self, "strain", strain.copy())
        object.__setattr__(self, "stress", stress.copy())

    def young_modulus(self, n=5):
        """Slope of a least-squares line through the first ``n`` samples.

        Parameters
        ----------
        n : int, optional
            Number of leading samples to fit, between 2 and the number of
            samples. They should all be elastic.

        Returns
        -------
        float
            Young's modulus [Pa].
        """
        _integer("n", n, 2)
        if n > len(self.strain):
            raise ValueError("n cannot exceed the number of samples")
        return float(np.polyfit(self.strain[:n], self.stress[:n], 1)[0])

    def toughness(self):
        """Area under the curve [J/m^3]."""
        from engineeringkit.materials.systems.elasticity import toughness

        return toughness(self.strain, self.stress)

    def offset_yield_strength(self, offset=0.002, n=5):
        """Stress where the curve crosses the elastic line shifted by ``offset`` strain.

        Parameters
        ----------
        offset : float, optional
            Plastic strain offset [-], positive (0.2% by default).
        n : int, optional
            Samples used for the elastic slope, as in :meth:`young_modulus`.

        Returns
        -------
        float
            Offset yield strength [Pa].

        Raises
        ------
        ValueError
            If the curve never crosses the offset line.
        """
        _positive(offset=offset)
        E = self.young_modulus(n)
        gap = self.stress - E * (self.strain - offset)
        crossing = np.flatnonzero((gap[:-1] > 0) & (gap[1:] <= 0))
        if not crossing.size:
            raise ValueError("the curve does not reach the offset line")
        i = crossing[0]
        t = gap[i] / (gap[i] - gap[i + 1])
        return float(self.stress[i] + t * (self.stress[i + 1] - self.stress[i]))

    def summary(self) -> dict:
        """Young's modulus [Pa], ultimate strength [Pa], strain at the end [-] and toughness [J/m^3]."""
        return {
            "young_modulus_Pa": self.young_modulus(min(5, len(self.strain))),
            "ultimate_strength_Pa": float(np.max(self.stress)),
            "final_strain": float(self.strain[-1]),
            "toughness_J_m3": self.toughness(),
        }

    def plot(self, ax=None):
        """Plot stress [Pa] against strain; return ``(fig, ax)``."""
        from engineeringkit.materials.visualizers.plots import plot_stress_strain

        ax = plot_stress_strain(self.strain, self.stress, ax=ax)
        return ax.figure, ax


@dataclass(frozen=True)
class CrackGrowthResult:
    """Fatigue crack length against cycles, from :func:`~engineeringkit.materials.crack_growth_life`."""

    crack_length: np.ndarray
    """ndarray: Crack lengths, increasing [m]."""

    cycles: np.ndarray
    """ndarray: Cycles needed to grow from the first crack length to each one."""

    @property
    def life(self) -> float:
        """float: Cycles to grow from the initial to the final crack length."""
        return float(self.cycles[-1])

    def cycles_to(self, a):
        """Cycles to reach a crack length, interpolated along the growth history.

        Parameters
        ----------
        a : float or array_like
            Crack length [m], between the initial and final lengths.

        Returns
        -------
        float or ndarray
            Load cycles [-].
        """
        a = np.asarray(a, float)
        if np.any(a < self.crack_length[0]) or np.any(a > self.crack_length[-1]):
            raise ValueError("a must lie between the initial and final crack lengths")
        return np.interp(a, self.crack_length, self.cycles)

    def summary(self) -> dict:
        """Initial and final crack lengths [m], life [cycles] and the cycles spent doubling the initial crack."""
        a0, af = float(self.crack_length[0]), float(self.crack_length[-1])
        return {
            "initial_crack_m": a0,
            "final_crack_m": af,
            "life_cycles": self.life,
            "cycles_to_double": float(self.cycles_to(2 * a0)) if 2 * a0 <= af else self.life,
        }

    def plot(self, ax=None):
        """Plot crack length [mm] against cycles; return ``(fig, ax)``."""
        from engineeringkit.materials.visualizers.plots import plot_crack_growth

        ax = plot_crack_growth(self.cycles, self.crack_length, ax=ax)
        return ax.figure, ax

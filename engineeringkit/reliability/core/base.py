"""Result containers for engineeringkit.reliability.

Most reliability models are closed-form probabilities that return plain
floats or arrays. The three that produce structured output, the
Kaplan-Meier estimator, the Hasofer-Lind first-order reliability method and
event-tree quantification, return the frozen dataclasses below with the
usual ``summary()`` and, where a picture helps, ``plot()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["KaplanMeierResult", "HasoferLindResult", "EventTreeResult"]


@dataclass(frozen=True)
class KaplanMeierResult:
    """Output of :func:`engineeringkit.reliability.kaplan_meier`."""

    times: np.ndarray
    """ndarray: Distinct failure times [s], increasing."""

    survival: np.ndarray
    """ndarray: Product-limit survival estimate just after each failure time [-]."""

    at_risk: np.ndarray
    """ndarray of int: Units still on test just before each failure time."""

    failures: np.ndarray
    """ndarray of int: Failures observed at each failure time."""

    std_error: np.ndarray
    """ndarray: Greenwood standard error of the survival estimate [-]."""

    n: int
    """int: Number of units on test."""

    censored: int
    """int: Number of units removed from test without failing."""

    def survival_at(self, t):
        """Step-function survival estimate, one before the first failure.

        Parameters
        ----------
        t : float or array_like
            Time(s) [s].

        Returns
        -------
        float or ndarray
            Survival probability [-].
        """
        index = np.searchsorted(self.times, np.asarray(t, dtype=float), side="right")
        return np.concatenate(([1.0], self.survival))[index]

    def summary(self) -> dict:
        """Units on test, failures, censored units and the last survival estimate."""
        return {
            "n": self.n,
            "failures": int(self.failures.sum()),
            "censored": self.censored,
            "final_survival": float(self.survival[-1]) if self.survival.size else 1.0,
        }

    def plot(self, ax=None, label="Kaplan-Meier"):
        """Draw the survival step function; return ``(fig, ax)``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created if omitted.
        label : str, optional
            Legend label of the curve.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """
        from engineeringkit.reliability.visualizers.plots import plot_survival

        ax = plot_survival(np.concatenate(([0.0], self.times)), np.concatenate(([1.0], self.survival)), ax=ax, label=label)
        return ax.figure, ax


@dataclass(frozen=True)
class HasoferLindResult:
    """Output of :func:`engineeringkit.reliability.hasofer_lind_index`."""

    beta: float
    """float: Hasofer-Lind reliability index, the distance to the limit state in standard normal space [-]."""

    design_point: np.ndarray
    """ndarray: Most probable failure point in the original variables (same units as ``mean``)."""

    alpha: np.ndarray
    """ndarray: Unit sensitivity vector; ``alpha[i]**2`` is variable ``i``'s share of the uncertainty [-]."""

    iterations: int
    """int: HL-RF iterations used."""

    @property
    def failure_probability(self) -> float:
        r"""First-order failure probability :math:`\Phi(-\beta)` [-]."""
        from engineeringkit.reliability.utils.normal import failure_probability_from_index

        return float(failure_probability_from_index(self.beta))

    def summary(self) -> dict:
        """Reliability index, first-order failure probability and iteration count."""
        return {"beta": self.beta, "failure_probability": self.failure_probability, "iterations": self.iterations}


@dataclass(frozen=True)
class EventTreeResult:
    """Output of :func:`engineeringkit.reliability.event_tree_probability`."""

    outcomes: np.ndarray
    """ndarray of bool, shape ``(2**m, m)``: ``True`` where a safety function fails in a sequence."""

    frequencies: np.ndarray
    """ndarray: Frequency of each sequence [1/s] (or probability if the initiator was a probability)."""

    def frequency(self, failed):
        """Total frequency of the sequences in which the functions listed in ``failed`` all fail.

        Parameters
        ----------
        failed : sequence of int
            Indices of the safety functions that must fail.

        Returns
        -------
        float
            Summed sequence frequency [1/s].
        """
        failed = list(failed)
        for i in failed:
            if not isinstance(i, (int, np.integer)) or not 0 <= i < self.outcomes.shape[1]:
                raise ValueError("failed must list valid safety-function indices")
        return float(self.frequencies[np.all(self.outcomes[:, failed], axis=1)].sum())

    def summary(self) -> dict:
        """Number of sequences, total frequency and the frequency with every function successful."""
        return {
            "sequences": int(self.frequencies.size),
            "total_frequency": float(self.frequencies.sum()),
            "all_success_frequency": float(self.frequencies[0]),
        }

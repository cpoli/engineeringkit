"""Result containers for engineeringkit.reliability.

Most reliability models are closed-form probabilities that return plain
floats or arrays. Those that produce structured output, the Kaplan-Meier
estimator, Weibull fits to life data, the Hasofer-Lind first-order
reliability method, Monte Carlo simulation, event-tree quantification, FMEA
ranking and Bayesian updating of a failure rate, return the frozen
dataclasses below with the usual ``summary()`` and, where a picture helps,
``plot()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "KaplanMeierResult",
    "WeibullFitResult",
    "HasoferLindResult",
    "MonteCarloResult",
    "EventTreeResult",
    "FMEAResult",
    "FailureRatePosterior",
]


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
class WeibullFitResult:
    """Output of :func:`~engineeringkit.reliability.weibull_rank_regression` and :func:`~engineeringkit.reliability.weibull_mle`."""

    eta: float
    """float: Fitted scale, the characteristic life (63.2% failed) [s]."""

    beta: float
    """float: Fitted shape [-]."""

    method: str
    """str: ``"rank regression on y"``, ``"rank regression on x"`` or ``"maximum likelihood"``."""

    failure_times: np.ndarray
    """ndarray: Failure times in increasing order [s]; censored units are not listed."""

    median_ranks: np.ndarray
    """ndarray: Median-rank estimate of the failure probability at each failure time [-]."""

    units: int
    """int: Number of units on test, censored ones included."""

    def b_life(self, p):
        r"""Time by which a fraction ``p`` of the population has failed, :math:`\eta(-\ln(1-p))^{1/\beta}`.

        Parameters
        ----------
        p : float or array_like
            Failure fraction [-], in (0, 1); ``0.1`` gives the B10 life.

        Returns
        -------
        float or ndarray
            Life [s].
        """
        from engineeringkit._validation import interval

        interval("p", p, 0, 1, closed=False)
        return self.eta * (-np.log1p(-np.asarray(p, dtype=float))) ** (1 / self.beta)

    def summary(self) -> dict:
        """Scale, shape, B10 life, units on test and failures."""
        return {
            "eta": self.eta,
            "beta": self.beta,
            "b10_life": float(self.b_life(0.1)),
            "n": self.units,
            "failures": int(self.failure_times.size),
        }

    def plot(self, ax=None):
        r"""Draw the data and the fitted line on Weibull probability axes; return ``(fig, ax)``.

        The axes are :math:`\ln t` against :math:`\ln(-\ln(1 - F))`, on which
        a Weibull distribution is a straight line of slope :math:`\beta`.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created if omitted.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """
        from engineeringkit.reliability.visualizers.plots import plot_weibull_probability

        ax = plot_weibull_probability(self.failure_times, self.median_ranks, self.eta, self.beta, ax=ax, label=self.method)
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
class MonteCarloResult:
    """Output of :func:`engineeringkit.reliability.monte_carlo_failure_probability`."""

    failure_probability: float
    """float: Estimated failure probability [-]."""

    std_error: float
    """float: Standard error of the estimate [-]."""

    samples: int
    """int: Number of samples drawn."""

    failures: int
    """int: Number of samples that fell in the failure domain."""

    importance_sampling: bool
    """bool: ``True`` if samples were drawn around a shifted centre and weighted."""

    @property
    def cov(self) -> float:
        """Coefficient of variation of the estimate, standard error over estimate [-] (infinite with no failures)."""
        return self.std_error / self.failure_probability if self.failure_probability > 0 else float("inf")

    @property
    def beta(self) -> float:
        r"""Reliability index :math:`-\Phi^{-1}(p_f)` equivalent to the estimate [-]."""
        from scipy import special

        return float(-special.ndtri(self.failure_probability))

    def confidence_interval(self, confidence=0.95):
        """Normal-approximation confidence interval on the failure probability, clipped at zero.

        Parameters
        ----------
        confidence : float, optional
            Two-sided confidence level [-], in (0, 1).

        Returns
        -------
        tuple of float
            ``(lower, upper)`` [-].
        """
        from scipy import special

        from engineeringkit._validation import interval

        interval("confidence", confidence, 0, 1, closed=False)
        half = float(special.ndtri(0.5 + confidence / 2)) * self.std_error
        return max(self.failure_probability - half, 0.0), self.failure_probability + half

    def summary(self) -> dict:
        """Estimate, standard error, coefficient of variation, sample and failure counts."""
        return {
            "failure_probability": self.failure_probability,
            "std_error": self.std_error,
            "cov": self.cov,
            "samples": self.samples,
            "failures": self.failures,
        }


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


@dataclass(frozen=True)
class FMEAResult:
    """Output of :func:`engineeringkit.reliability.fmea_ranking`, failure modes in decreasing risk priority."""

    modes: tuple
    """tuple of str: Failure-mode names, highest risk priority number first."""

    severity: np.ndarray
    """ndarray of int: Severity ranking of each mode, 1 to 10."""

    occurrence: np.ndarray
    """ndarray of int: Occurrence ranking of each mode, 1 to 10."""

    detection: np.ndarray
    """ndarray of int: Detection ranking of each mode, 1 (certain detection) to 10."""

    rpn: np.ndarray
    """ndarray of int: Risk priority number, severity x occurrence x detection."""

    def summary(self) -> dict:
        """Number of modes, the top mode, its RPN and the total RPN."""
        return {
            "modes": len(self.modes),
            "top_mode": self.modes[0],
            "top_rpn": int(self.rpn[0]),
            "total_rpn": int(self.rpn.sum()),
        }

    def plot(self, ax=None):
        """Draw a Pareto chart of the risk priority numbers; return ``(fig, ax)``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created if omitted.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        positions = np.arange(len(self.modes))
        ax.bar(positions, self.rpn, color="C0")
        ax.set_xticks(positions, self.modes, rotation=30, ha="right")
        share = ax.twinx()
        share.plot(positions, 100 * np.cumsum(self.rpn) / self.rpn.sum(), "C3-o")
        share.set(ylabel="cumulative share of RPN [%]", ylim=(0, 105))
        ax.set(ylabel="risk priority number", title="FMEA Pareto chart")
        return ax.figure, ax


@dataclass(frozen=True)
class FailureRatePosterior:
    """Output of :func:`engineeringkit.reliability.bayesian_failure_rate`: a gamma distribution on the failure rate."""

    shape: float
    """float: Posterior gamma shape, prior shape plus the failures observed [-]."""

    rate: float
    """float: Posterior gamma rate, prior rate plus the time on test [s]."""

    prior_shape: float
    """float: Prior gamma shape [-]."""

    prior_rate: float
    """float: Prior gamma rate [s]."""

    @property
    def mean(self) -> float:
        """Posterior mean failure rate, shape over rate [1/s]."""
        return self.shape / self.rate

    @property
    def prior_mean(self) -> float:
        """Prior mean failure rate [1/s]."""
        return self.prior_shape / self.prior_rate

    def upper_bound(self, confidence=0.9):
        """One-sided upper credible bound on the failure rate.

        Parameters
        ----------
        confidence : float, optional
            Posterior probability that the rate lies below the bound [-], in (0, 1).

        Returns
        -------
        float
            Failure rate [1/s].
        """
        from scipy import special

        from engineeringkit._validation import interval

        interval("confidence", confidence, 0, 1, closed=False)
        return float(special.gammaincinv(self.shape, confidence)) / self.rate

    def credible_interval(self, confidence=0.9):
        """Equal-tailed credible interval on the failure rate.

        Parameters
        ----------
        confidence : float, optional
            Posterior probability inside the interval [-], in (0, 1).

        Returns
        -------
        tuple of float
            ``(lower, upper)`` [1/s].
        """
        from scipy import special

        from engineeringkit._validation import interval

        interval("confidence", confidence, 0, 1, closed=False)
        tail = (1 - confidence) / 2
        return float(special.gammaincinv(self.shape, tail)) / self.rate, float(special.gammaincinv(self.shape, 1 - tail)) / self.rate

    def summary(self) -> dict:
        """Prior and posterior mean rates, posterior shape and rate, and the 90% upper bound."""
        return {
            "prior_mean": self.prior_mean,
            "mean": self.mean,
            "shape": self.shape,
            "rate": self.rate,
            "upper_90": self.upper_bound(0.9),
        }

    def plot(self, ax=None):
        """Draw the prior and posterior densities of the failure rate; return ``(fig, ax)``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created if omitted.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """
        import matplotlib.pyplot as plt
        from scipy import stats

        if ax is None:
            _, ax = plt.subplots()
        top = max(
            stats.gamma.ppf(0.999, self.shape, scale=1 / self.rate), stats.gamma.ppf(0.99, self.prior_shape, scale=1 / self.prior_rate)
        )
        lam = np.linspace(0, top, 400)
        ax.plot(lam, stats.gamma.pdf(lam, self.prior_shape, scale=1 / self.prior_rate), "--", label="prior")
        ax.plot(lam, stats.gamma.pdf(lam, self.shape, scale=1 / self.rate), label="posterior")
        ax.set(xlabel="failure rate [1/s]", ylabel="probability density")
        ax.legend()
        return ax.figure, ax

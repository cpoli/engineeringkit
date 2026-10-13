"""Inspectable Euler, RK4 and adaptive integration of SI state models."""

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class IntegrationResult:
    """Time [s] and state history with shape (samples, states)."""

    t: np.ndarray
    history: np.ndarray
    method: str
    evaluations: int

    def summary(self) -> dict:
        return {"method": self.method, "samples": len(self.t), "evaluations": self.evaluations, "final_state": self.history[-1].tolist()}

    def plot(self, labels=None, ax=None):
        """Plot state histories; labels should include their physical units."""
        from engineeringkit.visualizers import plot_history

        return plot_history(self.t, self.history, labels=labels, ax=ax)


def integrate_ode(rhs: Callable, y0, t, method: str = "RK45", *, rtol: float = 1e-8, atol: float = 1e-10) -> IntegrationResult:
    """Solve dy/dt = rhs(t, y) on a strictly increasing finite time grid.

    ``euler`` and ``rk4`` expose explicit fixed-step methods using the supplied
    grid. Other names select SciPy adaptive methods. Explicit methods are for
    non-stiff teaching problems; step size controls stability and accuracy.
    State units follow the caller's model; time is in seconds.
    """
    t = np.asarray(t, dtype=float)
    y0 = np.atleast_1d(np.asarray(y0, dtype=float))
    if t.ndim != 1 or len(t) < 2 or not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError("t must contain at least two finite, strictly increasing times")
    if y0.ndim != 1 or not y0.size or not np.all(np.isfinite(y0)):
        raise ValueError("y0 must be a nonempty finite state vector")
    if not np.isfinite(rtol) or not np.isfinite(atol) or rtol <= 0 or atol <= 0:
        raise ValueError("tolerances must be finite and positive")
    evaluations = 0

    def checked_rhs(time, state):
        nonlocal evaluations
        evaluations += 1
        derivative = np.asarray(rhs(time, state), dtype=float)
        if derivative.shape != y0.shape or not np.all(np.isfinite(derivative)):
            raise ValueError("rhs must return a finite derivative matching y0")
        return derivative

    if method.lower() in {"euler", "rk4"}:
        history = np.empty((len(t), len(y0)))
        history[0] = y0
        for i, dt in enumerate(np.diff(t)):
            ti, yi = t[i], history[i]
            k1 = checked_rhs(ti, yi)
            if method.lower() == "euler":
                history[i + 1] = yi + dt * k1
            else:
                k2 = checked_rhs(ti + dt / 2, yi + dt * k1 / 2)
                k3 = checked_rhs(ti + dt / 2, yi + dt * k2 / 2)
                k4 = checked_rhs(ti + dt, yi + dt * k3)
                history[i + 1] = yi + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            if not np.all(np.isfinite(history[i + 1])):
                raise RuntimeError("non-finite state; reduce step size")
    else:
        solution = solve_ivp(checked_rhs, (t[0], t[-1]), y0, t_eval=t, method=method, rtol=rtol, atol=atol)
        if not solution.success:
            raise RuntimeError(solution.message)
        history = solution.y.T
    return IntegrationResult(t.copy(), history, method, evaluations)


__all__ = ["IntegrationResult", "integrate_ode"]

"""Concrete control models: transfer functions, frequency response, stability, state space, optimal control and PID."""

from engineeringkit.control_systems.systems.frequency import bode, nyquist, nyquist_encirclements, stability_margins
from engineeringkit.control_systems.systems.optimal import dlqr, kalman_gain, lqr
from engineeringkit.control_systems.systems.pid import PID, ziegler_nichols
from engineeringkit.control_systems.systems.root_locus import root_locus
from engineeringkit.control_systems.systems.sampled import discretize_zoh
from engineeringkit.control_systems.systems.second_order import second_order_metrics
from engineeringkit.control_systems.systems.stability import (
    hurwitz_determinants,
    lyapunov_equation,
    lyapunov_stable,
    routh_first_column,
    routh_stable,
)
from engineeringkit.control_systems.systems.state_space import (
    controllability_matrix,
    feedforward_gain,
    is_controllable,
    is_observable,
    observability_matrix,
    observer_gain,
    state_feedback_poles,
    state_space,
)
from engineeringkit.control_systems.systems.transfer import (
    dc_gain,
    feedback,
    impulse_response,
    poles,
    step_response,
    transfer_function,
    zeros,
)

__all__ = [
    "transfer_function",
    "feedback",
    "poles",
    "zeros",
    "dc_gain",
    "step_response",
    "impulse_response",
    "bode",
    "nyquist",
    "nyquist_encirclements",
    "stability_margins",
    "routh_first_column",
    "routh_stable",
    "hurwitz_determinants",
    "lyapunov_equation",
    "lyapunov_stable",
    "root_locus",
    "second_order_metrics",
    "PID",
    "ziegler_nichols",
    "state_space",
    "controllability_matrix",
    "observability_matrix",
    "is_controllable",
    "is_observable",
    "state_feedback_poles",
    "observer_gain",
    "feedforward_gain",
    "discretize_zoh",
    "lqr",
    "dlqr",
    "kalman_gain",
]

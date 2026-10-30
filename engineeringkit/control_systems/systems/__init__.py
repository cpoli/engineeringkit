"""Concrete control models: transfer functions, frequency response, stability, state space, optimal control and PID."""

from engineeringkit.control_systems.systems.compensators import (
    design_lag,
    design_lead,
    design_lead_lag,
    lag_compensator,
    lead_compensator,
    lead_lag_compensator,
)
from engineeringkit.control_systems.systems.delay import pade_delay, smith_predictor
from engineeringkit.control_systems.systems.frequency import (
    bode,
    nichols,
    nichols_closed_loop,
    nyquist,
    nyquist_encirclements,
    stability_margins,
)
from engineeringkit.control_systems.systems.nonlinear import (
    dead_zone_describing_function,
    describing_function,
    limit_cycle,
    phase_portrait,
    relay_describing_function,
    saturation_describing_function,
)
from engineeringkit.control_systems.systems.optimal import dlqr, kalman_gain, lqr
from engineeringkit.control_systems.systems.pid import PID, ziegler_nichols
from engineeringkit.control_systems.systems.root_locus import root_locus
from engineeringkit.control_systems.systems.sampled import discretize_tustin, discretize_zoh, jury_stable
from engineeringkit.control_systems.systems.second_order import second_order_metrics
from engineeringkit.control_systems.systems.sensitivity import (
    complementary_sensitivity,
    sensitivity,
    sensitivity_peak,
)
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
from engineeringkit.control_systems.systems.steady_state import (
    error_constants,
    system_type,
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
    "nichols",
    "nichols_closed_loop",
    "routh_first_column",
    "routh_stable",
    "hurwitz_determinants",
    "lyapunov_equation",
    "lyapunov_stable",
    "root_locus",
    "second_order_metrics",
    "pade_delay",
    "smith_predictor",
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
    "discretize_tustin",
    "jury_stable",
    "lqr",
    "dlqr",
    "kalman_gain",
    "lead_compensator",
    "lag_compensator",
    "lead_lag_compensator",
    "design_lead",
    "design_lag",
    "design_lead_lag",
    "system_type",
    "error_constants",
    "sensitivity",
    "complementary_sensitivity",
    "sensitivity_peak",
    "describing_function",
    "relay_describing_function",
    "saturation_describing_function",
    "dead_zone_describing_function",
    "limit_cycle",
    "phase_portrait",
]

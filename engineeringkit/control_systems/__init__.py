"""engineeringkit.control_systems: feedback, stability, frequency response, state space and optimal control.

Black's feedback amplifier and Minorsky's PID controller with
Ziegler-Nichols tuning; Maxwell's characteristic roots, the Routh and
Hurwitz tests and Lyapunov's direct method; Nyquist's encirclement
criterion, Bode plots, stability margins, the Nichols chart and Bode's
sensitivity functions; lead, lag and lead-lag compensators, system type and
error constants; Evans's root locus; Hazen's second-order servo
specifications; Padé dead-time approximation and Smith's predictor;
zero-order-hold and Tustin discretisation and Jury's test; PID anti-windup;
describing functions, limit cycles and the phase plane; Kalman's
controllability and observability, pole placement and Luenberger observers;
and the LQR, discrete LQR and Kalman-Bucy Riccati designs. Models are
continuous-time and SISO unless stated; SciPy supplies the LTI objects,
Riccati/Lyapunov solvers, pole placement and step responses.
"""

from engineeringkit import __version__
from engineeringkit.control_systems.core.base import CompensatorResult, PhasePortraitResult, ResponseResult
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
from engineeringkit.control_systems.utils.matrices import state_matrices, weight_matrix
from engineeringkit.control_systems.utils.polynomials import as_polynomial, num_den

__all__ = [
    "__version__",
    "ResponseResult",
    "CompensatorResult",
    "PhasePortraitResult",
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
    "state_matrices",
    "weight_matrix",
    "as_polynomial",
    "num_den",
]

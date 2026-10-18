"""engineeringkit.reliability: lifetimes, redundancy, structural reliability, risk and maintenance.

Halley's life table and the bathtub hazard; Poisson failure counts;
exponential and Weibull lifetimes with Epstein-Sobel confidence bounds and
the Kaplan-Meier estimator for censored data; series, parallel and
k-out-of-n systems with Birnbaum importance; Freudenthal's stress-strength
interference and the Hasofer-Lind index; Watson's fault trees, WASH-1400
event trees and IEC 61508 safety integrity levels; Markov availability,
Barlow-Proschan age replacement and Duane reliability growth. Models are
hand-written closed forms; SciPy supplies special functions, chi-square
quantiles and the matrix exponential.
"""

from engineeringkit import __version__
from engineeringkit.reliability.core.base import EventTreeResult, HasoferLindResult, KaplanMeierResult
from engineeringkit.reliability.systems.architecture import (
    birnbaum_importance_series,
    k_of_n_reliability,
    parallel_reliability,
    series_reliability,
)
from engineeringkit.reliability.systems.growth import (
    reliability_growth_duane,
)
from engineeringkit.reliability.systems.lifetimes import (
    bathtub_hazard,
    expected_failures,
    exponential_failure_probability,
    exponential_hazard,
    exponential_mttf,
    exponential_reliability,
    failure_rate_from_mtbf,
    kaplan_meier,
    life_table_hazard,
    mtbf_confidence_bounds,
    poisson_failures_probability,
    weibull_cdf,
    weibull_hazard,
    weibull_mttf,
    weibull_reliability,
)
from engineeringkit.reliability.systems.maintenance import (
    age_replacement_cost_rate,
    availability,
    maintainability,
    markov_availability,
    markov_state_probabilities,
    renewal_interval_cost_rate,
    repair_rate,
    unavailability,
)
from engineeringkit.reliability.systems.risk import (
    event_tree_probability,
    fault_tree_probability,
    safety_integrity_level,
    safety_integrity_pfd,
)
from engineeringkit.reliability.systems.structural import (
    hasofer_lind_index,
    stress_strength_reliability,
)
from engineeringkit.reliability.utils.normal import failure_probability_from_index, reliability_index_from_probability

__all__ = [
    "__version__",
    "KaplanMeierResult",
    "HasoferLindResult",
    "EventTreeResult",
    "exponential_reliability",
    "exponential_failure_probability",
    "exponential_mttf",
    "exponential_hazard",
    "failure_rate_from_mtbf",
    "weibull_reliability",
    "weibull_cdf",
    "weibull_hazard",
    "weibull_mttf",
    "bathtub_hazard",
    "life_table_hazard",
    "expected_failures",
    "poisson_failures_probability",
    "mtbf_confidence_bounds",
    "kaplan_meier",
    "series_reliability",
    "parallel_reliability",
    "k_of_n_reliability",
    "birnbaum_importance_series",
    "stress_strength_reliability",
    "hasofer_lind_index",
    "fault_tree_probability",
    "event_tree_probability",
    "safety_integrity_pfd",
    "safety_integrity_level",
    "availability",
    "unavailability",
    "maintainability",
    "repair_rate",
    "markov_availability",
    "markov_state_probabilities",
    "age_replacement_cost_rate",
    "renewal_interval_cost_rate",
    "reliability_growth_duane",
    "failure_probability_from_index",
    "reliability_index_from_probability",
]

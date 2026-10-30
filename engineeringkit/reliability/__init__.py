"""engineeringkit.reliability: lifetimes, life data, redundancy, structural reliability, risk and maintenance.

Halley's life table and the bathtub hazard; Poisson failure counts;
exponential, Weibull, lognormal and gamma lifetimes with Epstein-Sobel
confidence bounds and the Kaplan-Meier estimator for censored data; Weibull
parameters by median-rank regression and maximum likelihood; Arrhenius,
Eyring and inverse-power-law accelerated life testing; Bayesian updating
of a failure rate; series, parallel, k-out-of-n, cold-standby and general
block-diagram systems with minimal cut sets and Birnbaum importance;
Freudenthal's stress-strength interference, the Hasofer-Lind index and
Monte Carlo; FMEA, Watson's fault trees, WASH-1400 event trees and IEC
61508 safety integrity levels; Markov availability, Barlow-Proschan age
replacement and Duane reliability growth. Models are
hand-written closed forms; SciPy supplies special functions, chi-square
quantiles and the matrix exponential.
"""

from engineeringkit import __version__
from engineeringkit.reliability.core.base import (
    EventTreeResult,
    FailureRatePosterior,
    FMEAResult,
    HasoferLindResult,
    KaplanMeierResult,
    MonteCarloResult,
    WeibullFitResult,
)
from engineeringkit.reliability.systems.acceleration import (
    arrhenius_acceleration_factor,
    arrhenius_fit,
    eyring_acceleration_factor,
    inverse_power_law_acceleration_factor,
)
from engineeringkit.reliability.systems.architecture import (
    birnbaum_importance_series,
    esary_proschan_bounds,
    k_of_n_reliability,
    minimal_cut_sets,
    minimal_path_sets,
    network_reliability,
    parallel_reliability,
    series_reliability,
    standby_mttf,
    standby_reliability,
)
from engineeringkit.reliability.systems.bayesian import (
    bayesian_failure_rate,
    gamma_prior_from_bounds,
)
from engineeringkit.reliability.systems.growth import (
    reliability_growth_duane,
)
from engineeringkit.reliability.systems.life_data import (
    adjusted_ranks,
    median_ranks,
    weibull_mle,
    weibull_rank_regression,
)
from engineeringkit.reliability.systems.lifetimes import (
    bathtub_hazard,
    expected_failures,
    exponential_failure_probability,
    exponential_hazard,
    exponential_mttf,
    exponential_reliability,
    failure_rate_from_mtbf,
    gamma_cdf,
    gamma_hazard,
    gamma_mttf,
    gamma_reliability,
    kaplan_meier,
    life_table_hazard,
    lognormal_cdf,
    lognormal_hazard,
    lognormal_mttf,
    lognormal_reliability,
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
    failure_mode_criticality,
    fault_tree_probability,
    fmea_ranking,
    risk_priority_number,
    safety_integrity_level,
    safety_integrity_pfd,
)
from engineeringkit.reliability.systems.structural import (
    hasofer_lind_index,
    monte_carlo_failure_probability,
    stress_strength_reliability,
)
from engineeringkit.reliability.utils.normal import failure_probability_from_index, reliability_index_from_probability

__all__ = [
    "__version__",
    "KaplanMeierResult",
    "WeibullFitResult",
    "HasoferLindResult",
    "MonteCarloResult",
    "EventTreeResult",
    "FMEAResult",
    "FailureRatePosterior",
    "exponential_reliability",
    "exponential_failure_probability",
    "exponential_mttf",
    "exponential_hazard",
    "failure_rate_from_mtbf",
    "weibull_reliability",
    "weibull_cdf",
    "weibull_hazard",
    "weibull_mttf",
    "lognormal_reliability",
    "lognormal_cdf",
    "lognormal_hazard",
    "lognormal_mttf",
    "gamma_reliability",
    "gamma_cdf",
    "gamma_hazard",
    "gamma_mttf",
    "median_ranks",
    "adjusted_ranks",
    "weibull_rank_regression",
    "weibull_mle",
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
    "minimal_path_sets",
    "minimal_cut_sets",
    "network_reliability",
    "esary_proschan_bounds",
    "standby_reliability",
    "standby_mttf",
    "stress_strength_reliability",
    "hasofer_lind_index",
    "monte_carlo_failure_probability",
    "risk_priority_number",
    "fmea_ranking",
    "failure_mode_criticality",
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
    "arrhenius_acceleration_factor",
    "arrhenius_fit",
    "eyring_acceleration_factor",
    "inverse_power_law_acceleration_factor",
    "bayesian_failure_rate",
    "gamma_prior_from_bounds",
    "reliability_growth_duane",
    "failure_probability_from_index",
    "reliability_index_from_probability",
]

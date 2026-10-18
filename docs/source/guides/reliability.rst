Reliability: distributions, independence and rate units
=======================================================

Times are seconds; failure and repair rates are 1/s. Use consistent conversion
before mixing hour-based field data with SI parameters. Probabilities must lie
in [0,1]. Exponential reliability accepts zero rate (perfect survival), while
MTTF requires positive rate because the zero-rate mean is infinite.

Weibull and numerical endpoints
-------------------------------

R(t)=exp(-(t/eta)^beta), eta>0, beta>0. The hazard at t=0 is infinity for beta<1,
1/eta for beta=1, and zero for beta>1. Infinity here is a genuine distribution
limit, not invalid input. CDF helpers use expm1 to preserve very small failure
probabilities. Poisson probabilities use log-gamma to avoid overflow in large
factorials; k=0 with zero expected failures has probability one.

.. code-block:: python

    from engineeringkit.reliability import weibull_mttf, weibull_reliability
    survival = weibull_reliability([0, 3, 6], eta=3, beta=2)
    mean = weibull_mttf(3, 2)  # approximately 2.65868 s

Tests compare Weibull survival with SciPy's distribution, integrate survival
to recover MTTF, enumerate independent component states for system reliability,
and check k-of-n against a binomial tail. Birnbaum sensitivity is checked with
a finite difference rather than only its product formula.

System scope
------------

Series/parallel formulas assume independent component failures and perfect
switching for parallel operation. k-of-n additionally assumes identical component
reliability. Common-cause failures, standby dynamics and repair during a mission
require different models. Availability=MTBF/(MTBF+MTTR) is steady-state average
availability for alternating failure/repair cycles, not mission survival.

The PFDavg helper lambda_du*T/2 assumes rare dangerous undetected failures and
perfect proof tests; this implementation requires lambda_du*T<=0.1. It does not
certify a safety integrity level, which also needs architectural constraints
and systematic capability. renewal_interval_cost_rate is the constant-hazard
baseline preventive_cost/T+failure_cost*lambda; age_replacement_cost_rate is
Barlow and Proschan's renewal-reward ratio for a Weibull part. The Duane
helper returns cumulative MTBF K*t^alpha; K has units s^(1-alpha), not a
dimensionless gain, and the instantaneous MTBF is K*t^alpha/(1-alpha).

Life data
---------

life_table_hazard returns the average hazard -ln(S(t+dt)/S(t))/dt over each
interval of a survivor table, so a constant hazard is recovered exactly at any
interval width. kaplan_meier treats a unit censored at the same time as a
failure as still at risk at that failure (the usual convention) and gives
Greenwood's standard error. mtbf_confidence_bounds uses chi-square quantiles
with 2r (failure-terminated) or 2r+2 (time-terminated) degrees of freedom for
the lower bound and 2r for the upper, which is infinite with no failures; it
assumes exponential lives.

Structural reliability and risk
-------------------------------

stress_strength_reliability assumes independent normal strength and stress.
hasofer_lind_index assumes independent normal variables, runs the HL-RF
iteration with a central-difference gradient, and returns a negative index
when the mean point already lies in the failure region; Phi(-beta) is a
first-order approximation that is exact only for a linear limit state.
fault_tree_probability combines gates directly when no basic event repeats and
otherwise enumerates every state of the basic events (at most 20), so it is
exact for independent basic events. event_tree_probability keeps every branch
and assumes independent branch points; real event trees prune branches and
condition later functions on earlier ones.

Repairable systems
------------------

markov_availability is the two-state model with constant failure and repair
rates, starting up. markov_state_probabilities accepts any generator matrix
(nonnegative off-diagonal rates, zero row sums) and uses SciPy's matrix
exponential; the up states must be summed by the caller.

Reference: `NIST Weibull reliability handbook <https://itl.nist.gov/div898/handbook/apr/section1/apr162.htm>`_.

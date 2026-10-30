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
reliability. Common-cause failures and repair during a mission require different
models.

A general block diagram is a mapping ``{name: (node, node)}``: each component is
an undirected link, and the system works if working links join ``source`` to
``sink``. minimal_path_sets lists the simple paths; minimal_cut_sets finds the
minimal sets meeting every path, by size (at most 20 components on paths);
network_reliability enumerates every component state and is exact for
independent components. esary_proschan_bounds returns the cut-set lower and
path-set upper bounds; the lower bound is exact for a series system and very
close for reliable components, the upper bound exact for a parallel one.

standby_reliability is cold standby of identical units: spares do not age while
waiting, and each changeover succeeds with probability p. With p = 1 the system
life is gamma (Erlang) with shape n, and standby_mttf is n/lambda. Warm standby
(spares that age) and unequal units are not modelled. Availability=MTBF/(MTBF+MTTR) is steady-state average
availability for alternating failure/repair cycles, not mission survival.

The PFDavg helper lambda_du*T/2 assumes rare dangerous undetected failures and
perfect proof tests; this implementation requires lambda_du*T<=0.1. It does not
certify a safety integrity level, which also needs architectural constraints
and systematic capability. renewal_interval_cost_rate is the constant-hazard
baseline preventive_cost/T+failure_cost*lambda; age_replacement_cost_rate is
Barlow and Proschan's renewal-reward ratio for a Weibull part. The Duane
helper returns cumulative MTBF K*t^alpha; K has units s^(1-alpha), not a
dimensionless gain, and the instantaneous MTBF is K*t^alpha/(1-alpha).

Lognormal and gamma lives
-------------------------

lognormal_* take mu and sigma of ln t (mu is the log of the median life in
seconds, so changing time units shifts mu, not sigma). The hazard is zero at
t = 0, peaks and then falls toward zero; it is formed in log space and stays
finite where the survival function underflows. gamma_* take a shape k and a
scale theta [s]; an integer k is k exponential stages in sequence. Their hazard
at t = 0 is infinite, 1/theta or zero for k below, at or above one, and tends
to 1/theta for every k.

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

median_ranks gives the exact beta-distribution median (Johnson) or Benard's
(i-0.3)/(n+0.4); they differ by less than 0.003 for n >= 5. adjusted_ranks
implements Johnson's reverse-rank formula O = (r*O_prev + n + 1)/(r + 1) and
ranks a failure before a censored unit at a tie. weibull_rank_regression fits
ln(-ln(1-F)) against ln t by least squares, regressing y on x by default or x
on y (``regress="x"``), and needs two failures at different times.
weibull_mle solves Cohen's profile equation for beta by bracketing and Brent's
method, after scaling times by the largest, and uses every censored time; it
needs a failure before the longest time on test. Both return a
WeibullFitResult with b_life(p) and a probability plot. Neither gives
confidence bounds on eta and beta, and both assume a two-parameter Weibull
(no location shift).

Accelerated testing and Bayesian rates
--------------------------------------

Activation energies are in joules: multiply electronvolts by
``scipy.constants.eV``. Temperatures are absolute (kelvin). The acceleration
factor is life in service over life on test, above one for a harsher test.
arrhenius_fit is a least-squares line of ln L against 1/T. Eyring's factor is
Arrhenius's times T_test/T_use; inverse_power_law_acceleration_factor is
(V_test/V_use)^n for any stress unit. Every model assumes the test does not
change the failure mechanism.

bayesian_failure_rate is the conjugate gamma-Poisson update for a constant
rate: shape a + r, rate b + T, where b has units of time. With a, b -> 0 its
upper bound equals the chi-square bound with 2r degrees of freedom, and with
a = 1, b -> 0 the time-terminated 2r + 2 bound. gamma_prior_from_bounds solves
for a by bisection on the ratio of the upper percentile to the median.

Structural reliability and risk
-------------------------------

stress_strength_reliability assumes independent normal strength and stress.
hasofer_lind_index assumes independent normal variables, runs the HL-RF
iteration with a central-difference gradient, and returns a negative index
when the mean point already lies in the failure region; Phi(-beta) is a
first-order approximation that is exact only for a linear limit state.
monte_carlo_failure_probability calls the limit state once on an array of
shape (m, samples), so the same ``x[0] - x[1]`` function serves both methods.
Crude sampling has standard error sqrt(p(1-p)/N); with ``center`` (the
Hasofer-Lind design point) it samples a unit normal shifted there and weights
by the density ratio. The confidence interval is the normal approximation and
is unreliable with only a handful of failures.

risk_priority_number checks that each ranking is an integer from 1 to 10;
the RPN is ordinal, so equal RPNs are not equal risks and high severities need
review on their own. fmea_ranking breaks ties by severity and then by name.
failure_mode_criticality is the MIL-STD-1629A product beta*alpha*lambda*t.
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

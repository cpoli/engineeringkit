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
certify a safety integrity level. The scheduled-maintenance cost helper is a
simple preventive_cost/T+failure_cost*lambda baseline, not an age-replacement
renewal model. The Duane helper returns cumulative MTBF K*t^alpha; K has units
s^(1-alpha), not a dimensionless gain.

Reference: `NIST Weibull reliability handbook <https://itl.nist.gov/div898/handbook/apr/section1/apr162.htm>`_.

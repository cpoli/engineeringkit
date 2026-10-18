Breakthroughs in Reliability
============================


.. include:: /_generated/nav/reliability.rst

.. epigraph::

   "Probabilistic logics and the synthesis of reliable organisms from
   unreliable components."
   -- John von Neumann, title of his 1956 paper

Reliability engineering asks how long things last and how to make systems
that last longer than their parts. Its mathematics began with life
insurance: the first life tables were estimates of a hazard rate. Poisson
supplied the law of rare events, and after the Second World War aircraft
structures, vacuum-tube electronics and missiles turned these ideas into a
discipline. Engineers learned to fit lifetime distributions to censored
test data, to combine components into series, parallel and voting systems,
to trace accidents through fault and event trees, and to schedule repair
and replacement. The entries trace the ideas behind
:mod:`engineeringkit.reliability`. Each links to the implementing function
and to a gallery example devoted to it.

.. contents:: Timeline
   :local:
   :depth: 1

1693 -- Halley's Life Table and the Hazard Rate
-----------------------------------------------

Edmond Halley obtained five years of births and burials from the city of
Breslau, a population with little migration, and tabulated how many of
1000 one-year-olds would still be alive at each later age. He used the
table to price life annuities. The fraction of those alive at one age who
die before the next is the probability of failing in that interval given
survival so far, which reliability engineers call the *hazard rate*.
Halley's hazard falls through childhood, is lowest around age ten and then
rises steadily: the bathtub shape that also describes many machines, with
early failures, a long period of random ones and wear-out.

*Implementation:* :func:`engineeringkit.reliability.life_table_hazard` and
:func:`engineeringkit.reliability.bathtub_hazard`. The example computes the
hazard from Halley's table and fits a bathtub curve through it.

*References:* E. Halley, "An estimate of the degrees of the mortality of
mankind, drawn from curious tables of the births and funerals at the city
of Breslaw," *Phil. Trans. R. Soc.* 17 (1693) 596-610.

.. minigallery:: ../../examples/reliability/lifetimes/plot_01_halley_life_table.py

1837 -- Poisson's Law of Rare Events
------------------------------------

Siméon Denis Poisson, studying the verdicts of juries, derived the limit of
the binomial distribution when there are many trials and each success is
rare: the number of events in a fixed exposure has probability
:math:`(\lambda T)^k e^{-\lambda T}/k!`. A fleet of units that fail
independently at a constant rate is exactly this situation, with
:math:`\lambda T` expected failures over a total operating time
:math:`T`. The distribution sizes spare-parts kits and underlies every
test of a constant failure rate.

*Implementation:* :func:`engineeringkit.reliability.poisson_failures_probability`
and :func:`engineeringkit.reliability.expected_failures`. The example
simulates failures with exponential gaps, checks the count distribution and
sizes a spares kit.

*References:* S. D. Poisson, *Recherches sur la probabilité des jugements
en matière criminelle et en matière civile* (Bachelier, Paris, 1837),
§81.

.. minigallery:: ../../examples/reliability/lifetimes/plot_02_poisson_failures.py

1947 -- Freudenthal's Stress-Strength Interference
--------------------------------------------------

Alfred Freudenthal argued that a structure's safety factor, mean strength
over mean load, cannot measure its safety: two designs with the same factor
fail at very different rates if their loads and strengths scatter
differently. Safety is the probability that strength exceeds load, which
depends on how far the two distributions overlap. For independent normal
variables it is :math:`\Phi(\beta)` with
:math:`\beta = (\mu_S - \mu_s)/\sqrt{\sigma_S^2 + \sigma_s^2}`. The paper
started probabilistic structural design.

*Implementation:* :func:`engineeringkit.reliability.stress_strength_reliability`.
The example checks the closed form by Monte Carlo and shows two designs
with the same safety factor and different reliabilities.

*References:* A. M. Freudenthal, "The safety of structures," *Trans. ASCE*
112 (1947) 125-159.

.. minigallery:: ../../examples/reliability/structural/plot_01_freudenthal_interference.py

1951 -- Weibull's Distribution
------------------------------

Waloddi Weibull proposed a distribution, :math:`R(t) = e^{-(t/\eta)^\beta}`,
for data as varied as the yield strength of steel, the size of fly ash and
the fatigue life of bearings. Its shape parameter sets the hazard: falling
for :math:`\beta < 1`, constant for :math:`\beta = 1` and rising above.
Because :math:`\ln(-\ln R)` is linear in :math:`\ln t`, data can be fitted
with a ruler on Weibull probability paper. The distribution became the
standard model of component life.

*Implementation:* :func:`engineeringkit.reliability.weibull_reliability`,
:func:`engineeringkit.reliability.weibull_cdf`,
:func:`engineeringkit.reliability.weibull_hazard` and
:func:`engineeringkit.reliability.weibull_mttf`. The example fits a small
life test on Weibull axes and compares the mean life with the gamma-function
formula.

*References:* W. Weibull, "A statistical distribution function of wide
applicability," *J. Appl. Mech.* 18 (1951) 293-297.

.. minigallery:: ../../examples/reliability/lifetimes/plot_03_weibull_distribution.py

1953 -- Epstein and Sobel's Exponential Life Testing
----------------------------------------------------

Benjamin Epstein and Milton Sobel worked out how to demonstrate a mean life
from a test that is stopped before every unit has failed. For exponential
lives the *total time on test*, summed over failed and surviving units,
carries all the information. The MTBF estimate is that time divided by the
number of failures :math:`r`, and twice the total time over the true MTBF is
chi-square with :math:`2r` degrees of freedom, which gives exact confidence
bounds. Their results became the basis of the US military's reliability
demonstration plans.

*Implementation:* :func:`engineeringkit.reliability.exponential_reliability`,
:func:`engineeringkit.reliability.exponential_mttf`,
:func:`engineeringkit.reliability.failure_rate_from_mtbf` and
:func:`engineeringkit.reliability.mtbf_confidence_bounds`. The example
repeats a failure-terminated test many times and checks the coverage of the
90% bounds.

*References:* B. Epstein and M. Sobel, "Life testing," *J. Amer. Statist.
Assoc.* 48 (1953) 486-502.

.. minigallery:: ../../examples/reliability/lifetimes/plot_04_epstein_sobel_life_testing.py

1956 -- Von Neumann's Redundancy
--------------------------------

John von Neumann asked how a computer, or a brain, can work when each of
its components is unreliable. Active parallel redundancy, :math:`1 - (1 -
r)^n`, survives if any copy survives. For logic, where a wrong output is as
bad as none, he used majority voting: three copies and a voter work with
probability :math:`3r^2 - 2r^3`. Voting on voters drives the error toward
zero, but only if each part is right more than half the time. The paper
founded the theory of fault-tolerant computing.

*Implementation:* :func:`engineeringkit.reliability.parallel_reliability`.
The example compares parallel and triple-modular redundancy and iterates
the majority vote.

*References:* J. von Neumann, "Probabilistic logics and the synthesis of
reliable organisms from unreliable components," in C. E. Shannon and J.
McCarthy (eds.), *Automata Studies* (Princeton University Press, 1956)
43-98.

.. minigallery:: ../../examples/reliability/systems/plot_01_von_neumann_redundancy.py

1957 -- The AGREE Report: Series Systems and MTBF
-------------------------------------------------

In the early 1950s much of the US forces' electronic equipment was out of
service at any one time, mostly because of vacuum-tube failures. The
Advisory Group on Reliability of Electronic Equipment (AGREE) recommended
that reliability be specified as a number in contracts, predicted from
parts counts and demonstrated by test. Its model was a series system of
parts with constant failure rates: the system's failure rate is the sum of
the parts' rates and its MTBF the reciprocal, so reliability falls with
every part added.

*Implementation:* :func:`engineeringkit.reliability.series_reliability`,
:func:`engineeringkit.reliability.exponential_hazard` and
:func:`engineeringkit.reliability.exponential_failure_probability`. The
example makes a parts-count prediction for a receiver and checks it
against the product of the parts' reliabilities.

*References:* Advisory Group on Reliability of Electronic Equipment,
*Reliability of Military Electronic Equipment* (Office of the Assistant
Secretary of Defense, Washington, 4 June 1957).

.. minigallery:: ../../examples/reliability/systems/plot_02_agree_series_systems.py

1958 -- Kaplan and Meier's Product-Limit Estimator
--------------------------------------------------

Field and clinical data are incomplete: some units are still working when
the data are collected and others are withdrawn for reasons other than
failure. Edward Kaplan and Paul Meier estimated the survival function by
multiplying the conditional survival :math:`1 - d_j/n_j` at each failure
time, where :math:`n_j` counts only the units still under observation. The
estimator uses every censored unit for as long as it was observed. Their
paper became one of the most cited in statistics.

*Implementation:* :func:`engineeringkit.reliability.kaplan_meier`, which
returns a :class:`~engineeringkit.reliability.KaplanMeierResult` with
Greenwood's standard errors. The example censors a Weibull sample at random
and recovers the true survival curve.

*References:* E. L. Kaplan and P. Meier, "Nonparametric estimation from
incomplete observations," *J. Amer. Statist. Assoc.* 53 (1958) 457-481.

.. minigallery:: ../../examples/reliability/lifetimes/plot_05_kaplan_meier.py

1961 -- Birnbaum, Esary and Saunders's k-out-of-n Systems
---------------------------------------------------------

Zygmunt Birnbaum, James Esary and Sam Saunders put system reliability on a
general footing. A *structure function* says which combinations of working
components keep the system working, and for independent components the
system's reliability follows from it. The :math:`k`-out-of-:math:`n`
system, which works if at least :math:`k` of its :math:`n` identical
components do, runs from parallel (:math:`k = 1`) to series
(:math:`k = n`), and its reliability is a binomial tail.

*Implementation:* :func:`engineeringkit.reliability.k_of_n_reliability`.
The example checks the binomial formula by enumerating component states and
computes the mean life of each architecture.

*References:* Z. W. Birnbaum, J. D. Esary and S. C. Saunders,
"Multi-component systems and structures and their reliability,"
*Technometrics* 3 (1961) 55-77.

.. minigallery:: ../../examples/reliability/systems/plot_03_k_out_of_n.py

1961 -- Watson's Fault Tree
---------------------------

The US Air Force asked Bell Telephone Laboratories whether the Minuteman
launch-control system could launch a missile without a valid order. H. A.
Watson answered with a new diagram that starts from the undesired top event
and works down through AND and OR gates to basic failures with known
probabilities. Boeing extended the method to the whole Minuteman system,
and David Haasl presented it at the first System Safety Symposium in 1965.
Fault trees became standard in aerospace and nuclear safety. An event that
feeds more than one branch makes the gate inputs dependent, so the tree
must then be evaluated over the states of its basic events.

*Implementation:* :func:`engineeringkit.reliability.fault_tree_probability`.
The example evaluates a tree with a repeated event exactly, checks it by
Monte Carlo and shows the error of combining the gates as if independent.

*References:* H. A. Watson, *Launch Control Safety Study* (Bell Telephone
Laboratories, Murray Hill, 1961); D. F. Haasl, "Advanced concepts in fault
tree analysis," *System Safety Symposium* (Seattle, 1965).

.. minigallery:: ../../examples/reliability/risk/plot_01_watson_fault_tree.py

1964 -- Duane's Reliability Growth
----------------------------------

James Duane, at General Electric, plotted the failure records of several
development programmes, from hydraulic devices to aircraft generators. On
log-log paper the cumulative MTBF, total test time over total failures,
rose along a straight line, :math:`\theta_c = K t^\alpha`. The
instantaneous MTBF is :math:`\theta_c/(1 - \alpha)`, and the line tells a
programme how much more testing a reliability target needs. Larry Crow
later recast the model as a non-homogeneous Poisson process (the AMSAA
model).

*Implementation:* :func:`engineeringkit.reliability.reliability_growth_duane`.
The example simulates a power-law failure process, fits Duane's line and
predicts the test time needed to reach a target MTBF.

*References:* J. T. Duane, "Learning curve approach to reliability
monitoring," *IEEE Trans. Aerospace* 2 (1964) 563-566.

.. minigallery:: ../../examples/reliability/growth/plot_01_duane_growth.py

1965 -- Barlow and Proschan's Age Replacement
---------------------------------------------

Richard Barlow and Frank Proschan's *Mathematical Theory of Reliability*
used renewal theory to decide when a part should be replaced before it
fails. Under an age-replacement policy a part is replaced at failure or at
age :math:`T`, whichever comes first. The long-run cost rate is the
expected cost per cycle over the expected cycle length, and it has a finite
minimum only if the hazard rises and a failure costs more than a planned
replacement. With a constant hazard a new part is no better than the old
one, so planned replacement never pays.

*Implementation:* :func:`engineeringkit.reliability.age_replacement_cost_rate`
and :func:`engineeringkit.reliability.renewal_interval_cost_rate`. The
example finds the optimum replacement age of a Weibull part, checks the
optimality condition and confirms the cost rate by simulation.

*References:* R. E. Barlow and F. Proschan, *Mathematical Theory of
Reliability* (Wiley, New York, 1965), ch. 4.

.. minigallery:: ../../examples/reliability/maintenance/plot_01_barlow_proschan_replacement.py

1965 -- Markov Models of Repairable Systems
-------------------------------------------

A unit that fails and is repaired at constant rates is a continuous-time
Markov chain, and so is any system of such units. Barlow and Proschan gave
the availability of a single unit, which decays from one to the steady
state :math:`\mu/(\lambda + \mu)` with time constant
:math:`1/(\lambda + \mu)`, and treated larger systems through their
transition-rate matrices. The state probabilities are the matrix
exponential of that matrix applied to the initial state, and setting their
derivatives to zero gives the long-run balance equations.

*Implementation:* :func:`engineeringkit.reliability.markov_availability`,
:func:`engineeringkit.reliability.markov_state_probabilities`,
:func:`engineeringkit.reliability.availability`,
:func:`engineeringkit.reliability.unavailability`,
:func:`engineeringkit.reliability.repair_rate` and
:func:`engineeringkit.reliability.maintainability`. The example checks the
closed form against the matrix exponential and solves two pumps with one
repair crew.

*References:* R. E. Barlow and F. Proschan, *Mathematical Theory of
Reliability* (Wiley, New York, 1965), ch. 5.

.. minigallery:: ../../examples/reliability/maintenance/plot_02_markov_availability.py

1969 -- Birnbaum's Component Importance
---------------------------------------

Which component should be improved first? Z. W. Birnbaum defined a
component's importance as the rate at which the system's reliability
changes with that component's reliability,
:math:`I_i = \partial R_s/\partial R_i`. It equals the probability that the
component is critical, that the system works if it works and fails if it
fails. In a series system the importance is the product of the other
components' reliabilities, so the weakest part matters most.

*Implementation:* :func:`engineeringkit.reliability.birnbaum_importance_series`.
The example checks the derivative numerically and shows that improving the
most important component gains the most.

*References:* Z. W. Birnbaum, "On the importance of different components
in a multicomponent system," in P. R. Krishnaiah (ed.), *Multivariate
Analysis II* (Academic Press, New York, 1969) 581-592.

.. minigallery:: ../../examples/reliability/systems/plot_04_birnbaum_importance.py

1974 -- Hasofer and Lind's Reliability Index
--------------------------------------------

C. Allin Cornell's 1969 reliability index, the mean of the safety margin
over its standard deviation, depends on how the failure condition is
written. Abraham Hasofer and Niels Lind removed the ambiguity by
transforming to standard normal variables and measuring the shortest
distance from the origin to the failure surface. The nearest point, the
*design point*, is the most likely way to fail, and
:math:`\Phi(-\beta)` approximates the failure probability. The index is the
basis of the first-order reliability method used to calibrate modern
structural codes.

*Implementation:* :func:`engineeringkit.reliability.hasofer_lind_index`,
which returns a :class:`~engineeringkit.reliability.HasoferLindResult`,
with :func:`engineeringkit.reliability.failure_probability_from_index` and
:func:`engineeringkit.reliability.reliability_index_from_probability`. The
example shows that Cornell's index changes with the formulation and
Hasofer and Lind's does not, and checks the design point by Monte Carlo.

*References:* A. M. Hasofer and N. C. Lind, "Exact and invariant
second-moment code format," *J. Eng. Mech. Div. ASCE* 100 (1974) 111-121.

.. minigallery:: ../../examples/reliability/structural/plot_02_hasofer_lind.py

1975 -- WASH-1400: Event Trees
------------------------------

The Reactor Safety Study, led by Norman Rasmussen for the US Nuclear
Regulatory Commission, was the first attempt to put a number on the risk
of a nuclear power plant. It worked forward from each initiating event,
such as a pipe break, and asked in turn whether each safety system
worked. Every path through the resulting event tree is an accident
sequence whose frequency is the initiator's frequency times its branch
probabilities, and fault trees supplied the branch probabilities. The
study's uncertainty bounds were criticised, but the 1979 accident at Three
Mile Island followed a small-break sequence of the kind it had ranked high,
and probabilistic risk assessment became regulatory practice.

*Implementation:* :func:`engineeringkit.reliability.event_tree_probability`,
which returns an :class:`~engineeringkit.reliability.EventTreeResult`. The
example quantifies a loss-of-coolant tree and checks the core damage
frequency against inclusion and exclusion.

*References:* U.S. Nuclear Regulatory Commission, *Reactor Safety Study: An
Assessment of Accident Risks in U.S. Commercial Nuclear Power Plants*
(WASH-1400, NUREG-75/014, 1975).

.. minigallery:: ../../examples/reliability/risk/plot_02_wash1400_event_tree.py

1998 -- IEC 61508 and Safety Integrity Levels
---------------------------------------------

Programmable electronics took over process trips, railway signalling and
machine guards, and there was no agreed way to say how safe a
software-based protection system was. The International Electrotechnical
Commission's standard IEC 61508 grades each safety function by the risk
reduction it must achieve, in four safety integrity levels. For a function
that acts on rare demands the measure is its average probability of failure
on demand. A dangerous fault that no self-test detects stays hidden until
the next proof test, so a single channel averages
:math:`\lambda_{DU} T/2` over a test interval :math:`T`.

*Implementation:* :func:`engineeringkit.reliability.safety_integrity_pfd`
and :func:`engineeringkit.reliability.safety_integrity_level`. The example
checks the approximation against the exact time average and shows how
common-cause failures limit a redundant pair.

*References:* IEC 61508, *Functional safety of electrical/electronic/
programmable electronic safety-related systems* (International
Electrotechnical Commission, Geneva, 1998-2000; 2nd ed. 2010), part 6,
annex B.

.. minigallery:: ../../examples/reliability/risk/plot_03_iec61508_sil.py

See Also
--------

- :doc:`/api/reliability`
- :doc:`/guides/reliability`
- :doc:`/history/materials_breakthroughs`
- :doc:`/history/structural_engineering_breakthroughs`

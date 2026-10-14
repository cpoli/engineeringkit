Breakthroughs in Reliability
============================

.. include:: /_generated/nav/reliability.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

Survival analysis: failures become distributions
------------------------------------------------

Reliability is the probability of surviving beyond a time. A constant hazard gives the exponential model, while Weibull shape permits increasing or decreasing hazard.

.. math::

    R(t)=e^{-(t/\eta)^\beta}

Use :py:obj:`engineeringkit.reliability.weibull_reliability` to explore the model.
Time and scale must share units; fitting a distribution and checking censoring are outside these algebraic helpers.

Redundancy: system architecture changes survival
------------------------------------------------

Series systems require every component to survive; parallel systems need at least one. Their simple product expressions make independence a critical modeling assumption.

.. math::

    R_s=\prod_i R_i,\qquad R_p=1-\prod_i(1-R_i)

Use :py:obj:`engineeringkit.reliability.parallel_reliability` to explore the model.
Common-cause failures violate independence and can sharply reduce the benefit of redundancy.

Repair: availability answers a different question
-------------------------------------------------

A repairable system can fail repeatedly yet be available for most of its lifetime. Steady availability combines mean operating time and mean repair time.

.. math::

    A=MTBF/(MTBF+MTTR)

Use :py:obj:`engineeringkit.reliability.availability` to explore the model.
This long-run alternating-renewal ratio is distinct from the probability of no failure during a mission.

Explore and verify
------------------

:doc:`/guides/reliability` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/reliability/plot_reliability_models.py ../../examples/reliability/plot_reliability_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://www.itl.nist.gov/div898/handbook/apr/apr_d.htm>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

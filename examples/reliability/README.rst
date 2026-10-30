Examples
========

.. include:: /_generated/nav/reliability.rst

This gallery walks through every public feature of
``engineeringkit.reliability``, one breakthrough per script: lifetime
distributions and censored data, Weibull fitting and Bayesian failure
rates, accelerated life testing, series, parallel, k-out-of-n and network
systems, stress-strength, first-order and Monte Carlo structural
reliability, FMEA, fault and event trees, availability and replacement,
and reliability growth.

:doc:`Units, assumptions and numerical validation </guides/reliability>`
gives the validity limits. Each script is self-contained and can be run
directly with ``python examples/reliability/<section>/<script>.py``.

Sections
--------

- **lifetimes** -- Halley's life table, Poisson's rare failures, Weibull's distribution, Epstein and Sobel's life testing, the Kaplan-Meier estimator and McAlister's lognormal law.
- **life_data** -- Johnson's median ranks, Cohen's censored maximum likelihood and Raiffa and Schlaifer's conjugate prior for a failure rate.
- **acceleration** -- Arrhenius's law and Eyring's activated complex.
- **systems** -- Von Neumann's redundancy, the AGREE report's series systems, k-out-of-n systems, Birnbaum's importance and Esary and Proschan's path and cut bounds.
- **structural** -- Freudenthal's stress-strength interference, Hasofer and Lind's reliability index and Metropolis and Ulam's Monte Carlo method.
- **risk** -- MIL-P-1629's FMEA, Watson's fault tree, WASH-1400's event trees and IEC 61508's safety integrity levels.
- **maintenance** -- Barlow and Proschan's age replacement and Markov availability.
- **growth** -- Duane's reliability growth.

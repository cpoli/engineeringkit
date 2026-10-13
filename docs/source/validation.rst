Numerical validation policy
===========================

Validation is layered. Many public scalar helpers reject nonfinite inputs and
impossible physical ranges before division, logarithms or square roots.
Geometric/distribution/matrix routines validate ordering, shape or integer
indices where needed. Domain guides distinguish enforced bounds from assumptions
that cannot be inferred from the arguments. Passing validation does not certify
that a physical model applies to a particular apparatus.

Independent benchmark matrix
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Domain
     - Independent evidence
   * - Structural
     - Beam equilibrium/curvature, support values, bar reaction balance, energy derivative
   * - Electrical
     - Closed-form transient current, charge integral, RC energy accounting, resonance
   * - Controls
     - Analytic scalar LQR/Kalman, matrix Riccati residuals and shape/weight failures
   * - Thermodynamics
     - Integrated P dV, T delta_S, four-state cycle balances and COP identities
   * - Fluids
     - Integrated parabolic profile, laminar Darcy equivalence, implicit Colebrook residual
   * - Heat transfer
     - Interface Q*R drops, thin-wall limit, equal-temperature and exchanger limits
   * - Materials
     - Isotropic modulus identity, elastic strain-energy triangle, fracture inverse, composite bounds
   * - Mechanical
     - Undamped energy, damping dissipation, critical solution, gear power, forced equilibrium
   * - Aerodynamics
     - Shock mass/momentum/energy balances, sonic and infinite-wing limits
   * - Signals
     - Parseval, periodogram integral, cutoff gain, sinc nodes, ADC error bounds
   * - Communications
     - Analytic FM phase, unsigned bit handling, linear/log link budgets and noise limits
   * - Reliability
     - Survival integral, SciPy distributions, exhaustive independent system states, sensitivities

Tests use analytical values, invariants, independent quadrature and established
SciPy distributions. Where a governing-equation residual is checked, it verifies
the numerical solve but does not independently prove the equation's physics.
Most assertions use tight tolerances for algebraic identities; integrations allow
appropriate quadrature/discretization error. Source coverage excludes tests and
is a diagnostic, not a proof of correctness. Coverage gaps and untested branches
are retained in the release report.

Behavior retained in v0.1
-------------------------

SeriesRLC returns current, correcting the previous charge/current mismatch.
FM integration starts at zero modulation phase and supports nonuniform grids.
Unsigned BPSK/QPSK input is converted before subtraction. Physical range errors
raise ValueError instead of propagating NaNs or divisions by zero. Colebrook
transition inputs are rejected. Rare-failure PFD and empirical Oswald estimates
now enforce their documented implementation ranges. NumPy 1.24 remains supported
by using SciPy trapezoid for stress-strain integration.

Coverage gate
-------------

The release suite covers all executable production statements. Coverage excludes
only test code; production validation, plotting, constants and unit helpers are
included. The report's fail_under=100 gate fails on any missed statement. This
is statement coverage; branch coverage is not configured. An executed statement
can still contain an error or rely on an invalid model assumption, so analytical
benchmarks and the domain guides remain essential.

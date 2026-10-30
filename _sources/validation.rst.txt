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
     - Beam equilibrium/curvature, support values, bar reaction balance, energy derivative, three-moment equation, moment distribution and frame stiffness against one another, method of joints against direct stiffness, mechanism work minimum, Engesser residual, Macaulay beams against handbook cases and superposition, Gauss-quadrature stress resultants of the flexure formula, Jourawski stress integrating to the shear force, Lamé equilibrium residual, stress-tensor invariants, rosette round trip, unit-load against stiffness deflection, slope-deflection against moment distribution, cable arc length by quadrature
   * - Electrical
     - Closed-form transient current, charge integral, RC and RL energy accounting, resonance, RLC step against SciPy, mesh against nodal analysis, diode load line against Lambert W, Newton-Raphson against Gauss-Seidel load flow, swing-curve energy and equal areas
   * - Controls
     - Analytic scalar LQR/Kalman, matrix Riccati residuals and shape/weight failures
   * - Thermodynamics
     - Integrated P dV, T delta_S, four-state cycle balances and COP identities, IF97 verification tables (regions 1 to 5) and Gibbs relations, van der Waals critical inflection and isenthalpic throttle, Curzon-Ahlborn power maximum, Rankine and Brayton plant energy balances and Çengel and Boles plants, Ericsson limit, Stirling and dual four-process balances, Peng-Robinson critical inflection, equal-area rule, Clapeyron slope and R-134a tables, Maxwell relation and Joule-Thomson throttle on IF97 steam, Dalton mixing entropy against isothermal expansions, psychrometric mass and energy balances, Hess heating values, flame-temperature energy balance and JANAF enthalpies
   * - Fluids
     - Integrated parabolic profile and wall pressure, laminar Darcy equivalence, Colebrook residual and limits, Hardy Cross against Newton, Pi-group null space
   * - Heat transfer
     - Interface Q*R drops, thin-wall limit, exchanger limits, Planck integral, plane-wall series and finite-difference energy, Blasius-Pohlhausen similarity solution, Graetz eigenvalue, exact crossflow series, effectiveness-NTU against F-LMTD, cylinder and sphere lumped and fixed-surface limits and energy balance, 2-D Fourier series and second-order convergence, shape-factor identities, two- and three-surface radiation networks
   * - Materials
     - Isotropic modulus identities, elastic strain-energy triangle, offset yield, Griffith-Irwin agreement, Paris life against closed form and quadrature, composite bounds, strain-life transition, Fick's second law, Callister and Shigley textbook values, Neuber's product, Halpin-Tsai limits, laminate ABD closed forms, Weibull size effect
   * - Mechanical
     - Undamped energy, damping dissipation, critical solution, gear power, forced equilibrium, simulated forced steady state, modal superposition against direct integration, double-pendulum energy, Rayleigh-Dunkerley bracket, Castigliano spring rate, integrated pressure cone, balancing residuals, finite-difference linkage and cam derivatives
   * - Aerodynamics
     - Normal and oblique shock mass/momentum/energy balances, sonic and infinite-wing limits, Prandtl-Meyer and Fanno ODE quadrature, Rayleigh-line conservation, exact elliptic lifting line, panel method against the exact cylinder and Kutta-Joukowski lift, built-up moment slope at the neutral point, integrated take-off roll and time-stepped propeller cruise
   * - Signals
     - Parseval, periodogram integral, cutoff gain, sinc nodes, ADC error bounds, FFT and Welch against NumPy and SciPy, Gibbs limit, Haar energy, Wiener and LMS system identification, Chebyshev and elliptic ripple bounds and orders, bilinear transform and prewarping, Parks-McClellan alternation, overlap-add against direct convolution, interpolation and decimation, Hilbert transform, DCT orthonormality, Levinson-Durbin against a Toeplitz solve, cepstral echo
   * - Communications
     - Analytic FM phase, unsigned bit handling, linear/log link budgets and noise limits
   * - Reliability
     - Survival integral, SciPy distributions, exhaustive independent system states, sensitivities, Monte Carlo fault trees, interference and Kaplan-Meier, Markov matrix exponential, beta-median ranks, SciPy and censored-likelihood Weibull fits, flood-fill network states, simulated standby, Monte Carlo against FORM, numerical gamma-Poisson posterior

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

Scope and remaining limitations
===============================

Version 0.1 strengthens physical input contracts, adds public scalar input/output
units and detailed validity guides, and independently benchmarks all twelve
domains. This does not establish that every original roadmap method is covered
or that these ideal models represent a particular real system.

Unresolved work includes wider benchmark coverage for legacy helpers, full-
project static typing, comprehensive historical primary references, property-
table thermodynamics, high-order digital filters in second-order sections,
conditioning diagnostics for near-singular trusses, and robust multi-crossover
stability margins. Some low-level SciPy wrappers delegate shape/range errors to
SciPy. Existing result objects and complex multidimensional helper APIs do not
yet all have full NumPy-style parameter docstrings. Raw numbers do not carry
units, so dimensional mistakes remain possible despite explicit documentation.

Model-specific unenforced assumptions (laminar development, Biot limits,
small-angle attached flow, far-field antennas, independent failures, linear
elasticity) are described in the domain guides. Probabilistic models are not
empirical reliability forecasts without calibrated data. Scalar conditional
parameters such as gamma, aspect ratio and distribution shape are not uniformly
vectorized; use vector arrays only where the model's operations support them.

Local verification records exactly which checks ran. CI configuration is not
evidence that the remote multi-version matrix has already passed. Releases are
published to PyPI only from a ``vX.Y.Z`` tag (``.github/workflows/release.yml``).

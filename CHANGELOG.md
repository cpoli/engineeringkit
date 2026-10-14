# Changelog

The release workflow takes a release's notes from its `## [X.Y.Z] - date`
section, so keep that heading format.

## [0.1.0] - Unreleased

- Add release infrastructure: `CITATION.cff`, a tag-driven release workflow
  (PyPI trusted publishing, GitHub Release, docs rebuild), a docs workflow
  that deploys `gh-pages` on every push to `main`, and a CI matrix over
  Python 3.10-3.14 on Linux and macOS.
- Reset the public release version to 0.1.0 at the maintainer's request.
- Reach 100% production statement coverage: 398 tests, zero missed statements.
- Enforce 100% statement coverage in the configured coverage report and CI.
- Add analytical/reference tests for all previously uncovered numerical helpers,
  result summaries, plots, animations, and physical-input error paths.
- Fix the SDOF animation interval argument and integer-input digital filtering.
- Remove unreachable duplicate validation guards without changing input contracts.
- Retain all twelve domains, validated models and detailed documentation.
- Add course material: four cross-domain tutorials, an exercise set with
  worked solutions for each of the twelve domains, and a quickstart
  notebook (Colab badge, served in JupyterLite). The docs build runs every
  tutorial and solution, and the test suite checks the numbers they print.

## Unpublished development milestones

Earlier 0.2–0.5 labels were development iterations, not external publications.
The initial public release is now 0.1.0. Full historical development notes are
available in earlier source versions; all implementation improvements are retained.

### Preserved development notes

### 0.5.0

- Add physical-domain/nonfinite input checks across the twelve engineering domains.
- Correct SeriesRLC charge/current mismatch and FM integration phase offset.
- Fix unsigned BPSK/QPSK mapping and validate binary input and matching shapes.
- Preserve tiny failure probabilities, large Poisson counts, and near-singular cycle/exchanger limits.
- Add Riccati matrix/weight, truss DOF/connectivity and stress-strain data validation.
- Add analytical, conservation, quadrature, distribution and invalid-domain benchmarks.
- Document scalar public API input/output units and twelve detailed model-validity guides.
- Add warnings-as-errors test policy and production-only coverage configuration.

### 0.4.0

- Implement shared Euler/RK4/adaptive integration with inspectable histories.
- Implement reusable state plots and displacement animations.
- Validate passive oscillator inputs and solver success; expose history/summary.
- Fix PID zero-integral cancellation, Routh zero rows and marginal stability.
- Validate root-locus polynomial degrees; fix complex SNR and mid-rise ADC levels.
- Handle zero m*L fin efficiency through its analytic limit.
- Format source, repair lint checks, include py.typed and expand Sphinx navigation/API.
- Add analytical, conservation, convergence, invalid-input and rendered-media tests.



### 0.3.0
- Completed all Phase 2 and Phase 3 engineering domains.
- Added heat transfer, materials, signal processing, mechanical engineering, aerodynamics, communications, and reliability implementations.
- Added domain-local numerical tests, documentation progression pages, and an all-domain executable gallery.


### 0.2.0
- Expanded all five Phase-1 domains with additional educational methods.
- Added direct-stiffness planar trusses, pole placement and Kalman estimation.
- Added RC/RL network helpers, Rankine/exergy relations, fluid similarity numbers and water hammer.
- Moved domain tests alongside subpackages and added history pages plus a runnable gallery.

### 0.1.0
- Initial educational foundation and five Phase-1 domains.

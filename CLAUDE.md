# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

EngineeringKit is an educational scientific-Python toolkit (sibling project to MathematicsKit, same
design philosophy): engineering concepts made computable, inspectable, visualizable, and understandable.
Public APIs take and return plain numbers/arrays in SI units — no unit-carrying dimensional types. NumPy/SciPy
are used for mature numerics; otherwise formulas/algorithms are hand-written where seeing the mechanics is
the pedagogical point.

Twelve domains live as sibling packages under `engineeringkit/`: `aerodynamics`, `communications`,
`control_systems`, `electrical_engineering`, `fluid_engineering`, `heat_transfer`, `materials`,
`mechanical_engineering`, `reliability`, `signal_processing`, `structural_engineering`, `thermodynamics`.
Plus shared infrastructure: `_validation.py` (input checks), `units/` (explicit SI conversion), `constants/`,
`integrators/` (ODE solvers), `visualizers/` (solver-independent plotting/animation helpers).

## Commands

```bash
# setup
python -m pip install -e '.[dev]'

# tests (MPLBACKEND=Agg avoids GUI backend issues/hangs from matplotlib)
MPLBACKEND=Agg python -m pytest
MPLBACKEND=Agg python -m pytest engineeringkit/structural_engineering/tests/test_structural.py
MPLBACKEND=Agg python -m pytest engineeringkit/structural_engineering/tests/test_structural.py -k test_name

# coverage (CI enforces fail_under=100 on production statements)
MPLBACKEND=Agg python -m pytest --cov=engineeringkit --cov-report=term-missing

# lint / format / types
python -m ruff check .
python -m ruff format --check .
python -m mypy   # whole package, settings in pyproject.toml; blocking in CI

# docs
python -m sphinx -W --keep-going -b html docs/source docs/_build/html

# build
python -m build
```

CI (`.github/workflows/ci.yml`) runs all of the above across Python 3.10-3.13; treat it as the
authoritative gate before considering a change done.

## Architecture conventions (apply across all twelve domains)

- **Domain layout is uniform.** Each domain is a single-module package (`engineeringkit/<domain>/__init__.py`)
  containing free functions for algebraic formulas plus `@dataclass(frozen=True)` model/result classes for
  stateful solvers (e.g. `SimplySupportedBeam`/`BeamResult`, `PlanarTruss`/`TrussResult`). Newer domains
  additionally use a `tests/` subpackage; older ones keep `test_<domain>.py`, `test_benchmarks.py`,
  `test_complete_api.py` as sibling files in the domain directory — follow whichever layout a given domain
  already uses.
- **Result objects share a shape.** Solver outputs are frozen dataclasses exposing `.summary()` (dict of
  scalar diagnostics) and `.plot()` (returns `(fig, ax)` or `(fig, axs)`; imports matplotlib lazily inside
  the method).
- **Validation is explicit and centralized in `_validation.py`.** Every public function validates its own
  inputs at the top via `finite`, `positive`, `nonnegative`, `interval`, or `integer` — imported aliased with
  a leading underscore (`from engineeringkit._validation import positive as _positive`). Don't rely on
  downstream NumPy/SciPy errors for input contracts; validate explicitly like the surrounding code does.
- **Docstrings are NumPy-style with a fixed Notes block.** Every public function's docstring documents each
  parameter's SI unit, the return unit, and ends with a Notes paragraph pointing to
  `:doc:`/guides/<domain>`` and stating inputs are raw SI numbers with no dimensional metadata. Match this
  pattern exactly when adding new public functions — the Sphinx docs and `docs/source/guides/` pages depend
  on it.
- **Three-tier test convention per domain**: unit/analytical tests (`test_<domain>.py`), independent
  conservation/analytical/distribution benchmark checks (`test_benchmarks.py`), and a sweep over the full
  public API surface (`test_complete_api.py`). When adding a model, add tests in the tier(s) that match —
  per `CONTRIBUTING.md`, every new model needs an analytical or benchmark test, not just a unit test.
- **100% statement coverage is enforced**, scoped to production code only (`tool.coverage.run: omit =
  ["*/test_*.py", "*/tests/*"]` in `pyproject.toml`). No `# pragma: no cover` escapes are used anywhere in
  this codebase — if new code can't be covered by a real test, reconsider the code rather than suppressing
  coverage.
- **`filterwarnings=["error"]`** is set in `pyproject.toml`, so any warning (including matplotlib/NumPy
  deprecation warnings) fails the test suite. Fix the root cause rather than filtering it.

## Documentation

Sphinx docs under `docs/source/` follow MathematicsKit's domain-first architecture (see `docs/README.md` for
the full map): a single registry, `SUBPACKAGES` in `docs/source/conf.py`, is the source of truth for domain
names/categories/descriptions and drives generated hub pages, cross-links, and API/history/example groupings.
When adding or renaming a domain, update that registry rather than hand-writing new hub/nav pages. Executed
example galleries live in `examples/<domain>/plot_<domain>_models.py` (+ `plot_<domain>_benchmark.py`) with a
`README.rst` header; Sphinx-Gallery executes them during the doc build, so example code must actually run.
Tutorials (`docs/source/tutorials/`) and per-domain exercises (`docs/source/exercises/`) are `.. plot::` blocks
whose printed output is written as `# ` comments under each print; `engineeringkit/tests/test_course_material.py`
re-runs them and checks those numbers, so update the comments when a model's output changes.

# Output Guidelines
- **Be Concise:** Provide direct code and answers first. Omit setup text, fluff, and conversational responses.
- **Code Generation:** Output only updated code blocks, functions, or unified diffs. Never output full files unless instructed.
- **Explanations:** Limit inline code comments. Explain high-level logic in 1–2 bullet points after the code block.
- **Format:** Use bulleted lists and tables for comparisons instead of paragraph blocks.
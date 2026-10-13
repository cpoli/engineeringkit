# Contributing

Install `pip install -e .[dev]`, add analytical/benchmark tests with every model, run `pytest` and `ruff check .`. Public APIs use SI units and NumPy-style documentation.

## Type checking

`python -m mypy` with no arguments checks the whole package, and a failure
fails CI. Library code is checked with `check_untyped_defs`,
`strict_equality` and `warn_unreachable`, so the bodies of unannotated
functions are checked too. Tests are checked at the signature level only.
All settings are in `[tool.mypy]` in `pyproject.toml`.

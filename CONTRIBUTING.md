# Contributing

Install `pip install -e .[dev]`, add analytical/benchmark tests with every model, run `pytest` and `ruff check .`. Public APIs use SI units and NumPy-style documentation.

## Type checking

`python -m mypy` with no arguments checks the whole package, and a failure
fails CI. Library code is checked with `check_untyped_defs`,
`strict_equality` and `warn_unreachable`, so the bodies of unannotated
functions are checked too. Tests are checked at the signature level only.
All settings are in `[tool.mypy]` in `pyproject.toml`.

## Range of validity

A model that still returns a number outside the range its docstring states
(a correlation used beyond the Reynolds or Prandtl numbers it was fitted
to) calls `validity(model, name, value, low, high)` from
`engineeringkit/_validation.py` after its `ValueError` checks. That issues
an `engineeringkit.ValidityWarning` (a `UserWarning`) and leaves the
result unchanged. State the range in the parameter or Notes text, say that
it warns, and test both sides of it: `pytest.warns(ValidityWarning)`
outside, and no warning at the ends. Inputs with no physical answer still
raise `ValueError`.

Since pytest runs with `filterwarnings = ["error"]`, a test or gallery
example that goes outside a range on purpose must say so with
`pytest.warns` or `warnings.catch_warnings`.

To silence the warning in your own code:

```python
import warnings

import engineeringkit as ek

warnings.filterwarnings("ignore", category=ek.ValidityWarning)  # everywhere

with warnings.catch_warnings():  # or only around one call
    warnings.simplefilter("ignore", ek.ValidityWarning)
    Nu = ek.heat_transfer.dittus_boelter_nusselt(500.0, 0.7)
```

## Array inputs

A public function whose parameters are all scalars (documented as
`float`, `float or array_like`, `complex`, `int`, `str` or `bool`) and
whose result is a number, a tuple of numbers or a result object of numbers
is a *scalar function*. Its `float` and `complex` parameters accept NumPy
arrays and broadcast elementwise:

- a scalar call returns a number (a Python or NumPy scalar), an array call
  returns an `ndarray` of the broadcast shape;
- a tuple result becomes a tuple of arrays, not an array of tuples, and a
  result object holds arrays in its fields;
- validation checks every element (`np.any(...)`, never a bare `if a < b`),
  and a branch on a value (`if Cr == 1`) becomes `np.where`;
- `int`, `str` and `bool` parameters (counts, orders, modes, flags) stay
  scalars.

A parameter may be **scalar-only** when each value needs its own root
find, quadrature or integration, or when the result describes one object
(a filter, a rotor, an absorber design). Say so in the first Notes
paragraph (``` ``Bi`` is scalar-only, since ...; ``Fo`` and ``x_star``
broadcast. ```) and name the parameters in `SCALAR_ONLY` in
`engineeringkit/tests/test_array_inputs.py`.

That test calls every scalar function once with 1-D arrays built from
scalar cases taken from the domain's own tests and compares the result with
the scalar calls. A new scalar function needs an entry in its `CASES`
table as soon as a domain test calls it. Functions added by D items still
in progress on 2026-10-23 are listed in `PENDING` (strict xfail): the item
that added each one vectorizes it or moves it to `SCALAR_ONLY`.

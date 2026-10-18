# engineeringkit

[![PyPI version](https://badgen.net/pypi/v/engineeringkit)](https://pypi.org/project/engineeringkit/) [![Python versions](https://badgen.net/pypi/python/engineeringkit)](https://pypi.org/project/engineeringkit/) [![License](https://img.shields.io/github/license/cpoli/engineeringkit)](https://github.com/cpoli/engineeringkit/blob/main/LICENSE) [![CI](https://github.com/cpoli/engineeringkit/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/cpoli/engineeringkit/actions/workflows/ci.yml) [![Docs](https://badgen.net/badge/docs/cpoli.github.io%2Fengineeringkit/blue)](https://cpoli.github.io/engineeringkit/) [![Ruff](https://badgen.net/badge/code%20style/ruff/purple)](https://github.com/astral-sh/ruff) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23452015.svg)](https://doi.org/10.5281/zenodo.23452015)

**Engineering concepts made computable → inspectable → visualizable → understandable.**

EngineeringKit is an educational scientific-Python toolkit following the design philosophy of MathematicsKit: explicit intermediate results, NumPy/SciPy for mature numerics, and hand-written formulas/algorithms where seeing the mechanics is the lesson. SI units are used by public APIs.

## Install
```bash
pip install engineeringkit
```

## Beam
```python
from engineeringkit.structural_engineering import SimplySupportedBeam

beam = SimplySupportedBeam(5.0, 200e9, 8e-6)
r = beam.solve_point_load(10_000)
print(r.summary())
r.plot()
```

## PID
```python
from scipy import signal
from engineeringkit.control_systems import PID

plant = signal.TransferFunction([1], [1, 2, 1])
r = PID(3, 1, 0.2).simulate(plant)
r.plot()
```

## RLC
```python
from engineeringkit.electrical_engineering import SeriesRLC

r = SeriesRLC(10, 0.1, 100e-6).step_response()
r.plot()
```

## Thermodynamics
```python
from engineeringkit.thermodynamics import OttoCycle

print(OttoCycle(10).summary())
```

## Pipe flow
```python
from engineeringkit.fluid_engineering import reynolds_number, colebrook_friction_factor

Re = reynolds_number(998, 2, 0.05, 1.002e-3)
f = colebrook_friction_factor(Re, 1e-4)
```

## Release scope (v0.1)

All twelve domains include executable educational models. Version 0.1 adds
physical input validation, independent conservation/analytical/distribution
benchmarks, public scalar API units, and twelve detailed model guides. It also
corrects RLC current response, FM phase integration and unsigned modulation,
and improves numerical limiting cases. See `docs/source/guides/` for model
assumptions and `VERIFICATION.md` for the actual test/build evidence.

Documentation: <https://cpoli.github.io/engineeringkit/>. To cite
engineeringkit, use `CITATION.cff` (GitHub's "Cite this repository").

## Inspect a numerical method

```python
import numpy as np
from engineeringkit.integrators import integrate_ode

result = integrate_ode(lambda t, y: -y, [1.0], np.linspace(0, 2, 101), "rk4")
print(result.summary())
fig, ax = result.plot(["normalized state"])
```

## Verify locally

```bash
python -m pip install -e '.[dev]'
MPLBACKEND=Agg python -m pytest
python -m ruff check .
python -m ruff format --check .
python -m sphinx -W --keep-going -b html docs/source docs/_build/html
python -m build
```

## Coverage requirement

Run the suite with production statement coverage:

```bash
MPLBACKEND=Agg python -m pytest --cov=engineeringkit --cov-report=term-missing
```

The report enforces 100% statement coverage. Test modules are excluded from the
coverage denominator; all production modules, private validation helpers,
constants, units, plotting and integration code are included. No production
statements are marked with coverage-ignore pragmas. This measures statement
execution, not 100% branch coverage or proof of scientific correctness.

## Documentation navigation

Documentation now follows MathematicsKit's domain-first architecture: Subpackages,
History, Examples and API. Each domain hub links all three learning paths; grouped
cards and cross-links are generated from one registry in `docs/source/conf.py`.
Executed galleries come from `examples/<domain>/`, and detailed SI/assumption
guides remain linked from the galleries and domain hubs. See `docs/README.md`.

The MathematicsKit-style Sphinx navigation connects twelve domain hubs to
History, Examples and API pages. History pages explain key equations and limits;
24 plotted gallery scripts include twelve explicit benchmark demonstrations.

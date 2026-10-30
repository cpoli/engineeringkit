# engineeringkit

[![PyPI version](https://badgen.net/pypi/v/engineeringkit)](https://pypi.org/project/engineeringkit/) [![Python versions](https://badgen.net/pypi/python/engineeringkit)](https://pypi.org/project/engineeringkit/) [![License](https://img.shields.io/github/license/cpoli/engineeringkit)](https://github.com/cpoli/engineeringkit/blob/main/LICENSE) [![CI](https://github.com/cpoli/engineeringkit/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/cpoli/engineeringkit/actions/workflows/ci.yml) [![Docs](https://badgen.net/badge/docs/cpoli.github.io%2Fengineeringkit/blue)](https://cpoli.github.io/engineeringkit/) [![Ruff](https://badgen.net/badge/code%20style/ruff/purple)](https://github.com/astral-sh/ruff) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23452015.svg)](https://doi.org/10.5281/zenodo.23452015)

**Engineering concepts made computable → inspectable → visualizable → understandable.**

EngineeringKit is an educational scientific-Python toolkit following the design philosophy of MathematicsKit: explicit intermediate results, NumPy/SciPy for mature numerics, and hand-written formulas/algorithms where seeing the mechanics is the lesson. SI units are used by public APIs.

![Panel-method streamlines over a NACA 2412, two chaotic double pendulums released 0.001 rad apart, and the Van der Pol limit cycle computed with EngineeringKit](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_hero.png)

## Install
```bash
pip install engineeringkit
```

## Subpackages

Domain subpackages, each with runnable examples linked below:

- [`engineeringkit.aerodynamics`](https://cpoli.github.io/engineeringkit/api/gallery/aerodynamics/) -- Airfoil and wing baselines, drag polars and compressible flow.

  ![Panel-method pressure distribution on a NACA 2412, lifting-line spanwise loading, and the oblique-shock theta-beta-Mach chart](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_aerodynamics.png)

- [`engineeringkit.communications`](https://cpoli.github.io/engineeringkit/api/gallery/communications/) -- Modulation and pulse shaping, noise, capacity, detection, PCM, source and channel coding, fading and link budgets.

  ![A noisy 16-QAM constellation, a raised-cosine eye diagram, and bit error rate curves](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_communications.png)

- [`engineeringkit.control_systems`](https://cpoli.github.io/engineeringkit/api/gallery/control_systems/) -- Classical stability, compensator design, PID with anti-windup, nonlinear loops, state feedback, LQR and estimation.

  ![PID step responses, a root locus, and a Nyquist plot](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_control_systems.png)

- [`engineeringkit.electrical_engineering`](https://cpoli.github.io/engineeringkit/api/gallery/electrical_engineering/) -- DC and AC circuits, transients, electronics, machines and power systems.

  ![Series RLC step responses, a Smith chart, and induction motor torque-speed curves](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_electrical_engineering.png)

- [`engineeringkit.fluid_engineering`](https://cpoli.github.io/engineeringkit/api/gallery/fluid_engineering/) -- Hydrostatics, flow measurement, momentum, pipe friction, drag, similarity, open channels, pumps and networks.

  ![The Moody diagram, pump operating points, and open-channel specific energy curves](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_fluid_engineering.png)

- [`engineeringkit.heat_transfer`](https://cpoli.github.io/engineeringkit/api/gallery/heat_transfer/) -- Conduction, convection, radiation, fins and heat exchangers.

  ![Steady 2-D conduction in a plate, Planck blackbody spectra, and heat exchanger effectiveness-NTU curves](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_heat_transfer.png)

- [`engineeringkit.materials`](https://cpoli.github.io/engineeringkit/api/gallery/materials/) -- Crystals and diffusion, elastic constants, hardness, yield criteria, composites and laminates, fatigue and notches, fracture, creep and dislocation plasticity.

  ![Ramberg-Osgood stress-strain curves, Tresca and von Mises yield surfaces, and Paris-law fatigue crack growth](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_materials.png)

- [`engineeringkit.mechanical_engineering`](https://cpoli.github.io/engineeringkit/api/gallery/mechanical_engineering/) -- Rigid-body motion, Lagrange's equations, machine elements, mechanisms, cams and multi-degree-of-freedom vibration.

  ![Vibration transmissibility, normal modes of a five-mass chain, and a four-bar coupler curve](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_mechanical_engineering.png)

- [`engineeringkit.reliability`](https://cpoli.github.io/engineeringkit/api/gallery/reliability/) -- Lifetimes and life data, accelerated testing, networks, structural reliability, FMEA, fault trees and maintenance.

  ![The bathtub hazard curve, a Weibull probability plot, and a Kaplan-Meier survival estimate](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_reliability.png)

- [`engineeringkit.signal_processing`](https://cpoli.github.io/engineeringkit/api/gallery/signal_processing/) -- Sampling and rate conversion, the FFT and DCT, spectra and AR models, IIR and FIR design, adaptive filters.

  ![The Gibbs phenomenon, Butterworth, Chebyshev and elliptic low-pass filters, and an STFT spectrogram of a chirp](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_signal_processing.png)

- [`engineeringkit.structural_engineering`](https://cpoli.github.io/engineeringkit/api/gallery/structural_engineering/) -- Members, sections, beams on any supports, columns, stress and strain at a point, pressure vessels, energy methods, trusses, frames, cables and arches.

  ![Pratt truss member forces, a two-span continuous beam's moment and deflection, and portal frame sway](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_structural_engineering.png)

- [`engineeringkit.thermodynamics`](https://cpoli.github.io/engineeringkit/api/gallery/thermodynamics/) -- Ideal and real gases, mixtures, the two laws, steady-flow devices, engine and refrigeration cycles, exergy, steam properties, psychrometrics and combustion.

  ![An ideal Rankine cycle on the steam dome, a dual-cycle p-v diagram, and air-standard cycle efficiencies](https://raw.githubusercontent.com/cpoli/engineeringkit/main/docs/source/_static/images/readme_thermodynamics.png)

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
make -C docs html
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

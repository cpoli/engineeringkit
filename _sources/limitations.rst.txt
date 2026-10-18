Scope and remaining limitations
===============================

engineeringkit is a teaching toolkit. Each domain covers the canonical
topics of an undergraduate course, and each history entry has a gallery
example that ends in a check of the result it illustrates. This page lists
what is still open. The plan for each item is in ``ROADMAP.md``.

Out of scope
------------

The only dependencies are NumPy, SciPy and Matplotlib, so the following
are not planned: CAD and meshing, general finite element or CFD solvers,
wrappers around CoolProp, python-control or SPICE, and unit-carrying types.
Inputs and outputs are raw SI numbers. The docstrings give the unit of
every parameter and return value, but nothing stops a value in the wrong
unit from being passed.

Open numerical gaps
-------------------

- **Steam above 1073.15 K.** ``steam_properties`` implements IAPWS-IF97
  regions 1 to 4. Region 5 (1073.15 K to 2273.15 K, for gas turbines) is
  not implemented and raises ``ValueError``. ``steam_properties`` takes
  scalars only, and close to the critical point, where :math:`c_p`
  diverges, its region 3 values carry the formulation's larger tolerance.
- **Real fluids other than water.** Van der Waals's equation, the
  Joule-Thomson coefficient and the Clausius-Clapeyron and Antoine vapour
  pressures are available, but there are no property tables for
  refrigerants or other working fluids. Refrigeration stops at Carnot and
  energy-balance coefficients of performance; there is no vapour-
  compression cycle.
- **SciPy errors.** Some thin wrappers around SciPy leave shape and range
  checks to SciPy, so a bad input can raise a SciPy error rather than one
  from ``_validation.py`` (``filter_signal`` with a two-dimensional ``b``,
  for example).

Documentation gaps
------------------

- Two history entries have no primary source: Jacobi's maximum power
  theorem (``electrical_engineering``) cites a textbook, and Thoma's 1924
  cavitation paper (``fluid_engineering``) is cited by title and year only.

Using the models
----------------

The models are ideal, and their assumptions are stated, not enforced:
laminar development, Biot-number limits, small-angle attached flow,
far-field antennas, independent failures, linear elasticity and so on.
Each domain guide lists them. A result can be computed outside a model's
range of validity without any warning. The probabilistic models in
``reliability`` are not forecasts for a real system unless they are fitted
to data from it.

Release status
--------------

Nothing has been released to PyPI yet, so the Colab and JupyterLite links
of the tutorials and quickstart notebook do not work until the first
release. Releases are published to PyPI only from a ``vX.Y.Z`` tag
(``.github/workflows/release.yml``). The CI configuration shows which
checks run, not that the multi-version matrix has passed on a given
commit.

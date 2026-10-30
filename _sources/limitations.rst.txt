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

- **Steam near the critical point.** ``steam_properties`` implements all
  five IAPWS-IF97 regions. Close to the critical point, where :math:`c_p`
  diverges, its region 3 values carry the formulation's larger tolerance.
- **Real fluids other than water** come from the Peng-Robinson equation
  with a constant ideal-gas heat capacity (``refrigerant_properties``, with
  R-134a built in). Vapour pressures are within about 2% and latent heats
  within about 3% of the reference tables. Liquid volumes are 5-15% off, and
  saturation stops at 0.45 and 0.9999 of the critical temperature. There are
  no reference-quality refrigerant equations and no fluid mixtures beyond
  ideal gases.
- **Combustion** is complete combustion with no dissociation, so flame
  temperatures above about 2000 K are upper bounds, and rich mixtures (CO and
  H2 in the products) are not covered.
- **SciPy errors.** Some thin wrappers around SciPy leave shape and range
  checks to SciPy, so a bad input can raise a SciPy error rather than one
  from ``_validation.py`` (``filter_signal`` with a two-dimensional ``b``,
  for example).

Documentation gaps
------------------

- Two history entries have no primary source: Jacobi's maximum power
  theorem (``electrical_engineering``) cites a textbook, and Thoma's 1924
  cavitation paper (``fluid_engineering``) is cited by title and year only.
- In ``electrical_engineering``, Evans and Wagner's 1926 stability paper
  (swing curves and the equal-area criterion) is cited by title, journal
  and year only and was not checked, and Steinmetz's 1897 induction motor
  paper was checked against bibliographies only.

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

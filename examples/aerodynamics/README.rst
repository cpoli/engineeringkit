Examples
========

.. include:: /_generated/nav/aerodynamics.rst

This gallery walks through every public feature of ``engineeringkit.aerodynamics``,
one breakthrough per script: dynamic pressure and force coefficients,
circulation lift, airfoil sections, finite wings, cruise performance,
boundary layers, compressible flow and the standard atmosphere.

:doc:`Units, assumptions and numerical validation </guides/aerodynamics>`
gives the physical validity limits. Each script is self-contained and can be
run directly with ``python examples/aerodynamics/<section>/<script>.py``.

Sections
--------

- **forces** -- Bernoulli's dynamic pressure, Cayley's lift/drag split,
  Lilienthal's coefficients and the Wrights' correction of Smeaton's coefficient.
- **circulation** -- the Kutta-Joukowski theorem on a spinning cylinder.
- **airfoils** -- NACA four-digit sections, thin-airfoil lift and pitching
  moment, and the Hess-Smith panel method.
- **wings** -- Prandtl's lifting line, Glauert's Fourier solution, induced
  drag and the Oswald drag polar.
- **performance** -- the jet and propeller Breguet range equations.
- **boundary_layer** -- Reynolds similarity, the Blasius layer and turbulent skin friction.
- **compressible** -- sound speed, the Mach cone, stagnation states,
  Prandtl-Glauert, normal and oblique shocks, the Prandtl-Meyer expansion,
  the de Laval nozzle, and Fanno and Rayleigh duct flow.
- **atmosphere** -- the International Standard Atmosphere.

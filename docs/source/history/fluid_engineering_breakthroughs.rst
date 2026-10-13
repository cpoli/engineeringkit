Breakthroughs in Fluid Engineering
==================================

.. include:: /_generated/nav/fluid_engineering.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

Bernoulli: pressure, speed and height exchange energy
-----------------------------------------------------

An ideal streamline energy balance explains pressure changes in nozzles and elevation changes. Continuity supplies the corresponding velocity relation.

.. math::

    p+\tfrac12\rho v^2+\rho gz=\mathrm{constant}

Use :py:obj:`engineeringkit.fluid_engineering.bernoulli_pressure` to explore the model.
Steady incompressible inviscid flow along a streamline is assumed; pumps and losses need extra terms.

Poiseuille: viscosity changes the answer
----------------------------------------

Viscous shear produces a parabolic velocity profile in a fully developed circular pipe. The fourth power of radius explains the sensitivity of flow to bore size.

.. math::

    Q=\pi\Delta p R^4/(8\mu L)

Use :py:obj:`engineeringkit.fluid_engineering.poiseuille_flow` to explore the model.
Newtonian laminar flow with no slip is assumed; entrance effects and turbulence are excluded.

Reynolds and pipe correlations: similarity meets experiment
-----------------------------------------------------------

Reynolds number compares inertial and viscous effects. Darcy loss and roughness correlations extend pipe calculations beyond the laminar solution.

.. math::

    Re=\rho vD/\mu,\qquad h_f=f(L/D)v^2/(2g)

Use :py:obj:`engineeringkit.fluid_engineering.colebrook_friction_factor` to explore the model.
Use Darcy rather than Fanning friction factors. Transition flow is not reliably captured by a turbulent correlation.

Explore and verify
------------------

:doc:`/guides/fluid_engineering` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/fluid_engineering/plot_fluid_engineering_models.py ../../examples/fluid_engineering/plot_fluid_engineering_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/lecture-notes/>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

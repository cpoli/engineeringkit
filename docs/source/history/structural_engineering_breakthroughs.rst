Breakthroughs in Structural Engineering
=======================================

.. include:: /_generated/nav/structural_engineering.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

Hooke: material response becomes measurable
-------------------------------------------

The proportional relation between stress and strain connects a material test to a structural displacement. It supplies the constitutive law that equilibrium alone cannot provide.

.. math::

    \sigma=E\varepsilon,\qquad \delta=FL/(EA)

Use :py:obj:`engineeringkit.structural_engineering.axial_deformation` to explore the model.
Uniform sections, small strains and linear elasticity are assumed. Yielding requires a different constitutive model.

Euler: stability differs from strength
--------------------------------------

A slender column can lose stability before its material reaches a strength limit. Effective length captures the ideal support conditions and reveals the strong length dependence.

.. math::

    P_{cr}=\pi^2 EI/(KL)^2

Use :py:obj:`engineeringkit.structural_engineering.euler_buckling_load` to explore the model.
The formula describes ideal elastic buckling; imperfections and inelastic behavior require additional analysis.

Energy methods: loads and displacements connect
-----------------------------------------------

Castigliano and reciprocal-work ideas organize displacement calculations around stored energy. They also motivate stiffness methods, where compatibility and equilibrium are assembled across members.

.. math::

    U=F^2 L/(2EA),\qquad \delta=\partial U/\partial F

Use :py:obj:`engineeringkit.structural_engineering.strain_energy_axial` to explore the model.
Linear elastic reciprocal behavior is assumed. The derivative must actually be evaluated before using castigliano_deflection.

Explore and verify
------------------

:doc:`/guides/structural_engineering` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/structural_engineering/plot_structural_engineering_models.py ../../examples/structural_engineering/plot_structural_engineering_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://ocw.mit.edu/courses/16-21-techniques-for-structural-analysis-and-design-spring-2005/pages/lecture-notes/>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

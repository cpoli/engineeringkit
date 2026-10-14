Breakthroughs in Materials
==========================

.. include:: /_generated/nav/materials.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

Elasticity and yielding: deformation has regimes
------------------------------------------------

Hooke elasticity describes recoverable small strain. Multiaxial yield criteria replace a single normal-stress threshold with an equivalent stress measure.

.. math::

    \sigma_{VM}=\sqrt{\sigma_x^2-\sigma_x\sigma_y+\sigma_y^2+3\tau_{xy}^2}

Use :py:obj:`engineeringkit.materials.von_mises_plane_stress` to explore the model.
The plane-stress criterion is intended for isotropic ductile yielding and does not model subsequent plastic deformation.

Griffith and stress intensity: cracks alter strength
----------------------------------------------------

An energy balance explains why a crack can grow even when nominal stress is modest. Stress intensity combines applied tension, crack size and a geometry correction.

.. math::

    K_I=Y\sigma\sqrt{\pi a}

Use :py:obj:`engineeringkit.materials.stress_intensity` to explore the model.
Linear-elastic mode-I fracture requires suitable constraint and small-scale yielding; the caller supplies the geometry factor.

Fatigue: load history becomes part of the model
-----------------------------------------------

Basquin describes an empirical stress-life relation, Paris describes stable crack growth, and Miner accumulates cycle fractions. Each answers a different fatigue question.

.. math::

    da/dN=C(\Delta K)^m,\qquad D=\sum_i n_i/N_i

Use :py:obj:`engineeringkit.materials.paris_crack_growth` to explore the model.
Calibration units and validity ranges matter. Miner damage omits sequence effects and Paris growth omits threshold and rapid-fracture regimes.

Explore and verify
------------------

:doc:`/guides/materials` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/materials/plot_materials_models.py ../../examples/materials/plot_materials_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://ocw.mit.edu/courses/16-21-techniques-for-structural-analysis-and-design-spring-2005/pages/lecture-notes/>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

Breakthroughs in Heat Transfer
==============================

.. include:: /_generated/nav/heat_transfer.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

Fourier: heat follows a temperature gradient
--------------------------------------------

A constitutive law for conduction makes thermal resistance networks possible. Geometry changes the resistance expression even when material conductivity is constant.

.. math::

    \dot Q=kA\Delta T/L

Use :py:obj:`engineeringkit.heat_transfer.plane_wall_resistance` to explore the model.
One-dimensional steady conduction and constant conductivity are assumed.

Transient models: internal gradients matter
-------------------------------------------

A lumped energy balance predicts exponential cooling when the body is nearly isothermal. The Biot number indicates whether internal conduction can support that approximation.

.. math::

    (T-T_\infty)/(T_i-T_\infty)=e^{-hAt/(\rho c_p V)}

Use :py:obj:`engineeringkit.heat_transfer.lumped_temperature` to explore the model.
Check the Biot number and characteristic length; a small value supports, but does not replace, modeling judgment.

Radiation and exchangers: different mechanisms need different models
--------------------------------------------------------------------

Stefan-Boltzmann radiation depends on absolute temperature to the fourth power. Heat-exchanger effectiveness relates finite stream capacities to the maximum possible heat transfer.

.. math::

    \dot Q=\epsilon\sigma A(T_s^4-T_{sur}^4)

Use :py:obj:`engineeringkit.heat_transfer.radiation_heat_rate` to explore the model.
The radiation helper assumes a surface exchanging with a large surrounding enclosure; exchanger formulas assume their specified flow arrangement.

Explore and verify
------------------

:doc:`/guides/heat_transfer` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/heat_transfer/plot_heat_transfer_models.py ../../examples/heat_transfer/plot_heat_transfer_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/pages/lecture-notes/>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

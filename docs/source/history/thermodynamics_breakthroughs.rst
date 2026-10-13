Breakthroughs in Thermodynamics
===============================

.. include:: /_generated/nav/thermodynamics.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

The first law: heat and work share an energy balance
----------------------------------------------------

Energy conservation links heat supplied, internal-energy change and work delivered. A consistent sign convention is essential when changing between closed and flowing systems.

.. math::

    Q=\Delta U+W_{by}

Use :py:obj:`engineeringkit.thermodynamics.first_law_heat` to explore the model.
The helper uses work done by the system as positive; kinetic and potential energy need separate terms when relevant.

Carnot and entropy: conservation is not enough
----------------------------------------------

The second law restricts the conversion of heat to work. A reversible engine between two reservoirs supplies an upper efficiency bound expressed using absolute temperatures.

.. math::

    \eta_C=1-T_c/T_h

Use :py:obj:`engineeringkit.thermodynamics.carnot_efficiency` to explore the model.
Both temperatures are kelvin; finite-rate heat transfer and irreversibility lower real efficiency.

Cycles and exergy: useful work depends on the surroundings
----------------------------------------------------------

Ideal Otto, Diesel and Brayton cycles relate efficiency to compression or pressure ratios. Exergy measures the work potential of energy relative to an environment.

.. math::

    B_Q=Q(1-T_0/T_s),\qquad B_{destroyed}=T_0S_{gen}

Use :py:obj:`engineeringkit.thermodynamics.exergy_heat` to explore the model.
Ideal-gas cycle formulas use constant heat-capacity ratios. Exergy accounting requires a consistent environmental reference.

Explore and verify
------------------

:doc:`/guides/thermodynamics` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/thermodynamics/plot_thermodynamics_models.py ../../examples/thermodynamics/plot_thermodynamics_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/pages/lecture-notes/>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

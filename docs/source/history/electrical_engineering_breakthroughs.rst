Breakthroughs in Electrical Engineering
=======================================

.. include:: /_generated/nav/electrical_engineering.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

Ohm and Kirchhoff: networks become solvable
-------------------------------------------

Component laws combined with current and voltage conservation make resistor networks quantitative. Equivalent sources simplify a network as seen from two terminals.

.. math::

    V=IR

Use :py:obj:`engineeringkit.electrical_engineering.ohms_law` to explore the model.
A linear resistance and the stated passive-network assumptions are required.

Energy storage: circuits acquire memory
---------------------------------------

A capacitor stores charge, so its voltage cannot be determined by an instantaneous resistance law alone. The RC time constant sets the transient scale.

.. math::

    V_C(t)=V_s+(V_0-V_s)e^{-t/(RC)}

Use :py:obj:`engineeringkit.electrical_engineering.rc_charging_voltage` to explore the model.
Ideal constant components and a constant applied source are assumed.

Phasors: differential equations become impedance algebra
--------------------------------------------------------

Complex impedance encodes both amplitude and phase in sinusoidal steady state. Inductive and capacitive reactances cancel at series resonance.

.. math::

    Z=R+j(\omega L-1/(\omega C))

Use :py:obj:`engineeringkit.electrical_engineering.series_rlc_impedance` to explore the model.
Phasors describe steady sinusoidal operation, not startup transients or nonlinear switching.

Explore and verify
------------------

:doc:`/guides/electrical_engineering` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/electrical_engineering/plot_electrical_engineering_models.py ../../examples/electrical_engineering/plot_electrical_engineering_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://ocw.mit.edu/courses/8-02-physics-ii-electricity-and-magnetism-spring-2007/resources/cha12ac_circuits/>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

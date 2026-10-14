Breakthroughs in Control Systems
================================

.. include:: /_generated/nav/control_systems.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

Feedback: error drives correction
---------------------------------

Negative feedback uses the measured output to adjust the input. Integral action can remove steady error, while saturation can accumulate controller state unless anti-windup is used.

.. math::

    u=K_p e+K_i\int e\,dt+K_d\dot e

Use :py:obj:`engineeringkit.control_systems.PID` to explore the model.
Controller gains alone do not guarantee stability; plant dynamics and actuator limits matter.

Classical analysis: stability is visible in poles and frequency
---------------------------------------------------------------

Routh tests, root loci and Bode plots expose different aspects of a linear system. Frequency-domain margins describe proximity to a feedback instability under their stated crossover conventions.

.. math::

    1+K G(s)=0

Use :py:obj:`engineeringkit.control_systems.root_locus` to explore the model.
Sampled margin estimates depend on frequency range and resolution; inspect the guide before interpreting missing crossovers.

State space: internal dynamics enter the design
-----------------------------------------------

Controllability determines whether inputs can influence all states. LQR balances state error against control effort through a Riccati equation.

.. math::

    \dot x=Ax+Bu,\qquad J=\int_0^\infty(x^TQx+u^TRu)\,dt

Use :py:obj:`engineeringkit.control_systems.lqr` to explore the model.
Continuous-time linear dynamics and suitable stabilizability and weighting assumptions are required.

Explore and verify
------------------

:doc:`/guides/control_systems` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/control_systems/plot_control_systems_models.py ../../examples/control_systems/plot_control_systems_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/pages/lecture-notes/>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

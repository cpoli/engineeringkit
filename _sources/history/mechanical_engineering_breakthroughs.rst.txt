Breakthroughs in Mechanical Engineering
=======================================

.. include:: /_generated/nav/mechanical_engineering.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

Newton: forces predict motion
-----------------------------

An equation of motion turns a free-body diagram into an initial-value problem. Mass, stiffness and damping distinguish inertia, stored energy and dissipation.

.. math::

    m\ddot x+c\dot x+kx=F(t)

Use :py:obj:`engineeringkit.mechanical_engineering.simulate_sdof` to explore the model.
The single-degree-of-freedom model omits spatial modes and nonlinear contact.

Vibration theory: resonance becomes predictable
-----------------------------------------------

Natural frequency and damping ratio provide a common description of mechanical oscillators. Frequency ratios explain why an isolator can amplify motion near resonance.

.. math::

    \omega_n=\sqrt{k/m},\qquad \zeta=c/(2\sqrt{km})

Use :py:obj:`engineeringkit.mechanical_engineering.transmissibility` to explore the model.
Steady harmonic response assumes a linear oscillator; the undamped resonance is singular.

Machine elements: motion and power are coupled
----------------------------------------------

Gear geometry links speed reduction to torque increase. Efficiency accounts for power losses without changing the tooth-count kinematics.

.. math::

    \omega_o=\omega_i N_i/N_o

Use :py:obj:`engineeringkit.mechanical_engineering.gear_output_torque` to explore the model.
These calculations omit backlash, tooth contact stresses and transient drivetrain loads.

Explore and verify
------------------

:doc:`/guides/mechanical_engineering` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/mechanical_engineering/plot_mechanical_engineering_models.py ../../examples/mechanical_engineering/plot_mechanical_engineering_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/pages/mechanical-vibration/>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

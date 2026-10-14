Fluid engineering exercises
===========================

.. include:: /_generated/nav/fluid_engineering.rst

Three problems built on the :doc:`fluid engineering gallery </api/gallery/fluid_engineering/index>`.
Try each one before opening its solution. The solutions run when the
documentation is built, and the test suite checks every number they print.

1. Colebrook or Swamee-Jain?
----------------------------

Water at 20 °C (:math:`\rho = 998` kg/m³, :math:`\mu = 1.0\times10^{-3}`
Pa s) flows through a 100 mm commercial-steel pipe (roughness 0.045 mm).
Find the head loss per 100 m at 2 m/s with the implicit Colebrook
equation. Then find the largest difference between Colebrook and the
explicit Swamee-Jain formula over velocities from 0.1 to 5 m/s.

*Builds on:* :doc:`/api/gallery/fluid_engineering/index` (the Moody diagram).

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.fluid_engineering import (
          colebrook_friction_factor,
          darcy_head_loss,
          reynolds_number,
          swamee_jain_friction_factor,
      )

      rho, mu, D, eps = 998.0, 1.0e-3, 0.1, 0.045e-3
      Re = reynolds_number(rho, 2.0, D, mu)
      f = colebrook_friction_factor(Re, eps / D)
      print(f"Re = {Re:.3g}, f = {f:.5f}, head loss {darcy_head_loss(f, 100.0, D, 2.0):.2f} m per 100 m")
      # Re = 2e+05, f = 0.01856, head loss 3.79 m per 100 m

      v = np.linspace(0.1, 5, 50)
      Re_v = reynolds_number(rho, v, D, mu)
      cb = np.array([colebrook_friction_factor(r, eps / D) for r in Re_v])
      sj = np.array([swamee_jain_friction_factor(r, eps / D) for r in Re_v])
      print(f"largest Swamee-Jain error {100 * np.max(abs(sj / cb - 1)):.2f}%")
      # largest Swamee-Jain error 0.65%
      assert np.max(abs(sj / cb - 1)) < 0.02

   Swamee and Jain's explicit formula stays within about 1% of Colebrook
   here, well inside the scatter of the roughness value itself.

2. A venturi meter
------------------

A venturi meter with a 150 mm inlet and a 75 mm throat carries water
(:math:`\rho = 998` kg/m³) and shows a pressure difference of 20 kPa
between inlet and throat. Find the flow rate and both velocities, and check
that they satisfy continuity and Bernoulli's equation.

*Builds on:* :doc:`/api/gallery/fluid_engineering/index`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.fluid_engineering import bernoulli_pressure, continuity_velocity, venturi_flow

      rho, dp = 998.0, 20e3
      A1, A2 = np.pi * 0.150**2 / 4, np.pi * 0.075**2 / 4
      Q = venturi_flow(dp, rho, A1, A2)
      v1 = Q / A1
      v2 = continuity_velocity(A1, v1, A2)
      print(f"Q = {1e3 * Q:.2f} L/s, v1 = {v1:.3f} m/s, v2 = {v2:.3f} m/s")
      # Q = 28.89 L/s, v1 = 1.635 m/s, v2 = 6.539 m/s
      assert np.isclose(v2, 4 * v1)
      assert np.isclose(bernoulli_pressure(dp, rho, v1, v2), 0.0, atol=1e-6)

   Halving the diameter quadruples the velocity, and the pressure drop is
   the increase in dynamic pressure. A real meter reads a few percent less
   flow; its discharge coefficient corrects for the friction losses.

3. Closing a valve too fast
---------------------------

Water flows at 2 m/s through a 600 m steel pipeline from a reservoir with
50 m of head. The pressure wave speed in the pipe is 1200 m/s. Find the
Joukowsky pressure rise if a valve at the end closes instantly, compare it
with the static pressure, and find the closure time below which the full
rise occurs.

*Builds on:* :doc:`/api/gallery/fluid_engineering/index`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.fluid_engineering import hydrostatic_pressure, water_hammer_pressure

      rho, a, v, length, head = 998.0, 1200.0, 2.0, 600.0, 50.0
      surge = water_hammer_pressure(rho, a, v)
      static = hydrostatic_pressure(rho, head)
      print(f"surge {surge / 1e6:.2f} MPa, static {static / 1e6:.3f} MPa: {surge / static:.1f} times")
      print(f"closures faster than 2L/a = {2 * length / a:.1f} s see the full surge")
      # surge 2.40 MPa, static 0.489 MPa: 4.9 times
      # closures faster than 2L/a = 1.0 s see the full surge
      assert np.isclose(surge / (rho * 9.80665), a * v / 9.80665)

   The surge is worth 245 m of water, five times the static head. A valve
   closed over many multiples of :math:`2L/a` lets the reflected wave from
   the reservoir cancel most of it.

Heat transfer exercises
=======================

.. include:: /_generated/nav/heat_transfer.rst

Three problems built on the :doc:`heat transfer gallery </api/gallery/heat_transfer/index>`.
Try each one before opening its solution. The solutions run when the
documentation is built, and the test suite checks every number they print.

1. Heat loss through an insulated wall
--------------------------------------

A 10 m² house wall is 200 mm of brick (:math:`k = 0.72` W/(m K)), 50 mm of
mineral wool (0.04) and 15 mm of plaster (0.22). The inside film
coefficient is 8 W/(m² K) and the outside one 25. Find the heat loss with
20 °C inside and -5 °C outside, the temperature at each interface, and the
share of the total resistance due to the insulation.

*Builds on:* :doc:`/api/gallery/heat_transfer/conduction/plot_01_fourier_law`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.heat_transfer import composite_wall, convection_resistance, plane_wall_resistance

      A, T_in, T_out = 10.0, 293.15, 268.15
      layers = [(0.015, 0.22), (0.05, 0.04), (0.2, 0.72)]              # plaster, wool, brick: inside to outside
      wall = composite_wall(T_in, T_out, layers, A=A, h_hot=8.0, h_cold=25.0)
      print(f"heat loss {wall.heat_rate:.1f} W")
      print("interfaces [C]:", np.round(wall.interface_temperatures - 273.15, 2))
      # heat loss 142.0 W
      # interfaces [C]: [20.   18.23 17.26 -0.49 -4.43 -5.  ]

      R = [convection_resistance(8.0, A)] + [plane_wall_resistance(L, k, A) for L, k in layers] + [convection_resistance(25.0, A)]
      assert np.isclose(wall.heat_rate, (T_in - T_out) / sum(R))
      print(f"insulation: {R[2] / sum(R):.0%} of the resistance")
      # insulation: 71% of the resistance

   Five centimetres of insulation carry more than two-thirds of the wall's thermal
   resistance. The plaster's inner surface is at 18.2 °C; the surface
   temperature, not just the air temperature, sets comfort and the risk of
   condensation.

2. Is a pin fin worth it?
-------------------------

An aluminium pin fin (:math:`k = 200` W/(m K)), 5 mm in diameter and 50 mm
long, stands on a surface at 80 °C in air at 25 °C with
:math:`h = 30` W/(m² K). Treat the tip as adiabatic. Find the fin
parameter :math:`mL`, its efficiency, its heat rate, and its
effectiveness: the heat it carries over what the bare base area would.

*Builds on:* :doc:`/api/gallery/heat_transfer/conduction/plot_02_fin_efficiency`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.heat_transfer import convection_heat_rate, fin_efficiency_adiabatic, fin_heat_rate_adiabatic, fin_parameter

      k, D, L, h, Tb, Tinf = 200.0, 0.005, 0.05, 30.0, 353.15, 298.15
      P, Ac = np.pi * D, np.pi * D**2 / 4
      m = fin_parameter(h, P, k, Ac)
      eta = fin_efficiency_adiabatic(m, L)
      q = fin_heat_rate_adiabatic(h, P, k, Ac, L, Tb, Tinf)
      bare = convection_heat_rate(h, Ac, Tb, Tinf)
      print(f"mL = {m * L:.3f}, efficiency {eta:.3f}, heat rate {q:.3f} W, effectiveness {q / bare:.1f}")
      # mL = 0.548, efficiency 0.911, heat rate 1.180 W, effectiveness 36.4
      assert np.isclose(q, eta * h * P * L * (Tb - Tinf))

   The fin carries more than 36 times the heat of the base it covers, at 91%
   of the ideal of an isothermal fin. Fins pay off when :math:`h` is low,
   as in air; in boiling water the same pin would barely help.

3. Quenching a steel ball
-------------------------

A steel ball 10 mm in diameter (:math:`k = 40` W/(m K),
:math:`\rho = 7800` kg/m³, :math:`c_p = 460` J/(kg K)) leaves a furnace at
900 K and is quenched in oil at 300 K with :math:`h = 500` W/(m² K).
Check that the lumped-capacitance model applies, then find how long the
ball takes to reach 400 K.

*Builds on:* :doc:`/api/gallery/heat_transfer/transient/plot_01_biot_number`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import matplotlib.pyplot as plt
      import numpy as np

      from engineeringkit.heat_transfer import biot_number, lumped_temperature

      k, rho, cp, D, h, Ti, Tinf = 40.0, 7800.0, 460.0, 0.010, 500.0, 900.0, 300.0
      A, V = np.pi * D**2, np.pi * D**3 / 6
      Bi = biot_number(h, V / A, k)
      tau = rho * cp * V / (h * A)
      t_400 = tau * np.log((Ti - Tinf) / (400.0 - Tinf))
      print(f"Bi = {Bi:.4f}, time constant {tau:.2f} s, reaches 400 K after {t_400:.1f} s")
      # Bi = 0.0208, time constant 11.96 s, reaches 400 K after 21.4 s
      assert Bi < 0.1 and np.isclose(lumped_temperature(Ti, Tinf, h, A, rho, cp, V, t_400), 400.0)

      t = np.linspace(0, 40, 201)
      fig, ax = plt.subplots(figsize=(5, 3.2))
      ax.plot(t, lumped_temperature(Ti, Tinf, h, A, rho, cp, V, t))
      ax.axhline(400.0, color="k", ls=":")
      ax.set(xlabel="time [s]", ylabel="temperature [K]")
      fig.tight_layout()

   The Biot number is far below 0.1, so the ball's centre and surface stay
   within a few kelvin of each other and one temperature describes it.

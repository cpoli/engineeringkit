Materials exercises
===================

.. include:: /_generated/nav/materials.rst

Three problems built on the :doc:`materials gallery </api/gallery/materials/index>`.
Try each one before opening its solution. The solutions run when the
documentation is built, and the test suite checks every number they print.
The tutorial :doc:`/tutorials/beam_to_fatigue_life` puts fatigue to work
on a whole design.

1. Tresca or von Mises?
-----------------------

A point on the surface of a shaft carries a bending stress of 80 MPa and a
torsional shear stress of 50 MPa. The steel yields at 250 MPa in tension.
Find the principal stresses, the Tresca and von Mises equivalent stresses,
and the safety factor against yield by each criterion.

*Builds on:* :doc:`/api/gallery/materials/yield/plot_01_tresca_criterion` and
:doc:`/api/gallery/materials/yield/plot_02_von_mises_criterion`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.materials import factor_of_safety, tresca_equivalent, von_mises_plane_stress
      from engineeringkit.structural_engineering import mohr_circle

      sx, sy, txy, Sy = 80e6, 0.0, 50e6, 250e6
      center, radius, s1, s2 = mohr_circle(sx, sy, txy)
      tresca = tresca_equivalent([s1, s2, 0.0])
      mises = von_mises_plane_stress(sx, sy, txy)
      print(f"principal stresses {s1 / 1e6:.2f} and {s2 / 1e6:.2f} MPa")
      print(f"Tresca {tresca / 1e6:.2f} MPa (n = {factor_of_safety(Sy, tresca):.3f}), "
            f"von Mises {mises / 1e6:.2f} MPa (n = {factor_of_safety(Sy, mises):.3f})")
      # principal stresses 104.03 and -24.03 MPa
      # Tresca 128.06 MPa (n = 1.952), von Mises 117.90 MPa (n = 2.120)
      assert np.isclose(tresca, np.hypot(sx, 2 * txy)) and np.isclose(mises, np.sqrt(sx**2 + 3 * txy**2))

   Tresca is always the more conservative of the two, by up to 15% in pure
   shear. Here its equivalent stress is 9% higher.

2. How long until the crack is critical?
----------------------------------------

A wide steel plate with :math:`K_{Ic} = 50` MPa√m has a 1 mm edge crack
(:math:`Y = 1.12`) and carries a stress that cycles from 0 to 200 MPa.
Find the critical crack length, then the number of cycles to grow the
crack to it, with Paris's law :math:`da/dN = C\,\Delta K^m`,
:math:`C = 10^{-11}` m/cycle per (MPa√m)³ and :math:`m = 3`. Integrate
numerically and check with the closed form.

*Builds on:* :doc:`/api/gallery/materials/fracture/plot_02_irwin_stress_intensity` and
:doc:`/api/gallery/materials/fracture/plot_03_paris_law`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import matplotlib.pyplot as plt
      import numpy as np

      from engineeringkit.materials import critical_crack_length, paris_crack_growth, stress_intensity

      Kic, ds, Y, a0 = 50e6, 200e6, 1.12, 1e-3
      C, m = 1e-11 * (1e-6) ** 3, 3.0              # converted to (Pa sqrt(m))^-3
      ac = critical_crack_length(Kic, ds, Y)
      a = np.geomspace(a0, ac, 2001)
      dN_da = 1 / paris_crack_growth(C, m, stress_intensity(Y, ds, a))
      N = np.concatenate(([0.0], np.cumsum(np.diff(a) * (dN_da[1:] + dN_da[:-1]) / 2)))
      closed = 2 / ((m - 2) * C * (Y * ds * np.sqrt(np.pi)) ** m) * (a0 ** (1 - m / 2) - ac ** (1 - m / 2))
      print(f"critical crack length {1e3 * ac:.2f} mm")
      print(f"cycles to failure: numerical {N[-1]:.0f}, closed form {closed:.0f}")
      # critical crack length 15.86 mm
      # cycles to failure: numerical 75680, closed form 75680
      assert np.isclose(N[-1], closed, rtol=1e-5)

      fig, ax = plt.subplots(figsize=(5, 3.2))
      ax.plot(N, 1e3 * a)
      ax.set(xlabel="cycles", ylabel="crack length [mm]")
      fig.tight_layout()

   The crack spends most of its life small: half the cycles take it from
   1 mm to 2.6 mm, a sixth of the critical length. That is why inspection intervals are set from the
   smallest crack the inspection can find, not from the critical size.

3. Specific stiffness of a composite
------------------------------------

A unidirectional carbon-epoxy has carbon fibres (230 GPa, 1800 kg/m³) at a
60% volume fraction in epoxy (3.5 GPa, 1200 kg/m³). Find its modulus along
and across the fibres with the rule of mixtures, its density, and its
specific stiffness :math:`E/\rho` compared with steel (200 GPa,
7850 kg/m³) and aluminium (70 GPa, 2700 kg/m³).

*Builds on:* :doc:`/api/gallery/materials/composites/plot_01_voigt_reuss_bounds`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.materials import density_rule_of_mixtures, rule_of_mixtures_E_longitudinal, rule_of_mixtures_E_transverse

      Ef, Em, rf, rm, Vf = 230e9, 3.5e9, 1800.0, 1200.0, 0.6
      E1 = rule_of_mixtures_E_longitudinal(Ef, Em, Vf)
      E2 = rule_of_mixtures_E_transverse(Ef, Em, Vf)
      rho = density_rule_of_mixtures(rf, rm, Vf)
      print(f"E along {E1 / 1e9:.1f} GPa, across {E2 / 1e9:.2f} GPa, density {rho:.0f} kg/m^3")
      for name, E, r in (("CFRP along", E1, rho), ("CFRP across", E2, rho), ("steel", 200e9, 7850.0), ("aluminium", 70e9, 2700.0)):
          print(f"{name:11s}: E/rho = {E / r / 1e6:5.1f} MJ/kg")
      # E along 139.4 GPa, across 8.55 GPa, density 1560 kg/m^3
      # CFRP along : E/rho =  89.4 MJ/kg
      # CFRP across: E/rho =   5.5 MJ/kg
      # steel      : E/rho =  25.5 MJ/kg
      # aluminium  : E/rho =  25.9 MJ/kg
      assert E2 < E1 and np.isclose(E1, Vf * Ef + (1 - Vf) * Em)

   Along the fibres the composite is 3.5 times as stiff per kilogram as
   steel or aluminium (which are nearly identical). Across them it is
   worse than either, which is why laminates stack plies in several
   directions.

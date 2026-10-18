Structural engineering exercises
================================

.. include:: /_generated/nav/structural_engineering.rst

Three problems built on the :doc:`structural engineering gallery </api/gallery/structural_engineering/index>`.
Try each one before opening its solution. The solutions run when the
documentation is built, and the test suite checks every number they print.
The tutorial :doc:`/tutorials/beam_to_fatigue_life` follows a beam on to
its fatigue life.

1. Sizing a floor beam for stiffness
------------------------------------

A simply supported timber floor beam (:math:`E = 11` GPa) spans 4.5 m,
is 75 mm wide and carries 3 kN/m. Find the smallest depth, to the next
10 mm, that keeps the midspan deflection within span/360. Check the
deflection against :math:`5wL^4/(384EI)` and report the bending stress.

*Builds on:* :doc:`/api/gallery/structural_engineering/beams/plot_02_navier_beam`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.structural_engineering import SimplySupportedBeam, bending_stress, rectangular_second_moment

      E, L, b, w = 11e9, 4.5, 0.075, 3e3
      I_needed = 5 * w * L**4 / (384 * E * L / 360)
      h = np.ceil(1e2 * (12 * I_needed / b) ** (1 / 3)) / 1e2
      I = rectangular_second_moment(b, h)
      result = SimplySupportedBeam(L, E, I).solve_uniform_load(w)
      summary = result.summary()
      print(f"depth {1e3 * h:.0f} mm: deflection {1e3 * summary['max_abs_deflection']:.2f} mm (limit {1e3 * L / 360:.2f} mm)")
      print(f"bending stress {bending_stress(summary['max_abs_moment'], h / 2, I) / 1e6:.2f} MPa")
      # depth 270 mm: deflection 11.84 mm (limit 12.50 mm)
      # bending stress 8.33 MPa
      assert np.isclose(summary["max_abs_deflection"], 5 * w * L**4 / (384 * E * I), rtol=1e-6)
      assert summary["max_abs_deflection"] <= L / 360

   Stiffness, not strength, sizes this beam: the stress is well below the
   strength of structural timber. Since the deflection over the span goes
   as :math:`(L/h)^3` for a given load, doubling the span would need a beam
   twice as deep.

2. One column, four end conditions
----------------------------------

A steel tube (:math:`E = 200` GPa, yield 355 MPa), 100 mm outside diameter
and 5 mm wall, is a 3 m column. Find its Euler buckling load for effective
length factors :math:`K` = 0.5 (fixed-fixed), 0.7 (fixed-pinned), 1
(pinned-pinned) and 2 (a flagpole), and compare each with the squash load
:math:`A\sigma_y`. Below what slenderness does Euler's formula no longer
apply?

*Builds on:* :doc:`/api/gallery/structural_engineering/stability/plot_01_euler_buckling`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.structural_engineering import euler_buckling_load, radius_of_gyration, slenderness_ratio

      E, Sy, L, Do, t = 200e9, 355e6, 3.0, 0.100, 0.005
      Di = Do - 2 * t
      A, I = np.pi * (Do**2 - Di**2) / 4, np.pi * (Do**4 - Di**4) / 64
      r = radius_of_gyration(I, A)
      squash = A * Sy
      limit = np.pi * np.sqrt(E / Sy)
      print(f"r = {1e3 * r:.2f} mm, squash load {squash / 1e3:.0f} kN, Euler applies above KL/r = {limit:.1f}")
      for K in (0.5, 0.7, 1.0, 2.0):
          P = euler_buckling_load(E, I, L, K)
          print(f"K = {K}: KL/r = {slenderness_ratio(L, K, r):5.1f}, Euler load {P / 1e3:6.0f} kN = {P / squash:.2f} x squash")
      # r = 33.63 mm, squash load 530 kN, Euler applies above KL/r = 74.6
      # K = 0.5: KL/r =  44.6, Euler load   1481 kN = 2.80 x squash
      # K = 0.7: KL/r =  62.4, Euler load    756 kN = 1.43 x squash
      # K = 1.0: KL/r =  89.2, Euler load    370 kN = 0.70 x squash
      # K = 2.0: KL/r = 178.4, Euler load     93 kN = 0.17 x squash
      assert np.isclose(euler_buckling_load(E, I, L, 2.0) * 16, euler_buckling_load(E, I, L, 0.5))

   Euler's load exceeds the squash load in the fixed-fixed and
   fixed-pinned cases. Those columns are below the limiting slenderness of
   74.6, so they yield before they buckle, and an inelastic column formula
   (Johnson's parabola, or Engesser's tangent modulus) takes over.

3. A truss by hand and by stiffness
-----------------------------------

A triangular roof truss spans 4 m with its apex 1.5 m high and carries
10 kN down at the apex. It is pinned at the left support and on rollers at
the right. Find the member forces by the method of joints, then check them
with :class:`~engineeringkit.structural_engineering.PlanarTruss`, using
steel members of 500 mm² area. Check the apex deflection by virtual work,
:math:`\delta = \sum F f L/(EA)`, where :math:`f` are the member forces
from a unit load.

*Builds on:* :doc:`/api/gallery/structural_engineering/trusses/plot_01_whipple_method_of_joints`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.structural_engineering import PlanarTruss

      nodes = np.array([[0.0, 0.0], [4.0, 0.0], [2.0, 1.5]])
      members = ((0, 1), (0, 2), (1, 2))
      E, A, P = 200e9, 500e-6, 10e3
      lengths = np.array([np.linalg.norm(nodes[j] - nodes[i]) for i, j in members])

      sin_t, cos_t = 1.5 / 2.5, 2.0 / 2.5
      rafter = -P / (2 * sin_t)                      # compression, from the apex
      tie = -rafter * cos_t                          # tension, from a support
      hand = np.array([tie, rafter, rafter])

      truss = PlanarTruss(nodes, members, E, A)
      result = truss.solve([0, 0, 0, 0, 0, -P], fixed_dofs=[0, 1, 3])
      print("member forces [kN]: by hand", np.round(hand / 1e3, 3), " by stiffness", np.round(result.member_forces / 1e3, 3))
      # member forces [kN]: by hand [ 6.667 -8.333 -8.333]  by stiffness [ 6.667 -8.333 -8.333]
      assert np.allclose(result.member_forces, hand)

      virtual_work = np.sum(hand * (hand / P) * lengths / (E * A))
      print(f"apex deflection {1e3 * -result.displacements[5]:.4f} mm, by virtual work {1e3 * virtual_work:.4f} mm")
      # apex deflection 0.5250 mm, by virtual work 0.5250 mm
      assert np.isclose(-result.displacements[5], virtual_work)

   The truss is statically determinate (3 members + 3 reactions = 2 x 3
   joints), so the stiffness method must reproduce the hand forces exactly,
   whatever the members' stiffness. The deflection, by contrast, depends on
   :math:`EA`.

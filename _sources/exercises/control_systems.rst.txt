Control systems exercises
=========================

.. include:: /_generated/nav/control_systems.rst

Three problems built on the :doc:`control systems gallery </api/gallery/control_systems/index>`.
Try each one before opening its solution. The solutions run when the
documentation is built, and the test suite checks every number they print.

1. The largest stable gain
--------------------------

A plant :math:`G(s) = 1/((s+1)(s+2)(s+3))` is under proportional control
with gain :math:`K`. Use Routh's criterion to find the range of :math:`K`
for which the closed loop is stable, and the frequency at which it
oscillates at the limit. Check with the closed-loop poles.

*Builds on:* :doc:`/api/gallery/control_systems/stability/plot_02_routh_stability_criterion`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.control_systems import routh_first_column, routh_stable

      den = np.polymul(np.polymul([1, 1], [1, 2]), [1, 3])           # s^3 + 6 s^2 + 11 s + 6
      print("characteristic polynomial without K:", den)
      # characteristic polynomial without K: [ 1  6 11  6]

      # The s^1 row of the Routh array is (6*11 - (6 + K)) / 6, which vanishes at K = 60.
      for K in (59.0, 61.0):
          print(f"K = {K:.0f}: first column {np.round(routh_first_column(den + [0, 0, 0, K]), 3)}, stable {routh_stable(den + [0, 0, 0, K])}")
      # K = 59: first column [ 1.     6.     0.167 65.   ], stable True
      # K = 61: first column [ 1.     6.    -0.167 67.   ], stable False

      poles = np.roots(den + [0, 0, 0, 60.0])
      print("poles at K = 60:", np.round(np.sort_complex(poles), 4) + 0.0)  # + 0.0 turns -0. into 0.
      # poles at K = 60: [-6.+0.j      0.-3.3166j  0.+3.3166j]
      assert np.isclose(np.abs(poles.imag).max(), np.sqrt(11))

   The loop is stable for :math:`0 < K < 60`. At :math:`K = 60` the
   auxiliary polynomial :math:`6s^2 + 66` gives a pair of poles at
   :math:`\pm j\sqrt{11}`, an oscillation at 3.317 rad/s.

2. A servomechanism to specification
------------------------------------

A position servo behaves as a standard second-order system. It must
overshoot by no more than 10% and settle (2%) within 2 s. Find the smallest
damping ratio and natural frequency that do it, then simulate the step
response and measure both.

*Builds on:* :doc:`/api/gallery/control_systems/time_response/plot_01_hazen_servomechanism`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.control_systems import second_order_metrics, step_response, transfer_function

      ln_os = np.log(0.10)
      zeta = -ln_os / np.sqrt(np.pi**2 + ln_os**2)
      wn = 4 / (zeta * 2.0)
      metrics = second_order_metrics(zeta, wn)
      print(f"zeta = {zeta:.4f}, wn = {wn:.3f} rad/s")
      print(f"predicted: overshoot {metrics['percent_overshoot']:.2f}%, 2% settling {metrics['settling_time_2pct']:.2f} s")
      # zeta = 0.5912, wn = 3.383 rad/s
      # predicted: overshoot 10.00%, 2% settling 2.00 s

      t, y = step_response(transfer_function([wn**2], [1, 2 * zeta * wn, wn**2]), T=np.linspace(0, 5, 5001))
      outside = np.flatnonzero(abs(y - 1) > 0.02)
      print(f"simulated: overshoot {100 * (y.max() - 1):.2f}%, settles at {t[outside[-1] + 1]:.2f} s")
      # simulated: overshoot 10.00%, settles at 1.75 s
      assert np.isclose(100 * (y.max() - 1), 10.0, atol=0.01)

   The overshoot formula is exact. The settling time :math:`4/(\zeta\omega_n)`
   is an envelope bound, so the simulated response settles a little
   sooner.

3. An LQR for a double integrator
---------------------------------

A mass on a frictionless track, :math:`\ddot x = u`, has state
:math:`(x, \dot x)`. Find the LQR gain for :math:`Q = I` and
:math:`R = 1`, and check it against the hand solution
:math:`K = [1, \sqrt 3]`. Use Lyapunov's test to confirm that the closed
loop is stable.

*Builds on:* :doc:`/api/gallery/control_systems/stability/plot_03_lyapunov_direct_method`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.control_systems import lqr, lyapunov_stable

      A = np.array([[0.0, 1.0], [0.0, 0.0]])
      B = np.array([[0.0], [1.0]])
      K, P, closed_poles = lqr(A, B, np.eye(2), np.eye(1))
      print("K =", np.round(K, 6), " closed-loop poles:", np.round(closed_poles, 4))
      print("Lyapunov stable:", lyapunov_stable(A - B @ K), "(open loop:", lyapunov_stable(A), ")")
      # K = [[1.       1.732051]]  closed-loop poles: [-0.866+0.5j -0.866-0.5j]
      # Lyapunov stable: True (open loop: False )
      assert np.allclose(K, [[1.0, np.sqrt(3)]])

   The closed-loop poles are :math:`(-\sqrt3 \pm j)/2`, a damping ratio of
   0.866: the LQR with equal weights gives a well-damped loop without any
   pole placement by hand.

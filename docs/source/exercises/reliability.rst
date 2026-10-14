Reliability exercises
=====================

.. include:: /_generated/nav/reliability.rst

Three problems built on the :doc:`reliability gallery </api/gallery/reliability/index>`.
Try each one before opening its solution. The solutions run when the
documentation is built, and the test suite checks every number they print.
The tutorial :doc:`/tutorials/beam_to_fatigue_life` applies the Weibull
model to a fatigue life.

1. Three Weibull shapes
-----------------------

Three components have the same Weibull scale :math:`\eta = 1000` h and
shapes :math:`\beta = 0.5`, 1 and 3. For each, find the mean time to
failure, the :math:`B_{10}` life (the age by which 10% have failed) and
whether the hazard rate is falling, constant or rising.

*Builds on:* :doc:`/api/gallery/reliability/index`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.reliability import weibull_hazard, weibull_mttf, weibull_reliability

      eta = 1000.0
      for beta in (0.5, 1.0, 3.0):
          b10 = eta * (-np.log(0.9)) ** (1 / beta)
          assert np.isclose(weibull_reliability(b10, eta, beta), 0.9)
          trend = np.sign(weibull_hazard(2 * eta, eta, beta) - weibull_hazard(eta, eta, beta))
          print(f"beta {beta}: MTTF {weibull_mttf(eta, beta):6.1f} h, B10 {b10:6.1f} h, hazard {['falling', 'constant', 'rising'][int(trend) + 1]}")
      # beta 0.5: MTTF 2000.0 h, B10   11.1 h, hazard falling
      # beta 1.0: MTTF 1000.0 h, B10  105.4 h, hazard constant
      # beta 3.0: MTTF  893.0 h, B10  472.3 h, hazard rising
      assert np.isclose(weibull_mttf(eta, 1.0), eta)

   The infant-mortality component (:math:`\beta = 0.5`) has the longest
   mean life but the shortest :math:`B_{10}`: a few fail very early and the
   survivors last. Burn-in screens such parts; replacement schedules only
   help when :math:`\beta > 1`.

2. Redundant pumps
------------------

A cooling system has three pumps, each failing at a constant rate of
:math:`10^{-4}` per hour. Compare the probability of surviving 1000 h and
the mean time to failure when the system needs all three, any two, or any
one. Check the MTTFs against :math:`1/(3\lambda)`, :math:`5/(6\lambda)`
and :math:`11/(6\lambda)`.

*Builds on:* :doc:`/api/gallery/reliability/index`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np
      from scipy.integrate import quad

      from engineeringkit.reliability import exponential_reliability, k_of_n_reliability

      lam = 1e-4
      for k, exact in ((3, 1 / 3), (2, 5 / 6), (1, 11 / 6)):
          R = lambda t, k=k: k_of_n_reliability(3, k, float(exponential_reliability(t, lam)))
          mttf = quad(R, 0, np.inf)[0]
          print(f"{k}-out-of-3: R(1000 h) = {R(1000.0):.5f}, MTTF {mttf:6.0f} h")
          assert np.isclose(mttf, exact / lam)
      # 3-out-of-3: R(1000 h) = 0.74082, MTTF   3333 h
      # 2-out-of-3: R(1000 h) = 0.97456, MTTF   8333 h
      # 1-out-of-3: R(1000 h) = 0.99914, MTTF  18333 h

   At 1000 h, two-out-of-three is far more reliable than one pump alone
   (0.9746 against 0.9048), yet its mean time to failure is *shorter* than
   a single pump's 10,000 h. Redundancy without repair helps early and
   hurts late, because once one pump fails the remaining pair is a series
   system.

3. Availability and a safety integrity level
--------------------------------------------

A process pump has an MTBF of 5000 h and an MTTR of 24 h. Find its
availability and expected downtime per year. Its high-pressure trip has a
dangerous undetected failure rate of :math:`2\times10^{-6}` per hour and
is proof-tested once a year. Find its average probability of failure on
demand and its SIL band; how often must it be tested to reach SIL 3?

*Builds on:* :doc:`/api/gallery/reliability/index`.

.. dropdown:: Solution

   .. plot::
      :context: reset

      import numpy as np

      from engineeringkit.reliability import availability, safety_integrity_pfd

      A = availability(5000.0, 24.0)
      print(f"availability {A:.5f}, downtime {8760 * (1 - A):.1f} h a year")
      # availability 0.99522, downtime 41.8 h a year

      def sil(pfd):
          return next(level for level, upper in ((4, 1e-4), (3, 1e-3), (2, 1e-2), (1, 1e-1)) if pfd < upper)

      lam_du = 2e-6
      for months in (12, 3, 1):
          pfd = safety_integrity_pfd(lam_du, 8760 * months / 12)
          print(f"proof test every {months:2d} months: PFDavg {pfd:.2e}, SIL {sil(pfd)}")
      # proof test every 12 months: PFDavg 8.76e-03, SIL 2
      # proof test every  3 months: PFDavg 2.19e-03, SIL 2
      # proof test every  1 months: PFDavg 7.30e-04, SIL 3
      assert np.isclose(safety_integrity_pfd(lam_du, 8760), lam_du * 8760 / 2)

   Testing monthly brings the PFD down to :math:`7.3\times10^{-4}`, just
   inside SIL 3. In practice a better instrument or a second, independent
   trip is usually cheaper than testing that often.

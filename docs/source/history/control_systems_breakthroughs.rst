Breakthroughs in Control Systems
================================

.. include:: /_generated/nav/control_systems.rst

.. epigraph::

   "A governor is a part of a machine by means of which the velocity of the
   machine is kept nearly uniform, notwithstanding variations in the
   driving-power or the resistance."
   -- James Clerk Maxwell, 1868

Feedback control began as brass mechanisms on steam engines and became a
mathematical discipline when engines, amplifiers and gun directors pushed it
to instability. Maxwell, Routh and Hurwitz found when a loop is stable;
Bell Labs engineers read stability from frequency responses; Kalman and
his contemporaries recast design in state space and as optimisation.
Relays and saturating actuators, which no linear model describes, were
handled in the phase plane and by describing functions. The entries trace the ideas behind :mod:`engineeringkit.control_systems`. Each
links to the implementing function and to a gallery example devoted to it.

.. contents:: Timeline
   :local:
   :depth: 1

1788 -- Watt's Centrifugal Governor
-----------------------------------

Matthew Boulton saw centrifugal pendulums regulating the gap between
millstones, and James Watt fitted one to the rotative steam engine. Two
flyballs on a spindle driven by the engine swing outward as it speeds up
and, through a linkage, close the steam valve. The throttle moves in
proportion to the speed error: proportional feedback, built in brass. A
proportional governor holds speed far better than a fixed throttle, but under
load it settles at a lower speed. That steady *droop*,
:math:`\Delta\omega = \tau_{load}/(b + K_p)`, shrinks with gain and
never vanishes.

*Implementation:* :func:`engineeringkit.control_systems.feedback`,
:func:`engineeringkit.control_systems.step_response`, and
:class:`engineeringkit.control_systems.PID`. The example applies a load torque
step to an engine at four governor gains, zero among them, and checks the
droop formula.

*References:* O. Mayr, *The Origins of Feedback Control* (MIT Press, 1970).

.. minigallery:: ../../examples/control_systems/feedback/plot_01_watt_centrifugal_governor.py

1868 -- Maxwell's *On Governors*
--------------------------------

Engineers added integral action to governors to remove the droop, and some
of them began to "hunt", oscillating with growing amplitude. James Clerk
Maxwell linearised the equations of the engine and governor and showed that
the motion is stable exactly when every root of the characteristic
polynomial has a negative real part. For a cubic
:math:`s^3 + a_1 s^2 + a_2 s + a_3` he found the condition
:math:`a_1 a_2 > a_3`, and he asked mathematicians for a general test on the
coefficients. It was the first mathematical analysis of a feedback system.

*Implementation:* :func:`engineeringkit.control_systems.poles`. The example
sweeps the integral gain of a governor model, finds the hunting limit from
the roots, and compares it with Maxwell's condition.

*References:* J. C. Maxwell, "On governors," *Proc. R. Soc. Lond.* 16 (1868)
270-283.

.. minigallery:: ../../examples/control_systems/stability/plot_01_maxwell_on_governors.py

1877 -- Routh's Stability Criterion
-----------------------------------

Edward Routh answered Maxwell's question in the essay that won the 1877
Adams Prize. A table built from the polynomial coefficients by repeated
cross-multiplication has a first column whose sign changes count the roots
in the right half-plane. No root needs to be found, so the stability of a
loop can be written as inequalities on its gains. For
:math:`s^3 + 3s^2 + 2s + K`, the closed loop of :math:`K/(s(s+1)(s+2))`, the
table gives :math:`0 < K < 6`.

*Implementation:* :func:`engineeringkit.control_systems.routh_first_column`
and :func:`engineeringkit.control_systems.routh_stable`. The example checks
the sign-change count against computed roots over a range of gains and on a
fifth-order polynomial.

*References:* E. J. Routh, *A Treatise on the Stability of a Given State of
Motion* (Macmillan, 1877).

.. minigallery:: ../../examples/control_systems/stability/plot_02_routh_stability_criterion.py

1892 -- Lyapunov's Direct Method
--------------------------------

Aleksandr Lyapunov's doctoral thesis at Kharkov defined stability for
nonlinear motion and gave a way to prove it without solving the equations.
If a function :math:`V(x)`, positive away from the equilibrium, decreases
along every trajectory, the equilibrium is stable. For linear dynamics
:math:`\dot x = Ax` the quadratic :math:`V = x^\mathsf{T}Px` works exactly
when :math:`A^\mathsf{T}P + PA = -Q` has a positive definite solution. The
method was little known in the West until the late 1950s, when it became
the basis of nonlinear and adaptive control.

*Implementation:* :func:`engineeringkit.control_systems.lyapunov_equation`
and :func:`engineeringkit.control_systems.lyapunov_stable`. The example
solves the Lyapunov equation for a damped oscillator and checks that
:math:`V` decreases along a simulated trajectory at the rate
:math:`-x^\mathsf{T}x`.

*References:* A. M. Lyapunov, *The General Problem of the Stability of
Motion* (Kharkov, 1892; English translation, *Int. J. Control* 55 (1992)
531-773).

.. minigallery:: ../../examples/control_systems/stability/plot_03_lyapunov_direct_method.py

1895 -- Hurwitz's Determinants
------------------------------

Aurel Stodola, designing water-turbine regulators in Zürich, needed a
stability condition he could apply to high-order models and asked his
colleague Adolf Hurwitz. Hurwitz, unaware of Routh's work, stated it as
determinants. Build the matrix :math:`H_{ij} = a_{2j-i}` from the
coefficients; every root lies in the open left half-plane exactly when all
its leading principal minors are positive. The two tests are equivalent,
and the combined result is now called the Routh-Hurwitz criterion.

*Implementation:*
:func:`engineeringkit.control_systems.hurwitz_determinants`. The example
tests the quartic family :math:`(s+1)^4 + K` and finds :math:`\Delta_3`
passing through zero at the boundary :math:`K = 4`.

*References:* A. Hurwitz, "Ueber die Bedingungen, unter welchen eine
Gleichung nur Wurzeln mit negativen reellen Theilen besitzt," *Math. Ann.*
46 (1895) 273-284.

.. minigallery:: ../../examples/control_systems/stability/plot_04_hurwitz_determinants.py

1922 -- Minorsky's PID Ship Steering
------------------------------------

Nicolas Minorsky studied how helmsmen steer and found that they act on the
heading error, on its rate of change and on any error that persists. His
paper on automatic steering for the U.S. Navy gave the three-term law
:math:`\delta = K_p e + K_i\int e\,dt + K_d\dot e` and analysed its
stability. Proportional and derivative action leave a constant heading
offset under a steady wind; the integral term removes it. Trials on the
battleship *New Mexico* in 1923 confirmed the design, though the crew
preferred to steer by hand.

*Implementation:* :class:`engineeringkit.control_systems.PID`. The example
steers Nomoto's ship model against a constant wind moment with PD and PID
control and checks the offset :math:`d/K_p`.

*References:* N. Minorsky, "Directional stability of automatically steered
bodies," *J. Am. Soc. Naval Eng.* 34 (1922) 280-309.

.. minigallery:: ../../examples/control_systems/feedback/plot_02_minorsky_pid_ship_steering.py

1927 -- Black's Negative-Feedback Amplifier
-------------------------------------------

Telephone repeaters built from vacuum tubes drifted and distorted, and the
errors compounded along a transcontinental line. On the Lackawanna ferry in
August 1927, Harold Black sketched an amplifier that fed a fraction
:math:`\beta` of its output back with inverted sign. With a large open-loop
gain :math:`A` the closed-loop gain :math:`A/(1 + A\beta)` is nearly
:math:`1/\beta`, fixed by a passive divider. Gain variations and
distortion fall by the factor :math:`1 + A\beta`, and the bandwidth grows
by the same factor. The patent office took nine years to grant it.

*Implementation:* :func:`engineeringkit.control_systems.feedback` and
:func:`engineeringkit.control_systems.dc_gain`. The example varies the tube gain
by a factor of three, shows the closed-loop gain barely moving, and checks
the gain-bandwidth trade.

*References:* H. S. Black, "Stabilized feedback amplifiers," *Bell Syst.
Tech. J.* 13 (1934) 1-18.

.. minigallery:: ../../examples/control_systems/feedback/plot_03_black_negative_feedback_amplifier.py

1932 -- Nyquist's Stability Criterion
-------------------------------------

Black's amplifiers sometimes "sang", oscillating on their own. Harry Nyquist
showed how to predict this from the open-loop frequency response, which can
be measured, rather than from a model. Trace :math:`L(j\omega)` over all
frequencies. The closed loop has :math:`Z = N + P` right-half-plane poles,
where :math:`N` counts the clockwise encirclements of :math:`-1` and
:math:`P` the unstable open-loop poles. The criterion covers loops that are
unstable in open loop and loops with delays, where Routh's test cannot be
used.

*Implementation:* :func:`engineeringkit.control_systems.nyquist` and
:func:`engineeringkit.control_systems.nyquist_encirclements`. The example
counts encirclements for a stable plant at two gains and for an unstable
plant stabilised by feedback, and checks each count against the closed-loop
poles.

*References:* H. Nyquist, "Regeneration theory," *Bell Syst. Tech. J.* 11
(1932) 126-147.

.. minigallery:: ../../examples/control_systems/frequency/plot_01_nyquist_stability_criterion.py

1934 -- Hazen's Theory of Servomechanisms
-----------------------------------------

Harold Hazen at MIT studied position-following systems, in which a motor
drives a load to follow a command shaft, and named them *servomechanisms*.
His analysis treated them as second-order dynamics,
:math:`\omega_n^2/(s^2 + 2\zeta\omega_n s + \omega_n^2)`, and related
the transient to two numbers. The damping ratio :math:`\zeta` alone fixes
the overshoot, :math:`e^{-\pi\zeta/\sqrt{1-\zeta^2}}`, and with the
natural frequency it sets the peak and settling times. These remain the
specifications in which servo performance is written.

*Implementation:*
:func:`engineeringkit.control_systems.second_order_metrics`. The example reads the overshoot and
peak time from simulated step and impulse responses and checks them against
the formulas.

*References:* H. L. Hazen, "Theory of servo-mechanisms," *J. Franklin Inst.*
218 (1934) 279-331.

.. minigallery:: ../../examples/control_systems/time_response/plot_01_hazen_servomechanism.py

1937 -- Andronov, Vitt and Khaikin's Phase Plane
------------------------------------------------

Relays, valves and clutches make control loops that are not linear
anywhere, so transfer functions do not describe them. Aleksandr Andronov's
school in Gorky studied such systems by eliminating time and drawing the
state of :math:`\ddot x = f(x, \dot x)` as a curve in the plane
:math:`(x, \dot x)`. Self-sustained oscillations, which they called
auto-oscillations, appear there as isolated closed orbits, Poincaré's limit
cycles. A relay servo :math:`\ddot x = -M\,\mathrm{sign}(x)` traces a
family of closed orbits made of parabolas. With rate feedback the switching
line tilts, and the motion reaches the line and slides along it to rest.
Aleksandr Vitt was arrested as the book appeared, and his name was left off
its first edition and restored in the second (1959).

*Implementation:* :func:`engineeringkit.control_systems.phase_portrait`.
The example checks that the relay orbits conserve
:math:`\dot x^2/2 + M|x|`, and that with rate feedback the motion
reaches the switching line and decays along it as :math:`e^{-t/\tau}`.

*References:* A. A. Andronov and S. E. Khaikin, *Theory of Oscillations*
(ONTI, Moscow, 1937; English translation, Princeton University Press,
1949); A. A. Andronov, A. A. Vitt and S. E. Khaikin, *Theory of
Oscillators* (2nd ed., 1959; Pergamon, 1966). The bibliographic details
were checked against library catalogues. The account of Vitt's name comes
from secondary sources, and the 1937 text itself was not consulted.

.. minigallery:: ../../examples/control_systems/nonlinear/plot_01_andronov_phase_plane.py

1940 -- Bode's Gain and Phase Margins
-------------------------------------

Hendrik Bode showed that for a minimum-phase loop the gain and phase are
not independent: the slope of the gain curve fixes the phase. Plotted
against log-frequency, both curves become nearly straight lines, and the
distance from instability becomes two numbers. The *gain margin* is the
factor by which the gain can change at a phase crossover, the *phase
margin* the extra lag the loop tolerates at a gain crossover. A loop whose
phase crosses :math:`-180^\circ` more than once is stable only between two
gains, and its margin is the crossover nearest instability.

*Implementation:* :func:`engineeringkit.control_systems.bode` and
:func:`engineeringkit.control_systems.stability_margins`. The example checks
both margins of :math:`2/(s(s+1)(s+2))` analytically, then finds the two
gain margins of a conditionally stable loop and confirms that each puts
closed-loop poles on the imaginary axis.

*References:* H. W. Bode, "Relations between attenuation and phase in
feedback amplifier design," *Bell Syst. Tech. J.* 19 (1940) 421-454; H. W.
Bode, *Network Analysis and Feedback Amplifier Design* (Van Nostrand, 1945).

.. minigallery:: ../../examples/control_systems/frequency/plot_02_bode_stability_margins.py

1942 -- Ziegler-Nichols Tuning
------------------------------

Process plants had thousands of pneumatic controllers and no models of what
they controlled. John Ziegler and Nathaniel Nichols of Taylor Instruments
gave a recipe from one experiment. Raise a proportional gain until the loop
oscillates steadily, note this ultimate gain :math:`K_u` and period
:math:`T_u`, and read P, PI or PID settings from a table. The ultimate point
is the gain margin and phase-crossover frequency of the loop. The rules aim
for a quarter-amplitude decay, which is aggressive by later standards, but
they are still the starting point of most tuning.

*Implementation:* :func:`engineeringkit.control_systems.ziegler_nichols`.
The example finds :math:`K_u = 8` and :math:`\omega_u = \sqrt3` for
:math:`1/(s+1)^3` and compares the P, PI and PID step responses.

*References:* J. G. Ziegler and N. B. Nichols, "Optimum settings for
automatic controllers," *Trans. ASME* 64 (1942) 759-768.

.. minigallery:: ../../examples/control_systems/feedback/plot_04_ziegler_nichols_tuning.py

1945 -- Bode's Sensitivity Integral
-----------------------------------

In *Network Analysis and Feedback Amplifier Design* Hendrik Bode
measured how much of a change in the plant survives in the closed loop,
the sensitivity :math:`S = 1/(1 + L)`. Its complement
:math:`T = L/(1 + L)` carries the reference and the sensor noise to the
output, and since :math:`S + T = 1` the two cannot both be small at one
frequency. Bode also proved that a stable loop rolling off at least as
:math:`1/s^2` has :math:`\int_0^\infty \ln|S(j\omega)|\,d\omega = 0`.
Sensitivity pushed down over one band rises over another, the *waterbed
effect*. Freudenberg and Looze (1985) showed that unstable open-loop poles
raise the integral to :math:`\pi\sum\mathrm{Re}\,p_i`, a cost that feedback
cannot avoid. The peak :math:`M_s = \max|S|` is the inverse distance from
:math:`-1` and bounds both stability margins.

*Implementation:* :func:`engineeringkit.control_systems.sensitivity`,
:func:`engineeringkit.control_systems.complementary_sensitivity` and
:func:`engineeringkit.control_systems.sensitivity_peak`. The example
integrates :math:`\ln|S|` for a stable plant and for two unstable ones,
checks the results 0, :math:`\pi/2` and :math:`2\pi`, and checks the phase
margin against its :math:`M_s` bound.

*References:* H. W. Bode, *Network Analysis and Feedback Amplifier Design*
(Van Nostrand, 1945); J. S. Freudenberg and D. P. Looze, "Right half plane
poles and zeros and design tradeoffs in feedback systems," *IEEE Trans.
Autom. Control* 30 (1985) 555-565.

.. minigallery:: ../../examples/control_systems/frequency/plot_03_bode_sensitivity_integral.py

1947 -- The Nichols Chart
-------------------------

Radar fire-control servos built at the MIT Radiation Laboratory had to be
designed for a closed-loop resonant peak, while the measurements were of
the open loop. Nathaniel Nichols plotted the open-loop gain in decibels
against its phase. On these axes the closed-loop gain :math:`|L/(1 + L)|`
of a unity-feedback loop is a fixed family of contours, so the closed-loop
response can be read straight off the open-loop curve. The highest contour
the locus touches is the resonant peak :math:`M_r`. A change of gain slides
the curve up or down unchanged, and a compensator reshapes it near the
critical point :math:`(-180^\circ, 0\ \mathrm{dB})`. The chart appeared in
the laboratory's *Theory of Servomechanisms*.

*Implementation:* :func:`engineeringkit.control_systems.nichols`,
:func:`engineeringkit.control_systems.nichols_closed_loop`,
:func:`engineeringkit.control_systems.visualizers.plot_nichols` and
:func:`engineeringkit.control_systems.design_lead`. The
example draws :math:`4/(s(s+1)(s+2))` before and after lead compensation,
checks that the peak read from the chart equals the closed-loop Bode peak,
and shows the lead lowering it from 14.7 dB to 5.5 dB.

*References:* H. M. James, N. B. Nichols and R. S. Phillips (eds.), *Theory
of Servomechanisms*, MIT Radiation Laboratory Series 25 (McGraw-Hill,
1947). The bibliographic details were checked against library catalogues.

.. minigallery:: ../../examples/control_systems/frequency/plot_04_nichols_chart.py

1947 -- Tustin's Bilinear Transform
-----------------------------------

Arnold Tustin, an engineer at Metropolitan-Vickers, described signals by
their sampled values, "time series", so that a linear system became a
recurrence that could be worked by hand or machine. Replacing :math:`s` by
:math:`(2/T)(z - 1)/(z + 1)` is the trapezoidal rule for integration. It
maps the left half-plane onto the inside of the unit circle, so a stable
model gives a stable recurrence, and the imaginary axis onto the circle
itself. The price is a compressed frequency axis: :math:`\omega` appears at
:math:`(2/T)\arctan(\omega T/2)`. Prewarping, which replaces :math:`2/T` by
:math:`\omega_0/\tan(\omega_0T/2)`, makes the discrete response exact at
one chosen frequency. It is still the usual way to put a continuous
controller design into a digital controller.

*Implementation:* :func:`engineeringkit.control_systems.discretize_tustin`.
The example discretises a sharp resonance at 80% of the Nyquist frequency.
It checks the pole map :math:`z = (c + p)/(c - p)` and the shift of the
peak to :math:`(2/T)\arctan(\omega_0T/2)`, and checks that prewarping puts
the peak back at :math:`\omega_0`.

*References:* A. Tustin, "A method of analysing the behaviour of linear
systems in terms of time series," *J. IEE Part IIA* 94 (1947) 130-142. The
volume and pages were checked against citing papers, but not against the
paper itself.

.. minigallery:: ../../examples/control_systems/sampled/plot_03_tustin_bilinear_transform.py

1948 -- Evans's Root Locus
--------------------------

Walter Evans at North American Aviation wanted to see how the closed-loop
poles move as one gain changes, without factoring a polynomial for every
gain. The poles are the roots of :math:`D(s) + KN(s) = 0`, and Evans gave
rules for sketching their paths from the open-loop poles at :math:`K = 0` to
the zeros or to infinity along asymptotes. Break-away points, angles of
departure and imaginary-axis crossings come from the rules. He sold a
plastic "Spirule" for drawing them.

*Implementation:* :func:`engineeringkit.control_systems.root_locus`. The
example traces the locus of :math:`K/(s(s+2)(s+4))` and checks the
break-away point :math:`-2 + 2/\sqrt3` and the axis crossing at :math:`\pm
j\sqrt8`, :math:`K = 48`.

*References:* W. R. Evans, "Graphical analysis of control systems," *Trans.
AIEE* 67 (1948) 547-551; W. R. Evans, "Control system synthesis by root
locus method," *Trans. AIEE* 69 (1950) 66-69.

.. minigallery:: ../../examples/control_systems/root_locus/plot_01_evans_root_locus.py

1950 -- Kochenburger's Describing Function
------------------------------------------

On-off "contactor" servos were cheap and strong, but they chattered or
settled into steady oscillations that linear theory could not predict.
Ralph Kochenburger replaced the nonlinearity by its describing function
:math:`N(A)`, the gain it offers to the fundamental of a sinusoid of
amplitude :math:`A`; for an ideal relay :math:`N = 4M/(\pi A)`. The loop
can then oscillate where :math:`L(j\omega) = -1/N(A)`, so Nyquist's
picture predicts the amplitude and frequency of the limit cycle and
whether it is stable. The result is approximate, because it assumes
:math:`L` filters out the harmonics. Åström and Hägglund's relay autotuner
(1984) uses the same balance to find the ultimate point for Ziegler-Nichols
tuning.

*Implementation:* :func:`engineeringkit.control_systems.describing_function`,
:func:`engineeringkit.control_systems.relay_describing_function`,
:func:`engineeringkit.control_systems.saturation_describing_function`,
:func:`engineeringkit.control_systems.dead_zone_describing_function` and
:func:`engineeringkit.control_systems.limit_cycle`. The example predicts
the limit cycles of a relay around :math:`1/(s+1)^3` and of a saturation
around :math:`10/(s(s+1)(s+2))`, then simulates both loops. The predicted
amplitudes are within 4% of the simulated ones and the periods within 2%.

*References:* R. J. Kochenburger, "A frequency response method for
analyzing and synthesizing contactor servomechanisms," *Trans. AIEE* 69
(1950) 270-284, :doi:`10.1109/T-AIEE.1950.5060149`. The volume, pages and DOI
were checked against citing works.

.. minigallery:: ../../examples/control_systems/nonlinear/plot_02_kochenburger_describing_function.py

1952 -- Ragazzini and Zadeh: Sampled-Data Systems
-------------------------------------------------

Radar tracking and the first digital computers in the loop read the plant
only at sampling instants and hold each command until the next. John
Ragazzini and Lotfi Zadeh at Columbia gave the z-transform analysis of
such loops. For a zero-order hold the continuous model is replaced
*exactly*, at the sampling instants, by
:math:`x_{k+1} = e^{AT}x_k + \int_0^T e^{A\tau}d\tau\,B\,u_k`. The
Euler step :math:`x_{k+1} = (I + AT)x_k + BTu_k` is only a first-order
approximation of it.

*Implementation:* :func:`engineeringkit.control_systems.discretize_zoh`. The
example compares the exact recursion and Euler's step with a fine continuous
simulation of a mass-spring-damper.

*References:* J. R. Ragazzini and L. A. Zadeh, "The analysis of sampled-data
systems," *Trans. AIEE* 71 (1952) 225-234.

.. minigallery:: ../../examples/control_systems/sampled/plot_01_ragazzini_zadeh_zero_order_hold.py

1957 -- Smith's Predictor for Dead Time
---------------------------------------

In chemical plants the output is often measured only after a transport
delay :math:`\theta`. The delay adds phase lag :math:`\omega\theta`
without changing the gain, so a controller tuned for the delay-free plant
oscillates once the delay is added. Otto Smith at Berkeley put a model of
the plant inside the controller. It feeds back the undelayed model output
and corrects it by the mismatch between the measurement and the delayed
model output. With a perfect model the closed loop is
:math:`CGe^{-s\theta}/(1 + CG)`: the delay-free design, delayed.

*Implementation:* :func:`engineeringkit.control_systems.pade_delay` and
:func:`engineeringkit.control_systems.smith_predictor`. The example
simulates a PI loop around :math:`e^{-2s}/(s+1)` with a true delay buffer,
shows it diverging without the predictor, and checks that with it the output
matches the delayed delay-free response.

*References:* O. J. M. Smith, "Closer control of loops with dead time,"
*Chem. Eng. Prog.* 53 (1957) 217-219.

.. minigallery:: ../../examples/control_systems/time_response/plot_02_smith_predictor.py

1960 -- Kalman's Controllability and Observability
--------------------------------------------------

A transfer function describes only the part of a system that the input
reaches and the output sees. Rudolf Kalman split the state of
:math:`\dot x = Ax + Bu,\ y = Cx` into the parts that are controllable and
observable. The state can be steered anywhere only if
:math:`[B, AB, \dots, A^{n-1}B]` has full rank, and reconstructed from the
output only if :math:`[C; CA; \dots; CA^{n-1}]` does. A mode that fails
either test vanishes from the transfer function in a pole-zero
cancellation. It still exists, and if it is unstable no output feedback can
fix it.

*Implementation:* :func:`engineeringkit.control_systems.is_controllable` and
:func:`engineeringkit.control_systems.is_observable`. The example shows two
tanks on one pump whose level difference ignores the input, and a mode
hidden by a cancellation along which the state moves while the output stays
exactly zero.

*References:* R. E. Kalman, "On the general theory of control systems,"
*Proc. 1st IFAC Congress* (Moscow, 1960) 481-492.

.. minigallery:: ../../examples/control_systems/state_space/plot_01_kalman_controllability_observability.py

1960 -- Kalman's Linear-Quadratic Regulator
-------------------------------------------

Kalman also asked for the state feedback :math:`u = -Kx` that minimises
:math:`\int_0^\infty (x^\mathsf{T}Qx + u^\mathsf{T}Ru)\,dt`. The
answer is :math:`K = R^{-1}B^\mathsf{T}P`, where :math:`P` solves an
algebraic Riccati equation, and the minimum cost from :math:`x_0` is
:math:`x_0^\mathsf{T}Px_0`. The weights replace pole locations as the
designer's knobs, and the same design works for any number of inputs. In
discrete time it is the steady state of Bellman's dynamic-programming
recursion.

*Implementation:* :func:`engineeringkit.control_systems.lqr` and
:func:`engineeringkit.control_systems.dlqr`. The example checks the
closed-form gain of the double integrator, compares the simulated cost with
:math:`x_0^\mathsf{T}Px_0`, and shows the discrete gain converging to the
continuous one as the sampling period shrinks.

*References:* R. E. Kalman, "Contributions to the theory of optimal
control," *Bol. Soc. Mat. Mexicana* 5 (1960) 102-119; R. Bellman, *Dynamic
Programming* (Princeton University Press, 1957).

.. minigallery:: ../../examples/control_systems/optimal/plot_01_kalman_lqr.py

1961 -- Jury's Stability Test for Sampled Systems
-------------------------------------------------

A sampled loop is stable when every root of its characteristic polynomial
in :math:`z` lies inside the unit circle, so Routh's test does not apply.
Eliahu Jury at Berkeley, with Blanchard, gave the discrete counterpart as a
table. Each row reduces the polynomial by one degree, and the polynomial is
stable exactly when :math:`|a_n| < |a_0|` at every step. A loop that is
stable for every gain in continuous time has a finite gain limit once it is
sampled, and the limit falls as the sampling period grows.

*Implementation:* :func:`engineeringkit.control_systems.jury_stable`. The
example checks Jury's verdict against the roots for the sampled loop of
:math:`K/(s(s+1))`, finds the critical gain by bisection and compares it
with :math:`(1 - a)/b_2` over a range of periods.

*References:* E. I. Jury and J. Blanchard, "A stability test for linear
discrete systems in table form," *Proc. IRE* 49 (1961) 1947-1948.

.. minigallery:: ../../examples/control_systems/sampled/plot_02_jury_stability_test.py

1961 -- The Kalman-Bucy Filter
------------------------------

A state estimator has to trade trust in the model against trust in the
measurement. Kalman's discrete filter (1960) and its continuous-time form
with Richard Bucy give the gain that minimises the estimation error for a
linear system driven by white process noise of spectral density :math:`Q`
and measured through white noise of density :math:`R`. The error covariance
solves the LQR Riccati equation transposed, a duality Kalman pointed out.
The filter navigated the Apollo spacecraft.

*Implementation:* :func:`engineeringkit.control_systems.kalman_gain`. The
example estimates a random walk from noisy measurements and checks the error
variance, :math:`\sqrt{QR}` at the optimal gain and :math:`(Q +
\ell^2R)/2\ell` at others, by simulation.

*References:* R. E. Kalman and R. S. Bucy, "New results in linear filtering
and prediction theory," *J. Basic Eng.* 83 (1961) 95-108.

.. minigallery:: ../../examples/control_systems/optimal/plot_02_kalman_bucy_filter.py

1964 -- Luenberger's Observer
-----------------------------

State feedback needs every state, but sensors measure only some. David
Luenberger, in his Stanford thesis, built the missing states from a copy of
the model driven by the same input and corrected by the output error,
:math:`\dot{\hat x} = A\hat x + Bu + L(y - C\hat x)`. The error then
obeys :math:`\dot e = (A - LC)e` whatever the input, and for an observable
pair the gain :math:`L` places those poles anywhere. Without noise in the
problem it is a design choice, not an optimum. The Kalman-Bucy filter is
the observer whose gain is chosen for the noise.

*Implementation:* :func:`engineeringkit.control_systems.observer_gain`. The
example estimates the velocity of a lightly damped mass from its position
and checks that the error follows :math:`e^{(A-LC)t}e_0` exactly.

*References:* D. G. Luenberger, "Observing the state of a linear system,"
*IEEE Trans. Mil. Electron.* 8 (1964) 74-80.

.. minigallery:: ../../examples/control_systems/state_space/plot_02_luenberger_observer.py

1967 -- Wonham's Pole Assignment
--------------------------------

For a single input, controllability lets state feedback place the poles of
:math:`A - BK` anywhere, with a unique gain. Murray Wonham proved the
multi-input case: the closed-loop poles can be assigned arbitrarily if and
only if :math:`(A, B)` is controllable, even when no single input reaches
every mode. The extra freedom in :math:`K` was later used to shape
eigenvectors and robustness. A precompensator then scales the reference so
the output settles on it.

*Implementation:*
:func:`engineeringkit.control_systems.state_feedback_poles` and
:func:`engineeringkit.control_systems.feedforward_gain`. The example
stabilises an inverted pendulum on a cart, checks the placed poles and the
unit steady state of the cart position, and places a complex pair with two
inputs that are each insufficient alone.

*References:* W. M. Wonham, "On pole assignment in multi-input controllable
linear systems," *IEEE Trans. Autom. Control* 12 (1967) 660-665.

.. minigallery:: ../../examples/control_systems/state_space/plot_03_wonham_pole_assignment.py

1978 -- Doyle: Guaranteed Margins for LQG? There Are None
---------------------------------------------------------

The LQR loop is robust. Its Nyquist plot stays outside the unit disc about
:math:`-1`, so it tolerates any gain in :math:`[1/2, \infty)` and at least
60° of phase change. It seemed natural that an LQR gain fed by a
Kalman-Bucy filter would keep these margins. John Doyle's two-page note,
whose abstract reads "There are none.", gave a plant for which the margins
of the LQG loop shrink towards zero as the weights grow. The result started
the robust-control work of the 1980s, from loop-transfer recovery to
:math:`H_\infty` design.

*Implementation:* :func:`engineeringkit.control_systems.lqr`,
:func:`engineeringkit.control_systems.kalman_gain`, and
:func:`engineeringkit.control_systems.stability_margins`. The example
computes the plant-input margins of Doyle's example over all crossovers for
growing weights and checks the LQR guarantees against the LQG gain range
closing on 1.

*References:* J. C. Doyle, "Guaranteed margins for LQG regulators," *IEEE
Trans. Autom. Control* 23 (1978) 756-757.

.. minigallery:: ../../examples/control_systems/optimal/plot_03_doyle_lqg_margins.py

See Also
--------

- :doc:`/api/control_systems`
- :doc:`/guides/control_systems`
- :doc:`/history/signal_processing_breakthroughs`
- :doc:`/history/communications_breakthroughs`

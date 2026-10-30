Control systems: matrix contracts and numerical stability
=========================================================

All systems are continuous-time. For xdot=A*x+B*u and y=C*x+D*u, A is (n,n),
B is (n,m), C is (p,n) and D is (p,m). Matrix units depend on the selected
state/input/output units: A entries have units x_i/(s*x_j), B entries
x_i/(s*u_j), and C entries y_i/x_j. Raw matrices do not carry unit metadata.
Explicit two-dimensional matrices are required, including SISO B and C.

LQR and Kalman estimation
-------------------------

LQR minimizes integral(x.T*Q*x+u.T*R*u)dt; Q is symmetric positive semidefinite
and R symmetric positive definite. Its return tuple is (K,P,closed_loop_poles).
P solves A.T*P+P*A-P*B*inv(R)*B.T*P+Q=0, and feedback is u=-K*x. Q/R units
must make the two quadratic terms consistent. Existence still requires the
appropriate stabilizability/detectability conditions; SciPy errors remain
visible if no stabilizing solution exists.

Kalman-Bucy uses process spectral covariance Q and measurement spectral
covariance R, not sampled-data covariance. The return is (L,P,estimator_poles).
The residual is A*P+P*A.T-P*C.T*inv(R)*C*P+Q=0. This function supplies a
steady-state gain; it does not run a discrete measurement filter.

.. code-block:: python

    from engineeringkit.control_systems import lqr
    K, P, poles = lqr([[0]], [[1]], [[4]], [[1]])
    # Scalar analytic optimum: K=2, P=2, pole=-2 [1/s].

Tests check scalar analytic optima, matrix Riccati residuals, symmetry, invalid
weights and matrix shapes, and list-to-array conversion.

Classical methods
-----------------

Polynomials are in descending powers of s [1/s]. Root locus requires a strictly
proper open-loop model, and row order follows gains; roots within a row are
not branch-tracked. Routh uses an auxiliary-polynomial derivative for zero
rows, and a positive epsilon for zero pivots. Strict stability also excludes
imaginary-axis poles using a relative root tolerance. Near multiple/marginal
roots need separate analysis. Stability margins are grid estimates read at every
crossover: gain margins at each phase crossing of -180+360k deg, phase margins
(wrapped into [-180, 180) deg) at each unity-gain crossing. The scalar
gain_margin is the one nearest instability on a log scale (it is below 1 for
the lower limit of a conditionally stable loop), phase_margin_deg the one of
smallest magnitude; the arrays list all of them. A grid that misses a crossing
misses its margin, negative loop gains are not handled, and time delays must be
replaced by pade_delay first.

Dead time and sampled loops
---------------------------

pade_delay(theta, n) returns the diagonal (n, n) Padé approximant of
exp(-s*theta) [theta in s]: all-pass, with n right-half-plane zeros, and with
accurate phase for omega*theta up to about n. smith_predictor(C, G, theta, n)
returns C/(1+C*G*(1-D)) with D the Padé factor; it assumes a perfect model and
cannot stabilise an unstable plant. jury_stable takes a polynomial in z
(descending powers) and tests abs(z) < 1 for every root without computing them;
roots on the unit circle count as unstable, with a 1e-12 relative tolerance.

discretize_tustin(G, dt, prewarp_frequency) substitutes s = c*(z-1)/(z+1)
with c = 2/dt [1/s], or c = w0/tan(w0*dt/2) to match G exactly at w0 [rad/s],
which must be below pi/dt. It needs a proper model, returns a monic
TransferFunction in z with dt set, keeps stable models stable, and shows a
continuous frequency w at (2/dt)*arctan(w/c).

PID derivative action is ideal, with no roll-off. When ki=0 no integrator is
included. Second-order metrics require 0<zeta<1; 4/(zeta*omega_n) is an
approximate 2% settling-time estimate, not exact.

PID.simulate_saturated(plant, T, u_min, u_max, anti_windup, tracking_time)
clips the actuator to [u_min, u_max] and integrates the loop with fixed RK4
steps on the grid T, so T must be fine. anti_windup=None lets the integrator
wind up; "clamping" stops it while the actuator is saturated in the direction
the error pushes; "back_calculation" adds (u - v)/Tt, with Tt defaulting to
sqrt(Ti*Td), or Ti without derivative action. The plant must be strictly
proper; derivative action acts on the measurement and needs relative degree
at least 2. The result carries the actuator signal u.

Compensators and steady-state error
-----------------------------------

lead_compensator(phi_m, w_m) returns (T*s+1)/(alpha*T*s+1) with
alpha=(1-sin phi_m)/(1+sin phi_m) and T=1/(w_m*sqrt(alpha)): unit DC gain,
peak lead phi_m (0<phi_m<90 deg) at w_m, gain 1/sqrt(alpha) there.
lag_compensator(beta, w_z) returns (s/w_z+1)/(beta*s/w_z+1), beta>1: unit DC
gain, high-frequency gain 1/beta. Both keep the error constants, which the
loop gain sets first.

design_lead(L, pm) follows Ogata: lead = pm - pm0 + extra_deg, placed where
abs(L) = sqrt(alpha). design_lag(L, pm) moves the crossover down to where the
phase is -180 + pm + extra_deg, with beta = abs(L) there and the zero
zero_ratio times lower. design_lead_lag(L, pm, wc) puts the lead's peak at
the chosen crossover wc and lets the lag remove the excess gain. Each returns
a CompensatorResult with the margins actually achieved, read from
stability_margins on the default grid; extra_deg (5 deg by default) absorbs
the error of the one-shot recipe, so check the result rather than assume it.
Designs that one stage cannot deliver (lead of 90 deg or more, a lag that
would have to raise the gain) raise ValueError.

system_type(L) counts exact zero coefficients: poles at s=0 minus zeros at
s=0. error_constants(L) returns Kp [-], Kv [1/s], Ka [1/s^2] and the
steady-state errors to a unit step, a unit-slope ramp and t^2/2, and raises
if the closed loop is unstable, since then no steady state exists.

Sensitivity and the Nichols chart
---------------------------------

sensitivity(L) = 1/(1+L) and complementary_sensitivity(L) = L/(1+L), with no
cancellation, so S+T=1. sensitivity_peak(L, w) reads Ms and Mt on a grid and
returns the guaranteed margins GM >= Ms/(Ms-1) and PM >= 2*arcsin(1/(2*Ms));
it does not check closed-loop stability. Bode's integral of ln(abs(S)) is zero for
a stable loop of relative degree two or more, and pi times the sum of the
unstable open-loop poles otherwise.

nichols(L, w) returns the open-loop gain [dB] and unwrapped phase [deg], the
closed-loop gain and phase, and the resonant peak on the grid.
nichols_closed_loop(gain_db, phase_deg) converts chart coordinates to the
closed-loop gain and phase; it is periodic in phase and grows without bound
at (-180 deg, 0 dB). plot_nichols draws the M-contours (default -20 to 6 dB).

Nonlinear loops
---------------

describing_function(f, A) takes the first Fourier harmonic of a memoryless,
vectorised f by the rectangle rule (samples points per period): spectrally accurate for
smooth f, first-order accurate across jumps. The closed forms are
relay_describing_function (4M/(pi*A), times exp(-j*arcsin(h/A)) with
hysteresis h<A), saturation_describing_function and
dead_zone_describing_function (gain slope, breakpoint limit or width, all in
input units). limit_cycle(L, N, (A_min, A_max)) solves N(A) = 1/abs(L) at every
phase crossover of L for a real N and reports amplitude, frequency, period
and Loeb's stability verdict, which assumes an open-loop stable L. It is a
first-harmonic estimate: check it by simulation when L filters harmonics
poorly.

phase_portrait(f, x0, v0, t) integrates xdot=v, vdot=f(x, v) with fixed RK4
steps on t; a discontinuous f such as a relay is first-order accurate at
each switch and chatters on a sliding line at the step size.

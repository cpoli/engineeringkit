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
roots need separate analysis. Stability margins are grid estimates using the
first detected crossover; they do not certify stability for multiple crossovers,
negative gains, time delays or a frequency grid that misses crossings.

PID derivative action is ideal, with no roll-off, saturation or anti-windup.
When ki=0 no integrator is included. Second-order metrics require 0<zeta<1;
4/(zeta*omega_n) is an approximate 2% settling-time estimate, not exact.

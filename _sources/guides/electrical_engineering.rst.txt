Electrical engineering: passive components and circuit responses
================================================================

Voltages are V, current is A, resistance is ohm, capacitance is F, inductance
is H and angular frequency is rad/s. RMS phasors are assumed for power helpers.
Negative voltage/current are signed quantities; passive resistances are
nonnegative, and denominators used to infer current/resistance must be nonzero.
Parallel-resistance helpers require strictly positive branches; represent an
ideal short circuit separately rather than passing a zero branch.

RLC current convention
-----------------------

SeriesRLC models current for an applied voltage. The transfer function is
I(s)/V(s)=s/(L*s^2+R*s+1/C), with zero initial energy. The previous numerator
1 represented charge, despite the response being named current. That error is
corrected here. The response to a 1 V step approaches zero current for R>0 and
has integrated charge C. User-supplied transient grids must be uniformly spaced
and strictly increasing, as required by SciPy step response. Lossless R=0 is allowed and oscillates without decay.

.. code-block:: python

    import numpy as np
    from engineeringkit.electrical_engineering import SeriesRLC
    result = SeriesRLC(2, 1, .25).step_response(np.linspace(0, 20, 20001))
    fig, ax = result.plot()  # current [A]

For the underdamped case the independent benchmark uses
exp(-R*t/(2L))*sin(wd*t)/(L*wd), wd^2=1/(LC)-(R/(2L))^2. Tests also integrate
current to obtain C, verify resonance produces purely resistive impedance,
and check source energy equals capacitor storage plus resistor dissipation
for RC charging. Capacitor impedance excludes zero frequency, where the
ideal capacitor is an open circuit.

Transmission-line impedance and velocity are lossless models using H/m and
F/m. Thévenin/Norton helpers transform supplied equivalent parameters; they
do not solve an arbitrary network. Power factor is bounded by [-1,1], allowing
signed real power. Balanced sinusoidal operation is assumed for three-phase
power; the result does not cover unbalanced phases or harmonic distortion.

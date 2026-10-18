Breakthroughs in Electrical Engineering
=======================================


.. include:: /_generated/nav/electrical_engineering.rst

.. epigraph::

   "The force of the current in a galvanic circuit is directly as the sum of
   all the tensions, and inversely as the entire reduced length of the
   circuit."
   -- Georg Simon Ohm, 1827

Circuit theory began as a branch of physics, with Ohm's wires and
Kirchhoff's junctions. The telegraph, the electric light and the telephone
then turned it into engineering: equivalent circuits, impedance and the
transmission line were invented by telegraph and power engineers to answer
practical questions. The entries trace the ideas behind
:mod:`engineeringkit.electrical_engineering`, from a thermocouple and a
galvanometer to a digital computer solving the power grid. Each links to the
implementing function and to a gallery example devoted to it.

.. contents:: Timeline
   :local:
   :depth: 1

1827 -- Ohm's Law
-----------------

Georg Simon Ohm, a schoolteacher in Cologne, replaced unsteady voltaic
piles with a thermocouple and measured the current through brass wires of
different lengths with a torsion galvanometer. His readings fitted
:math:`X = a/(b + x)`, an "exciting force" divided by the resistance of the
wire plus that of the rest of the circuit. *Die galvanische Kette* (1827)
framed the result by analogy with Fourier's heat flow and was received
coldly in Germany, until the Royal Society awarded Ohm its Copley Medal in
1841. In modern form, :math:`V = IR`, and resistances add in series while
conductances add in parallel.

*Implementation:* :func:`engineeringkit.electrical_engineering.ohms_law`,
:func:`engineeringkit.electrical_engineering.series_resistance` and
:func:`engineeringkit.electrical_engineering.parallel_resistance`. The
example rebuilds Ohm's experiment and recovers the source EMF and internal
resistance from the readings.

*References:* G. S. Ohm, *Die galvanische Kette, mathematisch bearbeitet*
(Berlin: T. H. Riemann, 1827).

.. minigallery:: ../../examples/electrical_engineering/dc_circuits/plot_01_ohm_law.py

c. 1840 -- Jacobi's Maximum Power Theorem
-----------------------------------------

Moritz Hermann von Jacobi built some of the first practical electric motors
in St Petersburg, one of which drove a boat on the Neva in 1839. Asking how
to get the most power from his batteries, he found that a source of EMF
:math:`V` and internal resistance :math:`R_s` delivers the most power,
:math:`V^2/4R_s`, to a load equal to its own resistance. The result is
often misread as a recipe for efficiency: at the match half the power is
lost inside the source. Power systems therefore run far from the match,
while radio receivers and amplifiers, where the signal matters more than
the energy, are matched deliberately.

*Implementation:* :func:`engineeringkit.electrical_engineering.maximum_power_transfer`.
The example sweeps the load and plots the delivered power against the
efficiency.

*References:* M. H. Jacobi's motor experiments of 1837-1840; for a modern
statement see C. K. Alexander and M. N. O. Sadiku, *Fundamentals of
Electric Circuits* (New York: McGraw-Hill), ch. 4.

.. minigallery:: ../../examples/electrical_engineering/networks/plot_01_jacobi_maximum_power.py

1843 -- Wheatstone's Bridge
---------------------------

Samuel Hunter Christie described the circuit in 1833. Charles Wheatstone
credited him and made it the standard instrument for measuring resistance
in his 1843 Bakerian lecture. Two voltage dividers are connected across one
source and a galvanometer joins their midpoints. When
:math:`R_1/R_2 = R_3/R_4` the galvanometer reads zero, whatever the source
voltage. A null needs no calibrated meter, only calibrated resistors, so
the bridge measured resistances, and later capacitances and inductances,
more accurately than any deflection instrument.

*Implementation:* :func:`engineeringkit.electrical_engineering.wheatstone_balanced`.
The example balances a bridge against an unknown resistor, using nodal
analysis for the galvanometer current.

*References:* C. Wheatstone, "An account of several new instruments and
processes for determining the constants of a voltaic circuit," *Phil.
Trans. R. Soc.* 133 (1843) 303-327.

.. minigallery:: ../../examples/electrical_engineering/dc_circuits/plot_02_wheatstone_bridge.py

1845 -- Kirchhoff's Circuit Laws and Nodal Analysis
---------------------------------------------------

Gustav Kirchhoff, then a 21-year-old student of Franz Neumann in
Königsberg, ended a paper on current flow in a plane conductor with two
laws for networks of wires: the currents meeting at a junction sum to
zero, and around any closed loop the EMFs equal the sum of the resistance
drops. In 1847 he showed how many independent equations a network gives,
an early result in graph theory. Writing the current law at every node
gives the linear system :math:`G\mathbf{v} = \mathbf{i}`. With one extra row
per voltage source (modified nodal analysis, 1975), it is how SPICE and
every circuit simulator still work.

*Implementation:* :func:`engineeringkit.electrical_engineering.nodal_analysis`,
returning a :class:`~engineeringkit.electrical_engineering.NodalResult`.
The example solves a two-source ladder, checks both laws and verifies
superposition.

*References:* G. Kirchhoff, "Ueber den Durchgang eines elektrischen Stromes
durch eine Ebene, insbesondere durch eine kreisförmige," *Ann. Phys.* 64
(1845) 497-514; C.-W. Ho, A. E. Ruehli and P. A. Brennan, "The modified
nodal approach to network analysis," *IEEE Trans. Circuits Syst.* 22
(1975) 504-509.

.. minigallery:: ../../examples/electrical_engineering/dc_circuits/plot_03_kirchhoff_nodal_analysis.py

1855 -- Kelvin's Telegraph Cable and the RC Time Constant
---------------------------------------------------------

Faraday had noticed that a submarine cable, a copper core in gutta-percha
under seawater, behaves like a long capacitor. William Thomson, later Lord
Kelvin, worked out the consequence. With inductance neglected, the cable
obeys a diffusion equation like Fourier's heat equation: a pulse arrives
smeared out, over a time proportional to :math:`rc\,\ell^2`. This *law of
squares* set the signalling speed of the first transatlantic cables. Its
building block is one resistor charging one capacitor, which approaches the
source voltage as :math:`1 - e^{-t/RC}`.

*Implementation:* :func:`engineeringkit.electrical_engineering.rc_time_constant`,
:func:`engineeringkit.electrical_engineering.rc_charging_voltage` and
:func:`engineeringkit.electrical_engineering.rl_time_constant`. The example
checks one RC section, then chains sections into a cable and recovers the
law of squares.

*References:* W. Thomson, "On the theory of the electric telegraph," *Proc.
R. Soc. London* 7 (1855) 382-399.

.. minigallery:: ../../examples/electrical_engineering/transients/plot_01_kelvin_rc_cable.py

1876-1887 -- Heaviside's Telegrapher's Equations and the Distortionless Line
----------------------------------------------------------------------------

Oliver Heaviside, a former telegraph operator who taught himself
mathematics, kept all four constants of a line: series resistance and
inductance, shunt conductance and capacitance per unit length. The
resulting telegrapher's equations carry waves that decay and travel
according to :math:`\gamma = \sqrt{(R + j\omega L)(G + j\omega C)}`. In 1887
he showed that if :math:`R/L = G/C` every frequency is attenuated equally and
travels at the same speed, so signals arrive undistorted. Telephone cables
had far too little inductance. The Post Office ignored him, but Pupin's and
Campbell's loading coils (1899-1900) added the missing inductance and made
long-distance telephony work.

*Implementation:* :func:`engineeringkit.electrical_engineering.propagation_constant`,
:func:`engineeringkit.electrical_engineering.distortionless_condition`,
:func:`engineeringkit.electrical_engineering.transmission_line_z0` and
:func:`engineeringkit.electrical_engineering.transmission_line_velocity`.
The example compares attenuation, velocity and a received pulse for an
unloaded and a distortionless line.

*References:* O. Heaviside, "On the extra current," *Phil. Mag.* 2 (1876)
135-145; O. Heaviside, "Electromagnetic induction and its propagation,"
*The Electrician* (1887), reprinted in *Electrical Papers*, vol. 2 (London:
Macmillan, 1892).

.. minigallery:: ../../examples/electrical_engineering/transmission_lines/plot_01_heaviside_telegrapher.py

1883 -- Thévenin's Equivalent Circuit
-------------------------------------

Léon Charles Thévenin, an engineer of the French Postes et Télégraphes,
proved that any network of EMFs and resistances, seen from two of its
terminals, acts as a single EMF in series with a single resistance.
Hermann von Helmholtz had stated the result in 1853, but it was Thévenin's
paper that engineers read. The EMF is the open-circuit voltage. The
resistance is the open-circuit voltage divided by the short-circuit
current, or equally the resistance seen with every source zeroed. The
theorem lets an engineer analyse what is connected to a network without
solving the network again for each load.

*Implementation:* :func:`engineeringkit.electrical_engineering.thevenin_from_source`,
with :func:`engineeringkit.electrical_engineering.nodal_analysis` for the
open- and short-circuit tests. The example reduces a two-source network and
checks the equivalent for every load.

*References:* L. Thévenin, "Extension de la loi d'Ohm aux circuits
électromoteurs complexes," *Annales Télégraphiques* 10 (1883) 222-224.

.. minigallery:: ../../examples/electrical_engineering/networks/plot_02_thevenin_equivalent.py

1886 -- Stanley's Transformer
-----------------------------

Gaulard and Gibbs's "secondary generator" and the closed-core transformer
of Zipernowsky, Bláthy and Déri (1885) showed that AC voltage could be
changed with coils on iron. William Stanley, for George Westinghouse,
connected the primaries in parallel across the line and in March 1886 lit
Great Barrington, Massachusetts, stepping 500 V down to 100 V at each
building. An ideal transformer scales voltage by the turns ratio, current
by its inverse and impedance by its square. Because line loss falls as the
square of the transmission voltage, the transformer won the "war of the
currents" for AC.

*Implementation:* :func:`engineeringkit.electrical_engineering.transformer_secondary`.
The example checks the power and impedance transformation and the
:math:`1/V^2` line loss.

*References:* W. Stanley, "Induction-coil," U.S. Patent 349,611 (1886).

.. minigallery:: ../../examples/electrical_engineering/ac_circuits/plot_01_stanley_transformer.py

1888 -- Tesla's Polyphase System
--------------------------------

Nikola Tesla's lecture to the AIEE described motors driven by currents
displaced in phase. Together they produce a rotating magnetic field that
turns an induction motor with no commutator or brushes. Westinghouse bought
the patents. Mikhail Dolivo-Dobrovolsky chose three phases for the
Lauffen-Frankfurt line of 1891. A balanced three-phase load draws constant
instantaneous power, :math:`\sqrt3 V_L I_L\cos\phi`, so generators and
motors run without torque pulsation, and the three return currents cancel,
so three wires do the work of six.

*Implementation:* :func:`engineeringkit.electrical_engineering.three_phase_power`.
The example checks the constant total power, the zero neutral current and
the rotating field.

*References:* N. Tesla, "A new system of alternate current motors and
transformers," *Trans. AIEE* 5 (1888) 308-327.

.. minigallery:: ../../examples/electrical_engineering/power/plot_01_tesla_three_phase.py

1893 -- Steinmetz's Phasors and Complex Impedance
-------------------------------------------------

Charles Proteus Steinmetz, at the International Electrical Congress in
Chicago, showed that a sinusoid :math:`A\cos(\omega t + \phi)` can be carried
as the complex number :math:`Ae^{j\phi}`. Differentiation becomes
multiplication by :math:`j\omega`, so an inductor is the impedance
:math:`j\omega L` and a capacitor :math:`1/j\omega C`. The differential
equations of an AC circuit become Ohm's and Kirchhoff's laws with complex
numbers. Steinmetz's method, taught at General Electric and in his
textbooks, made AC design routine.

*Implementation:* :func:`engineeringkit.electrical_engineering.inductor_impedance`,
:func:`engineeringkit.electrical_engineering.capacitor_impedance`,
:func:`engineeringkit.electrical_engineering.series_rlc_impedance`,
:func:`engineeringkit.electrical_engineering.phasor` and
:func:`engineeringkit.electrical_engineering.polar`, drawn with
:func:`engineeringkit.electrical_engineering.visualizers.plot_phasors`. The
example solves a series RLC circuit by phasors and checks it against a
time-domain simulation.

*References:* C. P. Steinmetz, "Complex quantities and their use in
electrical engineering," *Proc. International Electrical Congress, Chicago
1893* (New York: AIEE, 1894) 33-74.

.. minigallery:: ../../examples/electrical_engineering/ac_circuits/plot_02_steinmetz_phasors.py

1897 -- Lodge's Syntonic Tuning: Series Resonance
-------------------------------------------------

Spark transmitters radiated over a wide band, so every receiver heard every
station. Oliver Lodge's patent put an adjustable inductance in series with
the antenna capacitance at both ends of the link, tuning transmitter and
receiver to the same frequency, :math:`f_0 = 1/2\pi\sqrt{LC}`. At resonance
the reactances cancel and only the resistance limits the current. The
sharpness of the peak is the quality factor :math:`Q = \sqrt{L/C}/R`, and the
half-power bandwidth is :math:`f_0/Q`. Marconi's company bought the patent in
1912.

*Implementation:* :func:`engineeringkit.electrical_engineering.resonance_frequency`
and :class:`~engineeringkit.electrical_engineering.SeriesRLC` (its
``impedance``, ``quality_factor`` and ``step_response``). The example
measures the bandwidth and the rejection of a neighbouring station.

*References:* O. J. Lodge, "Improvements in syntonized telegraphy without
line wires," British Patent 11,575 (1897).

.. minigallery:: ../../examples/electrical_engineering/ac_circuits/plot_03_lodge_resonance.py

1918 -- Fortescue's Symmetrical Components
------------------------------------------

Charles LeGeyt Fortescue of Westinghouse showed that any three phasors are
the sum of three balanced sets: zero sequence (three equal phasors),
positive sequence (rotation a-b-c) and negative sequence (a-c-b). In a
symmetrical network the three sequences do not interact, so an unbalanced
problem splits into three single-phase ones. His 114-page paper became the
basis of fault analysis and protective relaying. A single line-to-ground
fault, for instance, connects the three sequence networks in series.

*Implementation:* :func:`engineeringkit.electrical_engineering.symmetrical_components`
and :func:`engineeringkit.electrical_engineering.phase_components`. The
example decomposes an unbalanced set and computes a line-to-ground fault.

*References:* C. L. Fortescue, "Method of symmetrical co-ordinates applied
to the solution of polyphase networks," *Trans. AIEE* 37 (1918)
1027-1140.

.. minigallery:: ../../examples/electrical_engineering/power/plot_02_fortescue_symmetrical_components.py

1926 -- Norton's Equivalent Circuit
-----------------------------------

Edward Lawry Norton of Bell Labs, in an internal memorandum on filter
design, and Hans Ferdinand Mayer of Siemens, in a paper on the vacuum-tube
amplifier, gave the dual of Thévenin's theorem in the same year. A linear
two-terminal network also acts as a current source :math:`V_{th}/R_{th}` in
parallel with the same resistance. Converting between the two forms (a
*source transformation*) is a standard hand-analysis step, and the Norton
form suits devices such as transistors that behave more like current
sources.

*Implementation:* :func:`engineeringkit.electrical_engineering.norton_from_thevenin`.
The example checks that the two forms are indistinguishable at the
terminals and combines two batteries in parallel.

*References:* E. L. Norton, "Design of finite networks for uniform
frequency characteristic," Bell Laboratories Technical Report
TM26-0-1860 (1926); H. F. Mayer, "Ueber das Ersatzschema der
Verstärkerröhre," *Telegraphen- und Fernsprech-Technik* 15 (1926) 335-337.

.. minigallery:: ../../examples/electrical_engineering/networks/plot_03_norton_equivalent.py

1927 -- Budeanu's Reactive Power and the Power Factor
-----------------------------------------------------

When the current lags the voltage, energy flows back and forth each half
cycle without doing work, yet still loads the cables. Constantin Budeanu
of Bucharest named this *reactive power*, :math:`Q = VI\sin\phi`, measured in
vars, a unit the IEC adopted in 1930. Real power :math:`P = VI\cos\phi`, the
apparent power :math:`S^2 = P^2 + Q^2` and the power factor :math:`P/S` became
the vocabulary of tariffs and power-factor correction. For distorted
waveforms Budeanu summed :math:`Q` over harmonics and called the remainder
distortion power. That definition is still debated, and IEEE 1459 (2000)
replaced it.

*Implementation:* :func:`engineeringkit.electrical_engineering.real_power`,
:func:`engineeringkit.electrical_engineering.reactive_power`,
:func:`engineeringkit.electrical_engineering.power_factor` and
:func:`engineeringkit.electrical_engineering.rms_sine`. The example
recovers the power triangle from waveforms and computes the distortion
power of a harmonic-rich load.

*References:* C. I. Budeanu, *Puissances réactives et fictives* (Bucharest:
Institut National Roumain pour l'Étude de l'Aménagement et de l'Utilisation
des Sources d'Énergie, 1927).

.. minigallery:: ../../examples/electrical_engineering/power/plot_03_budeanu_reactive_power.py

1939 -- Smith's Chart
---------------------

Phillip Hagar Smith of Bell Labs plotted the reflection coefficient
:math:`\Gamma = (Z_L - Z_0)/(Z_L + Z_0)` in the unit disc. Over it he drew
circles of constant normalised resistance and arcs of constant reactance.
Moving along a lossless line only rotates :math:`\Gamma` about the centre at
constant standing-wave ratio, so input impedances, standing waves and
matching networks could be read off a sheet of paper instead of worked out
with a slide rule. The chart outlived the slide rule: every vector network
analyser still displays it.

*Implementation:* :func:`engineeringkit.electrical_engineering.reflection_coefficient`,
:func:`engineeringkit.electrical_engineering.standing_wave_ratio` and
:func:`engineeringkit.electrical_engineering.input_impedance`, drawn with
:func:`engineeringkit.electrical_engineering.visualizers.plot_smith_chart`.
The example designs a single-stub match.

*References:* P. H. Smith, "Transmission line calculator," *Electronics* 12
(January 1939) 29-31.

.. minigallery:: ../../examples/electrical_engineering/transmission_lines/plot_02_smith_chart.py

1956 -- Ward and Hale's Digital Load Flow
-----------------------------------------

Power companies had studied their grids on AC network analysers, scale
models built from resistors, inductors and capacitors. J. B. Ward and
H. W. Hale of American Gas and Electric published the first practical
digital method. They described the network by its bus admittance matrix
and improved each bus voltage in turn from the power specified there, an
iteration now called Gauss-Seidel. It needs little memory but converges
slowly. Over-relaxation helps, and Newton-Raphson methods with sparse
matrices replaced it in the 1960s.

*Implementation:* :func:`engineeringkit.electrical_engineering.admittance_matrix`
and :func:`engineeringkit.electrical_engineering.gauss_seidel_load_flow`,
returning a :class:`~engineeringkit.electrical_engineering.LoadFlowResult`.
The example solves a five-bus network and checks it against SciPy's
nonlinear solver.

*References:* J. B. Ward and H. W. Hale, "Digital computer solution of
power-flow problems," *Trans. AIEE Part III* 75 (1956) 398-404.

.. minigallery:: ../../examples/electrical_engineering/power/plot_04_ward_hale_load_flow.py

See Also
--------

- :doc:`/api/electrical_engineering`
- :doc:`/guides/electrical_engineering`
- :doc:`/history/communications_breakthroughs`
- :doc:`/history/control_systems_breakthroughs`

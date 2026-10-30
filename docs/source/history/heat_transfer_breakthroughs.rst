Breakthroughs in Heat Transfer
==============================


.. include:: /_generated/nav/heat_transfer.rst

.. epigraph::

   "Heat, like gravity, penetrates every substance of the universe, its
   rays occupy all parts of space."
   -- Joseph Fourier, *Théorie analytique de la chaleur*, 1822

Heat moves by conduction, convection and radiation, and each became
computable at a different time. Conduction had its law and its equation by
1822. Radiation waited for the blackbody laws of the late nineteenth
century, and convection, which mixes heat transfer with fluid mechanics,
for the boundary-layer solutions and dimensionless correlations of the
twentieth. Phase change and heat
exchangers turned these pieces into the design rules of power stations,
engines and refrigerators. The entries trace the ideas behind
:mod:`engineeringkit.heat_transfer`. Each links to the implementing function
and to a gallery example devoted to it.

.. contents:: Timeline
   :local:
   :depth: 1

1701 -- Newton's Law of Cooling
-------------------------------

Isaac Newton wanted a temperature scale beyond the range of his linseed-oil
thermometer. He heated a block of iron red hot, set it in a steady wind,
and assumed that the heat it lost was proportional to its excess over the
air temperature. Equal intervals of time then multiply the excess by equal
factors, so the clock becomes a thermometer. The paper appeared anonymously
in the *Philosophical Transactions*. Today the law defines the convection
coefficient, :math:`\dot Q = hA(T_s - T_\infty)`, which hides all the fluid
mechanics in :math:`h`.

*Implementation:* :func:`engineeringkit.heat_transfer.convection_heat_rate`
and :func:`engineeringkit.heat_transfer.lumped_temperature`. The example
recovers the time constant from readings at equal intervals and checks the
energy balance of a cooling block.

*References:* [I. Newton], "Scala graduum caloris," *Phil. Trans. R. Soc.*
22 (1701) 824-829.

.. minigallery:: ../../examples/heat_transfer/convection/plot_01_newton_cooling.py

1804 -- Biot's Heated Bar and the Biot Number
---------------------------------------------

Jean-Baptiste Biot heated one end of iron bars and read thermometers along
them. The steady profile was set by a contest between conduction along the
bar and loss from its surface. Biot could not write a consistent
differential equation for it, because he treated conduction as an action
between distant slices; Fourier did in 1807. The ratio of surface
conductance to internal conductance, :math:`Bi = hL/k`, now carries Biot's
name and decides whether a body can be treated as one temperature.

*Implementation:* :func:`engineeringkit.heat_transfer.biot_number`. The
example compares the lumped model with the exact plane-wall solution at
Biot numbers of 0.01, 0.1 and 1.

*References:* J.-B. Biot, "Mémoire sur la propagation de la chaleur,"
*Bibliothèque britannique* 27 (1804) 310-329.

.. minigallery:: ../../examples/heat_transfer/transient/plot_01_biot_number.py

1822 -- Fourier's *Théorie analytique de la chaleur*
----------------------------------------------------

Joseph Fourier made the heat flux proportional to the temperature gradient,
:math:`q = -k\,\mathrm{d}T/\mathrm{d}x`, derived the heat equation from an
energy balance, and solved it with trigonometric series. The memoir of 1807
met objections from Lagrange and Laplace, and the full theory appeared as a
book in 1822. For steady flow through a slab the law is
:math:`\dot Q = kA\Delta T/L`, so layers behave like resistors in series.
Curved walls give logarithmic resistances and, with them, a critical radius
below which insulation increases heat loss.

*Implementation:* :func:`engineeringkit.heat_transfer.fourier_heat_rate`,
:func:`engineeringkit.heat_transfer.plane_wall_resistance`,
:func:`engineeringkit.heat_transfer.composite_wall`,
:func:`engineeringkit.heat_transfer.cylinder_resistance`,
:func:`engineeringkit.heat_transfer.critical_insulation_radius_cylinder`
and the conduction shape factors
(:func:`engineeringkit.heat_transfer.shape_factor_buried_cylinder`,
:func:`engineeringkit.heat_transfer.shape_factor_buried_sphere`,
:func:`engineeringkit.heat_transfer.shape_factor_eccentric_cylinders`,
:func:`engineeringkit.heat_transfer.shape_factor_parallel_cylinders`). The
example solves a furnace wall as a resistor chain, finds the critical
radius of wire insulation numerically, and works out the heat lost by a
buried pipe.

*References:* J. Fourier, *Théorie analytique de la chaleur* (Paris:
Firmin Didot, 1822).

.. minigallery:: ../../examples/heat_transfer/conduction/plot_01_fourier_law.py

1859 -- Kirchhoff's Law of Thermal Radiation
--------------------------------------------

Gustav Kirchhoff, working with Robert Bunsen on spectra, proved from the
second law that at each wavelength a body's emissivity equals its
absorptivity. If it did not, a body in an isothermal enclosure could heat
itself above its surroundings. The perfect absorber, the blackbody, became
the universal emitter whose spectrum depends on temperature alone, and
finding that spectrum occupied physics for the next forty years.

*Implementation:* :func:`engineeringkit.heat_transfer.radiation_heat_rate`
and :func:`engineeringkit.heat_transfer.radiation_coefficient`. The example
lets plates with :math:`\alpha = \varepsilon` and :math:`\alpha \ne \varepsilon`
come to equilibrium in an enclosure, and linearises radiation as a
heat-transfer coefficient.

*References:* G. Kirchhoff, "Über das Verhältniss zwischen dem
Emissionsvermögen und dem Absorptionsvermögen der Körper für Wärme und
Licht," *Ann. Phys. Chem.* 109 (1860) 275-301.

.. minigallery:: ../../examples/heat_transfer/radiation/plot_01_kirchhoff_emissivity.py

1862 -- Kelvin's Age of the Earth
---------------------------------

William Thomson treated the Earth as a semi-infinite solid that began molten
and has cooled from its surface ever since. The error-function solution
makes the surface gradient fall as :math:`1/\sqrt t`, so the gradient
measured in mines, about 1 °F per 50 ft, dates the Earth: about 98 million
years. Geologists and Darwin needed far longer. The calculation was right;
the model was not, because convection in the mantle (argued by John Perry
in 1895) and radioactive heat keep the interior hot.

*Implementation:*
:func:`engineeringkit.heat_transfer.semi_infinite_temperature`. The example
redoes Kelvin's arithmetic in SI units and checks the gradient against the
profile.

*References:* W. Thomson, "On the secular cooling of the earth," *Trans.
R. Soc. Edinburgh* 23 (1864) 157-169 (read 28 April 1862).

.. minigallery:: ../../examples/heat_transfer/transient/plot_02_kelvin_age_of_earth.py

1879-1884 -- Stefan and Boltzmann: Emission Grows as T⁴
-------------------------------------------------------

Josef Stefan noticed that John Tyndall's platinum wire radiated 11.7 times
as much at 1200 °C as at 525 °C, close to the ratio of the absolute
temperatures to the fourth power. Five years later his former student
Ludwig Boltzmann derived :math:`E_b = \sigma T^4` from thermodynamics and
Maxwell's radiation pressure. The value of :math:`\sigma` was explained
only in 1900, by integrating Planck's spectrum.

*Implementation:*
:func:`engineeringkit.heat_transfer.stefan_boltzmann_emissive_power` and
:func:`engineeringkit.heat_transfer.planck_spectral_emissive_power`. The
example reproduces Stefan's ratio and integrates Planck's law to recover
:math:`\sigma T^4`.

*References:* J. Stefan, "Über die Beziehung zwischen der Wärmestrahlung
und der Temperatur," *Sitzungsber. Akad. Wiss. Wien* 79 (1879) 391-428;
L. Boltzmann, "Ableitung des Stefan'schen Gesetzes, betreffend die
Abhängigkeit der Wärmestrahlung von der Temperatur aus der
electromagnetischen Lichttheorie," *Ann. Phys. Chem.* 22 (1884) 291-294.

.. minigallery:: ../../examples/heat_transfer/radiation/plot_02_stefan_boltzmann.py

1883 -- Graetz's Laminar Tube Flow
----------------------------------

Leo Graetz, measuring the thermal conductivity of liquids in Munich, ran
them through heated capillaries and needed the temperature of a fluid in
Poiseuille flow entering a tube whose wall is at a new temperature. He
solved it as a series of eigenfunctions of the radius. Far downstream only
the first survives, and the Nusselt number settles at 3.657; a uniform wall
flux gives 48/11. Near the inlet the thin thermal layer gives more, and the
profile is fully developed after the entry length :math:`0.05\,Re\,Pr\,D`.
Hausen (1943) fitted the mean over the entry region against the number
:math:`Re\,Pr\,D/L` that now carries Graetz's name; that paper was checked
only through Incropera et al.

*Implementation:* :func:`engineeringkit.heat_transfer.laminar_tube_nusselt`,
:func:`engineeringkit.heat_transfer.entry_length` and
:func:`engineeringkit.heat_transfer.hausen_nusselt`. The example marches
Graetz's problem down the tube by finite differences and recovers the
limit, the entry length and Hausen's average.

*References:* L. Graetz, "Über die Wärmeleitungsfähigkeit von
Flüssigkeiten," *Ann. Phys. Chem.* 18 (1883) 79-94, continued in 25 (1885)
337-357.

.. minigallery:: ../../examples/heat_transfer/convection/plot_04_graetz_tube_flow.py

1891 -- Stefan's Problem: A Moving Freezing Front
-------------------------------------------------

Josef Stefan asked how fast sea ice thickens in the polar winter. Latent
heat is released at the ice-water interface, whose position is part of the
unknown: the first moving-boundary problem. Stefan found that the thickness
grows as the square root of time. Franz Neumann had already given the exact
similarity solution in his Königsberg lectures,
:math:`s = 2\lambda\sqrt{\alpha t}` with
:math:`\lambda e^{\lambda^2}\operatorname{erf}\lambda = Ste/\sqrt\pi`. The
same mathematics governs casting, welding and frozen ground.

*Implementation:* :func:`engineeringkit.heat_transfer.neumann_stefan_front`
and :func:`engineeringkit.heat_transfer.stefan_number`. The example grows a
winter of sea ice, checks the energy balance at the front and compares
Stefan's quasi-steady estimate.

*References:* J. Stefan, "Über die Theorie der Eisbildung, insbesondere
über die Eisbildung im Polarmeere," *Ann. Phys. Chem.* 42 (1891) 269-286.

.. minigallery:: ../../examples/heat_transfer/phase_change/plot_01_stefan_problem.py

1893 -- Wien's Displacement Law
-------------------------------

Wilhelm Wien compressed blackbody radiation in an imaginary cavity with a
moving mirror. Doppler shifts and the second law forced the spectrum to
keep its shape: divided by :math:`T^5` it depends only on :math:`\lambda T`.
Its peak therefore sits at :math:`\lambda_{max}T = b = 2898` um K, in the
visible for the Sun and near 10 um for a person, which is why thermal
cameras work in that band.

*Implementation:*
:func:`engineeringkit.heat_transfer.blackbody_peak_wavelength`. The example
finds the peak of Planck's spectrum numerically and collapses spectra at
four temperatures onto one curve.

*References:* W. Wien, "Eine neue Beziehung der Strahlung schwarzer Körper
zum zweiten Hauptsatz der Wärmetheorie," *Sitzungsber. Preuss. Akad. Wiss.
Berlin* (1893) 55-62.

.. minigallery:: ../../examples/heat_transfer/radiation/plot_03_wien_displacement.py

1916 -- Nusselt's Film Condensation
-----------------------------------

Wilhelm Nusselt analysed steam condensing on a cold vertical wall. The
condensate forms a laminar film that drains under gravity; balancing weight
against viscous shear and latent heat against conduction across the film
gives a thickness growing as :math:`x^{1/4}` and an average coefficient
:math:`0.943[g\rho_l(\rho_l-\rho_v)k_l^3h_{fg}/(\mu_lL\Delta T)]^{1/4}`. It
was the first convection coefficient derived rather than measured, and is
still the basis of condenser design.

*Implementation:*
:func:`engineeringkit.heat_transfer.nusselt_film_condensation`. The example
integrates the local coefficient along the film and checks the mass and
energy balance of the condensate.

*References:* W. Nusselt, "Die Oberflächenkondensation des
Wasserdampfes," *Z. VDI* 60 (1916) 541-546, 569-575.

.. minigallery:: ../../examples/heat_transfer/phase_change/plot_02_nusselt_condensation.py

1921 -- Pohlhausen's Thermal Boundary Layer
-------------------------------------------

Blasius had reduced Prandtl's laminar boundary layer on a flat plate to
the equation :math:`f''' + \tfrac12 ff'' = 0` in 1908. Ernst Pohlhausen, in
the first volume of the *Zeitschrift für angewandte Mathematik und
Mechanik*, added the energy equation. With the same similarity variable it
becomes :math:`\theta'' + \tfrac12 Pr\,f\theta' = 0`, which integrates in
closed form once :math:`f` is known; the wall gradient is close to
:math:`0.332\,Pr^{1/3}`, so :math:`Nu_x = 0.332\,Re_x^{1/2}Pr^{1/3}`. With
the turbulent :math:`0.0296\,Re_x^{4/5}Pr^{1/3}` beyond transition, it
gives the flat-plate results every convection course starts from.

*Implementation:* :func:`engineeringkit.heat_transfer.flat_plate_nusselt`.
The example solves the Blasius and Pohlhausen equations, measures the
:math:`Pr^{1/3}` fit, and integrates the local coefficient over a plate
whose boundary layer turns turbulent.

*References:* E. Pohlhausen, "Der Wärmeaustausch zwischen festen Körpern
und Flüssigkeiten mit kleiner Reibung und kleiner Wärmeleitung," *Z.
Angew. Math. Mech.* 1 (1921) 115-121.

.. minigallery:: ../../examples/heat_transfer/convection/plot_05_pohlhausen_flat_plate.py

1922 -- Harper and Brown's Fin Efficiency
-----------------------------------------

Air-cooled aero engines shed their heat through cylinder fins. D. R.
Harper and W. B. Brown of the Bureau of Standards solved the fin equation
:math:`\theta'' = m^2\theta` for NACA, defined the fin efficiency, the heat
carried over that of an isothermal fin, and showed that a convecting tip is
accounted for by lengthening the fin by half its thickness. Their
efficiency curves are still how fins are sized.

*Implementation:* :func:`engineeringkit.heat_transfer.fin_parameter`,
:func:`engineeringkit.heat_transfer.fin_efficiency_adiabatic`,
:func:`engineeringkit.heat_transfer.fin_heat_rate_adiabatic` and
:func:`engineeringkit.heat_transfer.fin_heat_rate_convective_tip`. The example
solves the fin equation with a convecting tip numerically and checks the
corrected length.

*References:* D. R. Harper and W. B. Brown, "Mathematical equations for
heat conduction in the fins of air-cooled engines," NACA Report 158 (1922).

.. minigallery:: ../../examples/heat_transfer/conduction/plot_02_fin_efficiency.py

1928 -- Nusselt's Unit-Sphere View Factor
-----------------------------------------

The fraction of radiation leaving one diffuse surface that reaches another
is a fourfold integral. Wilhelm Nusselt reduced it to a drawing: project the
second surface onto a unit hemisphere, then down onto the base circle; the
covered fraction is the view factor. Closed forms for common geometries
followed, and Nusselt's construction lives on as the cosine-weighted
sampling of Monte Carlo ray tracers.

*Implementation:*
:func:`engineeringkit.heat_transfer.view_factor_parallel_plates` and
:func:`engineeringkit.heat_transfer.view_factor_coaxial_disks`. The example
turns the unit sphere into a Monte Carlo ray tracer and checks the closed
forms and reciprocity.

*References:* W. Nusselt, "Graphische Bestimmung des Winkelverhältnisses
bei der Wärmestrahlung," *Z. VDI* 72 (1928) 673.

.. minigallery:: ../../examples/heat_transfer/radiation/plot_04_nusselt_view_factor.py

1930 -- Dittus and Boelter: Turbulent Pipe Convection
-----------------------------------------------------

F. W. Dittus and L. M. K. Boelter at Berkeley correlated tests on
automobile radiator tubes with a power law in the Reynolds and Prandtl
numbers. McAdams's rounded version, :math:`Nu = 0.023Re^{0.8}Pr^n` with
:math:`n = 0.4` for heating and 0.3 for cooling, became the most quoted
equation in convection, though the original paper gave 0.0243 and 0.0265.
Colburn's analogy (1933) between friction and heat transfer reproduces the
0.023.

*Implementation:*
:func:`engineeringkit.heat_transfer.dittus_boelter_nusselt`. The example
sizes the coefficient for water in a tube, derives it from the Colburn
analogy and compares it with Gnielinski's correlation.

*References:* F. W. Dittus and L. M. K. Boelter, "Heat transfer in
automobile radiators of the tubular type," *Univ. Calif. Publ. Eng.* 2
(1930) 443-461; R. H. S. Winterton, "Where did the Dittus and Boelter
equation come from?", *Int. J. Heat Mass Transfer* 41 (1998) 809-810.

.. minigallery:: ../../examples/heat_transfer/convection/plot_02_dittus_boelter.py

1934 -- Nukiyama's Boiling Curve and Zuber's Critical Heat Flux
---------------------------------------------------------------

Shiro Nukiyama boiled water on electrically heated wires and found that the
heat flux reached a maximum a few tens of kelvin above saturation. At the
next increment of power the wire was blanketed by vapour and jumped to red
heat, or melted. The full curve, with its falling transition region and
minimum, was traced later by controlling temperature instead of power. In
1958 Novak Zuber explained the maximum as a Taylor-Helmholtz instability of
the vapour columns, :math:`q''_{max} = (\pi/24)h_{fg}\rho_v[\sigma g(\rho_l -
\rho_v)/\rho_v^2]^{1/4}`. The critical heat flux limits every boiler and
nuclear reactor core.

*Implementation:*
:func:`engineeringkit.heat_transfer.zuber_critical_heat_flux`. The example
sketches the boiling curve of water with standard correlations and finds
the temperature jump at constant power.

*References:* S. Nukiyama, "The maximum and minimum values of the heat Q
transmitted from metal to boiling water under atmospheric pressure," *J.
Soc. Mech. Eng. Japan* 37 (1934) 367-374 (transl. *Int. J. Heat Mass
Transfer* 9 (1966) 1419-1433); N. Zuber, "On the stability of boiling heat
transfer," *Trans. ASME* 80 (1958) 711-720.

.. minigallery:: ../../examples/heat_transfer/phase_change/plot_03_nukiyama_boiling_curve.py

1934 -- The Chilton-Colburn Analogy
-----------------------------------

Allan Colburn, at DuPont, showed in 1933 that turbulent heat transfer data
collapse when the Stanton number is multiplied by :math:`Pr^{2/3}`, and
that this :math:`j` factor equals half the friction factor. With Thomas
Chilton he extended it to mass transfer in 1934: equal :math:`j` factors
give the mass transfer coefficient from the heat transfer coefficient,
:math:`h_m = h/(\rho c_p Le^{2/3})`. One measurement then stands for three,
and evaporation, sublimation and absorption can be designed from heat
transfer data.

*Implementation:*
:func:`engineeringkit.heat_transfer.chilton_colburn_mass_transfer` and
:func:`engineeringkit.heat_transfer.lewis_number`. The example checks the
three :math:`j` factors of the laminar flat plate and computes a wet-bulb
temperature that does not depend on the air speed.

*References:* A. P. Colburn, "A method of correlating forced convection
heat transfer data and a comparison with fluid friction," *Trans. AIChE*
29 (1933) 174-210; T. H. Chilton and A. P. Colburn, "Mass transfer
(absorption) coefficients: prediction from data on heat transfer and fluid
friction," *Ind. Eng. Chem.* 26 (1934) 1183-1187.

.. minigallery:: ../../examples/heat_transfer/convection/plot_06_chilton_colburn_analogy.py

1940 -- Bowman, Mueller and Nagle: the LMTD Correction Factor
-------------------------------------------------------------

The log-mean temperature difference is exact only for pure counterflow or
parallel flow. R. A. Bowman, A. C. Mueller and W. M. Nagle collected the
analyses of multipass and crossflow exchangers made in the 1930s, Nagle's
and Bowman's among them, into one correction factor :math:`F` on the
counterflow LMTD, charted against the effectiveness :math:`P` of one stream
and the capacity ratio :math:`R`. Their charts are still in every
exchanger handbook. The effectiveness-NTU relations of the same
arrangements carry the same information in another form.

*Implementation:*
:func:`engineeringkit.heat_transfer.lmtd_correction_factor`,
:func:`engineeringkit.heat_transfer.exchanger_effectiveness_shell_tube` and
:func:`engineeringkit.heat_transfer.exchanger_effectiveness_crossflow`. The
example redraws the one-shell chart, checks the F-LMTD and
effectiveness-NTU methods against each other, and rates a two-shell oil
cooler.

*References:* R. A. Bowman, A. C. Mueller and W. M. Nagle, "Mean
temperature difference in design," *Trans. ASME* 62 (1940) 283-294.

.. minigallery:: ../../examples/heat_transfer/exchangers/plot_02_bowman_correction_factor.py

1947 -- Heisler's Transient Conduction Charts
---------------------------------------------

The cooling of a slab, cylinder or sphere is a Fourier series whose
eigenvalues solve a transcendental equation, :math:`\zeta\tan\zeta = Bi`
for the slab. M. P. Heisler, working on induction heating, kept only the
first term, valid once the Fourier number exceeds about 0.2, and charted the
centre temperature against :math:`Fo` for each :math:`1/Bi`. The Heisler
charts were the standard tool for transient conduction until computers made
the series trivial.

*Implementation:*
:func:`engineeringkit.heat_transfer.plane_wall_transient`,
:func:`engineeringkit.heat_transfer.cylinder_transient`,
:func:`engineeringkit.heat_transfer.sphere_transient` and
:func:`engineeringkit.heat_transfer.fourier_number`. The example redraws the
centre chart, measures the one-term error, times the quench of a steel
plate and compares the slab, cylinder and sphere.

*References:* M. P. Heisler, "Temperature charts for induction and
constant-temperature heating," *Trans. ASME* 69 (1947) 227-236.

.. minigallery:: ../../examples/heat_transfer/transient/plot_03_heisler_charts.py

1947 -- Crank and Nicolson: Numerical Conduction
------------------------------------------------

The explicit finite-difference scheme for the heat equation is stable only
for a mesh Fourier number :math:`\alpha\Delta t/\Delta x^2 \le 1/2`, which
forces tiny time steps; the fully implicit scheme is always stable but only
first-order accurate. John Crank and Phyllis Nicolson averaged the two and
obtained a scheme that is both unconditionally stable and second order in
time. It remains the default for diffusion problems, from heat conduction
to option pricing.

*Implementation:* :func:`engineeringkit.heat_transfer.conduction_1d` and,
for steady two-dimensional conduction with the same finite differences,
:func:`engineeringkit.heat_transfer.conduction_2d`. The example measures the
order of accuracy of the implicit and Crank-Nicolson schemes, watches the
explicit one blow up above :math:`r = 1/2`, and checks a 2-D plate against
its Fourier series.

*References:* J. Crank and P. Nicolson, "A practical method for numerical
evaluation of solutions of partial differential equations of the
heat-conduction type," *Proc. Cambridge Phil. Soc.* 43 (1947) 50-67.

.. minigallery:: ../../examples/heat_transfer/conduction/plot_03_crank_nicolson.py

1952 -- Rohsenow's Nucleate Boiling Correlation
-----------------------------------------------

In nucleate boiling the flux climbs steeply with the wall superheat, and
the data for different liquids and surfaces seemed to share nothing.
Warren Rohsenow at MIT treated each departing bubble as the agent of a
convection process and correlated them all with dimensionless groups, one
empirical surface-fluid constant :math:`C_{s,f}` and a cube law,
:math:`q'' \propto \Delta T_e^3`. The work circulated as an MIT report in
1951 and appeared in the *Transactions of the ASME* in 1952. Errors in the
flux of up to 100% are common, but because of the cube they mean only about
25% in the superheat.

*Implementation:*
:func:`engineeringkit.heat_transfer.rohsenow_boiling_heat_flux`. The
example boils water in a pan, shows the cube law and the effect of the
surface, and finds where the curve meets the critical heat flux.

*References:* W. M. Rohsenow, "A method of correlating heat-transfer data
for surface boiling of liquids," *Trans. ASME* 74 (1952) 969-976.

.. minigallery:: ../../examples/heat_transfer/phase_change/plot_04_rohsenow_nucleate_boiling.py

1955 -- Kays and London: Effectiveness-NTU
------------------------------------------

The log-mean temperature difference sizes an exchanger when all four end
temperatures are known, but rating a given exchanger with it needs
iteration. William Kays and Alexander London, testing compact exchangers at
Stanford, organised the results around the effectiveness
:math:`\varepsilon = \dot Q/\dot Q_{max}` as a function of
:math:`NTU = UA/C_{min}` and the capacity ratio, and published charts for
every flow arrangement in *Compact Heat Exchangers*.

*Implementation:*
:func:`engineeringkit.heat_transfer.exchanger_effectiveness_counterflow`,
:func:`engineeringkit.heat_transfer.exchanger_effectiveness_parallel`,
:func:`engineeringkit.heat_transfer.log_mean_temperature_difference`,
:func:`engineeringkit.heat_transfer.heat_exchanger_rate` and
:func:`engineeringkit.heat_transfer.overall_u`. The example rates an oil
cooler by both methods and integrates the stream temperatures along it.

*References:* W. M. Kays and A. L. London, *Compact Heat Exchangers*
(1955).

.. minigallery:: ../../examples/heat_transfer/exchangers/plot_01_kays_london_ntu.py

1956 -- Oppenheim's Radiation Network
-------------------------------------

Between gray surfaces, radiation is reflected back and forth, so every
surface's exchange depends on every other. Antoni Oppenheim, at Berkeley,
drew it as an electrical circuit. Each surface's blackbody emissive power
drives current through a surface resistance
:math:`(1 - \varepsilon)/(\varepsilon A)` to a node at its radiosity, and
the radiosity nodes are joined by space resistances :math:`1/(A_iF_{ij})`.
The network is linear in the radiosities. Insulated walls become floating
nodes, radiation shields add resistances in series, and the view factors
must obey reciprocity, :math:`A_iF_{ij} = A_jF_{ji}`, and summation,
:math:`\sum_jF_{ij} = 1`.

*Implementation:* :func:`engineeringkit.heat_transfer.radiation_enclosure`,
:func:`engineeringkit.heat_transfer.two_surface_radiation`,
:func:`engineeringkit.heat_transfer.radiation_shield_heat_flux` and
:func:`engineeringkit.heat_transfer.view_factor_reciprocity`. The example
builds a furnace's view factors from the two rules, solves its network by
hand and by matrix, and adds foil shields to a cryogenic line.

*References:* A. K. Oppenheim, "Radiation analysis by the network method,"
*Trans. ASME* 78 (1956) 725-735.

.. minigallery:: ../../examples/heat_transfer/radiation/plot_05_oppenheim_network.py

1972 -- Žukauskas's Tube Banks in Crossflow
-------------------------------------------

Boilers, air heaters and the shells of most exchangers are banks of tubes
in crossflow. The first row sees the free stream; deeper rows sit in the
wakes of those ahead, whose turbulence raises their heat transfer. Algirdas
Žukauskas, at the Lithuanian Academy of Sciences in Kaunas, gathered his
laboratory's measurements on aligned and staggered banks into one form,
:math:`\bar{Nu} = C_2CRe_{max}^mPr^{0.36}(Pr/Pr_s)^{1/4}`, with constants for
each layout and Reynolds range and a correction for banks of fewer than 20
rows. The constants used here are those tabulated by Incropera et al.; the
1972 review was located but its tables were not checked against them.

*Implementation:* :func:`engineeringkit.heat_transfer.zukauskas_nusselt`
and :func:`engineeringkit.heat_transfer.tube_bank_max_velocity`. The
example sizes an air heater, follows the air temperature through the bank,
and compares row counts and layouts.

*References:* A. Žukauskas, "Heat transfer from tubes in crossflow," *Adv.
Heat Transfer* 8 (1972) 93-160.

.. minigallery:: ../../examples/heat_transfer/convection/plot_07_zukauskas_tube_banks.py

1975 -- Churchill and Chu: Natural Convection over All Rayleigh Numbers
-----------------------------------------------------------------------

Free convection on a vertical plate was long described piecewise: a laminar
:math:`Ra^{1/4}` law, a turbulent :math:`Ra^{1/3}` law and a gap between
them. Stuart Churchill and Humbert Chu fitted one smooth expression from
creeping flow to turbulent plumes, using Churchill's method of blending
asymptotes. Their correlation is the standard for walls, radiators and
electronics enclosures.

*Implementation:*
:func:`engineeringkit.heat_transfer.churchill_chu_nusselt` and
:func:`engineeringkit.heat_transfer.rayleigh_number`. The example checks
both asymptotes and the laminar boundary-layer solution.

*References:* S. W. Churchill and H. H. S. Chu, "Correlating equations for
laminar and turbulent free convection from a vertical plate," *Int. J. Heat
Mass Transfer* 18 (1975) 1323-1329.

.. minigallery:: ../../examples/heat_transfer/convection/plot_03_churchill_chu.py

1977 -- Churchill and Bernstein: the Cylinder in Crossflow
----------------------------------------------------------

Heat transfer from a cylinder in crossflow was given by power laws fitted
range by range; Hilpert's of 1933 has five, each with its own constants.
Stuart Churchill and M. Bernstein applied the blending of asymptotes that
Churchill had used for free convection and wrote one expression from
creeping flow past the drag crisis. They described it as a slight
underestimate of most of the data it was fitted to.

*Implementation:*
:func:`engineeringkit.heat_transfer.churchill_bernstein_nusselt`. The
example compares it with Hilpert's table and calibrates a hot-wire
anemometer, recovering King's law of 1914.

*References:* S. W. Churchill and M. Bernstein, "A correlating equation for
forced convection from gases and liquids to a circular cylinder in
crossflow," *J. Heat Transfer* 99 (1977) 300-306.

.. minigallery:: ../../examples/heat_transfer/convection/plot_08_churchill_bernstein_cylinder.py

See Also
--------

- :doc:`/api/heat_transfer`
- :doc:`/guides/heat_transfer`
- :doc:`/history/thermodynamics_breakthroughs`
- :doc:`/history/fluid_engineering_breakthroughs`

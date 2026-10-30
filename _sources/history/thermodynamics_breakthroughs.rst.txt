Breakthroughs in Thermodynamics
===============================


.. include:: /_generated/nav/thermodynamics.rst

.. epigraph::

   "The energy of the universe is constant. The entropy of the universe
   tends to a maximum."
   -- Rudolf Clausius, 1865

Thermodynamics began as the theory of the steam engine and became a science
of energy itself. Carnot's question of how much work heat can give led,
through Joule's experiments, to the first law, and through Clausius and
Kelvin to the second law and entropy. Engineers turned the laws into cycles
for every kind of engine, refrigerator and heat pump. Equations of state and
vapour-pressure laws described the working fluids and their mixtures, and
exergy put all losses in one currency. The entries trace the ideas behind
:mod:`engineeringkit.thermodynamics`. Each links to the implementing
function and to a gallery example devoted to it.

.. contents:: Timeline
   :local:
   :depth: 1

1801 -- Dalton's Law of Partial Pressures
-----------------------------------------

John Dalton, a Manchester schoolteacher, read four essays on mixed gases,
steam and evaporation to the city's Literary and Philosophical Society in
October 1801. He argued that each gas in a mixture behaves as if it alone
filled the vessel, so the pressure of a mixture is the sum of the partial
pressures, :math:`p = \sum p_i`. The rule became the basis of every
ideal-gas mixture calculation, from moist air to combustion products, and
led Dalton towards his atomic theory.

*Implementation:* :func:`engineeringkit.thermodynamics.partial_pressures`,
:func:`engineeringkit.thermodynamics.mole_fractions`,
:func:`engineeringkit.thermodynamics.mass_fractions`,
:func:`engineeringkit.thermodynamics.mixture_molar_mass` and
:func:`engineeringkit.thermodynamics.entropy_of_mixing`. The example builds
dry air from its components and computes the entropy of mixing.

*References:* J. Dalton, "Experimental essays on the constitution of mixed
gases; on the force of steam or vapour from water and other liquids in
different temperatures, both in a Torricellian vacuum and in air; on
evaporation; and on the expansion of gases by heat," *Mem. Lit. Phil. Soc.
Manchester* 5 (1802) 535-602. The volume and pages follow the historical
literature; the volume itself could not be checked online.

.. minigallery:: ../../examples/thermodynamics/properties/plot_07_dalton_partial_pressures.py

1816 -- Stirling's Regenerator
------------------------------

Robert Stirling, a young minister of the Church of Scotland, patented a
hot-air engine whose working air passed back and forth through an
"economiser" of thin metal between the hot and cold spaces. The air left its
heat in the regenerator as it cooled and took it back as it returned, so only
the isothermal heat exchanges crossed the engine boundary. The ideal Stirling
cycle of two isotherms and two constant-volume processes then reaches
Carnot's efficiency, eight years before Carnot stated the bound.

*Implementation:* :func:`engineeringkit.thermodynamics.stirling_efficiency`.
The example draws the cycle and shows how much efficiency an imperfect
regenerator loses.

*References:* R. Stirling, "Improvements for diminishing the consumption of
fuel, and in particular an engine capable of being applied to the moving
(of) machinery on a principle entirely new," British Patent 4081 (1816). The
patent number and date follow the historical literature; the patent record
could not be checked online.

.. minigallery:: ../../examples/thermodynamics/cycles/plot_07_stirling_regenerator.py

1824 -- Carnot's *Réflexions*
-----------------------------

Sadi Carnot, a young French engineer officer, asked whether the work of a
steam engine had a limit. He compared heat falling from a hot body to a cold
one with water falling over a wheel, and argued that a reversible engine
between two temperatures is the best possible, whatever its working
substance. If any engine did better, it could drive a reversible engine
backwards and move heat from cold to hot for nothing. The little book was
almost ignored until Clapeyron and then Kelvin took it up. In modern form the
bound is :math:`\eta = 1 - T_c/T_h`.

*Implementation:* :func:`engineeringkit.thermodynamics.carnot_efficiency`,
:func:`engineeringkit.thermodynamics.heat_engine_work` and
:func:`engineeringkit.thermodynamics.thermal_efficiency`. The example runs
Carnot's argument with numbers and plots the bound against source
temperature.

*References:* S. Carnot, *Réflexions sur la puissance motrice du feu*
(Paris: Bachelier, 1824).

.. minigallery:: ../../examples/thermodynamics/cycles/plot_01_carnot_efficiency.py

1834 -- Clapeyron's Ideal Gas Law
---------------------------------

Émile Clapeyron restated Carnot's argument with a cycle drawn on pressure
and volume axes, the indicator diagram of the steam engineer. To compute the
work of a gas cycle he combined Boyle and Mariotte's law with Gay-Lussac's
law of expansion into :math:`pv = R(267 + t)`. With Regnault's better
coefficient and absolute temperature it became :math:`pV = nRT`, the
equation of state against which every real fluid is still measured.

*Implementation:* :func:`engineeringkit.thermodynamics.ideal_gas_pressure`
and :func:`engineeringkit.thermodynamics.ideal_gas_temperature`. The example
extrapolates a constant-volume gas thermometer to the absolute zero.

*References:* É. Clapeyron, "Mémoire sur la puissance motrice de la
chaleur," *J. École Polytechnique* 14 (1834) 153-190.

.. minigallery:: ../../examples/thermodynamics/properties/plot_01_clapeyron_ideal_gas.py

1834-1850 -- The Clausius-Clapeyron Relation
--------------------------------------------

In the same memoir Clapeyron applied Carnot's reasoning to an infinitesimal
cycle of boiling and condensing. The slope of the vapour-pressure curve came
out fixed by the latent heat and the change of volume,
:math:`\mathrm{d}p/\mathrm{d}T = h_{fg}/[T(v_g - v_f)]`. Clausius rederived
it in 1850 from the new mechanical theory of heat. Neglecting the liquid
volume and treating the vapour as an ideal gas, he integrated it to make
:math:`\ln p` linear in :math:`1/T`.

*Implementation:*
:func:`engineeringkit.thermodynamics.clausius_clapeyron_pressure`. The
example checks the exact slope against IAPWS-IF97 steam and finds the boiling
point on Everest.

*References:* É. Clapeyron (1834), as above; R. Clausius, "Über die
bewegende Kraft der Wärme," *Ann. Phys.* 79 (1850) 368-397, 500-524.

.. minigallery:: ../../examples/thermodynamics/properties/plot_02_clausius_clapeyron.py

1834 -- Perkins's Vapour-Compression Refrigerator
-------------------------------------------------

Jacob Perkins, an American inventor working in London, patented a closed
machine that evaporated a volatile liquid at low pressure to absorb heat,
compressed the vapour and condensed it under pressure, then let the liquid
back through a valve. Every domestic refrigerator and air conditioner still
works this way. The coefficient of performance, heat removed per unit of
work, is bounded by Carnot's :math:`T_c/(T_h - T_c)`.

*Implementation:* :func:`engineeringkit.thermodynamics.refrigeration_cop`,
:func:`engineeringkit.thermodynamics.carnot_cop_refrigerator`, and
:func:`engineeringkit.thermodynamics.vapor_compression_cycle` with
:class:`engineeringkit.thermodynamics.VaporCompressionResult` for the cycle on
a Peng-Robinson refrigerant. The examples run the ideal cycle with water as
refrigerant on IF97 properties and a textbook R-134a refrigerator.

*References:* J. Perkins, "Apparatus and means for producing ice, and in
cooling fluids," British Patent 6662 (1834).

.. minigallery:: ../../examples/thermodynamics/refrigeration/plot_01_perkins_vapor_compression.py ../../examples/thermodynamics/refrigeration/plot_03_perkins_r134a_cycle.py

1840 -- Hess's Law and the Heat of Combustion
---------------------------------------------

Germain Henri Hess, professor of chemistry in St Petersburg, measured the heat
of reactions run in one step and in several. The total was always the same:
the heat of a reaction depends only on its initial and final states. With
Lavoisier's conservation of mass this gives the energy released by burning a
fuel from tabulated enthalpies of formation. The same balance, with nothing
lost to the surroundings, fixes the temperature a flame can reach.

*Implementation:* :func:`engineeringkit.thermodynamics.heating_value`,
:func:`engineeringkit.thermodynamics.combustion_products`,
:func:`engineeringkit.thermodynamics.stoichiometric_air_fuel_ratio`,
:func:`engineeringkit.thermodynamics.fuel_molar_mass` and
:func:`engineeringkit.thermodynamics.adiabatic_flame_temperature`, with
product enthalpies from NIST's Shomate equations
(:mod:`engineeringkit.thermodynamics.utils.shomate`). The example finds the
heating values of three fuels and the flame temperature of methane.

*References:* G. H. Hess, "Thermochemische Untersuchungen," *Ann. Phys.*
126 (1840) 385-404 (Poggendorff's vol. 50), :doi:`10.1002/andp.18401260620`.

.. minigallery:: ../../examples/thermodynamics/combustion/plot_01_hess_heat_of_combustion.py

1843-1850 -- Joule and the First Law
------------------------------------

James Prescott Joule, a Manchester brewer's son, measured the heat produced
by electric currents, by compressing air, and finally by falling weights
that stirred water with a paddle wheel. Each time the same work gave the same
heat: 772 foot-pounds per British thermal unit, about 4.16 J/cal. Heat and work
became two forms of one conserved quantity, and Clausius wrote the balance
as :math:`Q = \Delta U + W`.

*Implementation:* :func:`engineeringkit.thermodynamics.first_law_heat`,
:func:`engineeringkit.thermodynamics.specific_internal_energy_ideal_gas` and
:func:`engineeringkit.thermodynamics.specific_enthalpy_ideal_gas`; for open
systems, :func:`engineeringkit.thermodynamics.steady_flow_work`,
:func:`engineeringkit.thermodynamics.nozzle_exit_velocity`,
:func:`engineeringkit.thermodynamics.mixing_chamber`,
:func:`engineeringkit.thermodynamics.throttle_steam` and the isentropic
efficiencies of turbines, compressors and nozzles. The examples redo the
paddle-wheel balance, convert Joule's number to SI and apply the balance to a
steam turbine, a nozzle, a throttle and a feedwater heater.

*References:* J. P. Joule, "On the calorific effects of magneto-electricity,
and on the mechanical value of heat," *Phil. Mag.* 23 (1843) 263-276, 347-355,
435-443; "On the mechanical equivalent of heat," *Phil. Trans. R. Soc.* 140
(1850) 61-82.

.. minigallery:: ../../examples/thermodynamics/laws/plot_01_joule_first_law.py ../../examples/thermodynamics/laws/plot_04_joule_steady_flow_devices.py

1852 -- Kelvin's Heat Pump
--------------------------

William Thomson saw that a reversed heat engine could warm a building with
heat taken from the colder outside air, delivering more heat than the work
spent on it. He proposed an air machine for the purpose and argued that it
would heat with less fuel than burning that fuel directly. The heat pump's
coefficient of performance is bounded by :math:`T_h/(T_h - T_c)`, one more
than the refrigerator's.

*Implementation:* :func:`engineeringkit.thermodynamics.heat_pump_cop` and
:func:`engineeringkit.thermodynamics.carnot_cop_heat_pump`. The example
works out Kelvin's air cycle as a reversed Brayton cycle and repeats his
economy argument.

*References:* W. Thomson, "On the economy of the heating or cooling of
buildings by means of currents of air," *Proc. Phil. Soc. Glasgow* 3 (1852)
269-272.

.. minigallery:: ../../examples/thermodynamics/refrigeration/plot_02_kelvin_heat_pump.py

1852 -- The Joule-Thomson Effect
--------------------------------

Joule and Thomson forced gases steadily through a porous plug and found
them slightly cooler on the far side. Air cooled by about a quarter of a
kelvin per atmosphere and carbon dioxide by about one, while hydrogen warmed.
The throttled flow keeps its enthalpy, so the effect measures how far a
real gas departs from ideal. Below an inversion temperature throttling
cools the gas, the basis of Linde's air liquefier of 1895.

*Implementation:*
:func:`engineeringkit.thermodynamics.joule_thomson_coefficient_vdw`. The
example compares nitrogen, carbon dioxide and hydrogen and integrates a
throttle valve.

*References:* J. P. Joule and W. Thomson, "On the thermal effects
experienced by air in rushing through small apertures," *Phil. Mag.* 4
(1852) 481-492.

.. minigallery:: ../../examples/thermodynamics/properties/plot_03_joule_thomson.py

1854-1865 -- Clausius's Inequality and Entropy
----------------------------------------------

Rudolf Clausius showed in 1854 that around any cycle the heat received,
divided by the temperature at which it crosses the boundary, sums to zero
for a reversible cycle and is negative otherwise. So :math:`\delta Q_{rev}/T`
is the change of a property of state, which in 1865 he named entropy, from
the Greek for transformation. The second law became a statement about
entropy: in an isolated system it never decreases.

*Implementation:* :func:`engineeringkit.thermodynamics.clausius_inequality`,
:func:`engineeringkit.thermodynamics.entropy_change_isothermal_ideal_gas`,
:func:`engineeringkit.thermodynamics.reversible_isothermal_work_ideal_gas`,
:func:`engineeringkit.thermodynamics.isentropic_pressure_volume` and
:func:`engineeringkit.thermodynamics.isentropic_temperature_ratio`. The
example builds an ideal-gas Carnot cycle and makes it irreversible.

*References:* R. Clausius, "Über eine veränderte Form des zweiten
Hauptsatzes der mechanischen Wärmetheorie," *Ann. Phys.* 93 (1854) 481-506;
"Über verschiedene für die Anwendung bequeme Formen der Hauptgleichungen der
mechanischen Wärmetheorie," *Ann. Phys.* 125 (1865) 353-400.

.. minigallery:: ../../examples/thermodynamics/laws/plot_02_clausius_entropy.py

1859 -- Rankine's Steam Cycle
-----------------------------

William Rankine's *Manual of the Steam Engine* was the first textbook of
engineering thermodynamics. It analysed the steam plant as a cycle of pump,
boiler, engine and condenser. Pumping the water as a liquid costs about 1% of
the engine's work, and the condenser lets the steam expand to a vacuum.
With superheat and reheat added, the Rankine cycle still produces most of the
world's electricity, from coal, nuclear and solar thermal plants.

*Implementation:*
:func:`engineeringkit.thermodynamics.rankine_ideal_efficiency`, with states
from :func:`engineeringkit.thermodynamics.steam_properties`, and
:func:`engineeringkit.thermodynamics.rankine_cycle` with
:class:`engineeringkit.thermodynamics.RankineResult` for plants with
superheat, reheat, an open feedwater heater and real turbines and pumps. The
examples compute a 10 MPa, 500 °C plant on the T-s diagram and reproduce
three textbook plants.

*References:* W. J. M. Rankine, *A Manual of the Steam Engine and Other
Prime Movers* (London: Griffin, 1859).

.. minigallery:: ../../examples/thermodynamics/cycles/plot_02_rankine_cycle.py ../../examples/thermodynamics/cycles/plot_08_rankine_reheat_regeneration.py

1871 -- Maxwell's Thermodynamic Relations
-----------------------------------------

In his textbook *Theory of Heat* James Clerk Maxwell drew an infinitesimal
cycle on the pressure-volume diagram. From it he read four relations between
the derivatives of pressure, volume, temperature and entropy, such as
:math:`(\partial s/\partial p)_T = -(\partial v/\partial T)_p`. They let an
engineer compute entropy and enthalpy changes, the difference of the specific
heats and the Joule-Thomson coefficient of any fluid from measurable
:math:`p`-:math:`v`-:math:`T` data. That is how every steam and refrigerant
table is built.

*Implementation:* :func:`engineeringkit.thermodynamics.cp_minus_cv` and
:func:`engineeringkit.thermodynamics.joule_thomson_coefficient`. The example
checks a Maxwell relation and both consequences on IAPWS-IF97 steam.

*References:* J. C. Maxwell, *Theory of Heat* (London: Longmans, Green,
1871), Ch. 9, "Four thermodynamic relations". The relations appear in
Chapter 9 of the later editions, which were checked; the 1871 first edition
could not be.

.. minigallery:: ../../examples/thermodynamics/laws/plot_03_maxwell_relations.py

1872 -- Brayton's Constant-Pressure Cycle
-----------------------------------------

George Brayton's "Ready Motor" compressed air, burned fuel in it at constant
pressure and expanded the hot gas in a second cylinder. Joule had proposed the
same cycle for an air engine. The engine lost out to Otto's, but the cycle of
isentropic compression, constant-pressure heating, isentropic expansion and
constant-pressure cooling is the cycle of every gas turbine and jet engine.
Its ideal efficiency depends only on the pressure ratio.

*Implementation:* :func:`engineeringkit.thermodynamics.brayton_efficiency`,
and :func:`engineeringkit.thermodynamics.brayton_cycle` with
:class:`engineeringkit.thermodynamics.BraytonResult` for gas turbines with
real components, regeneration, intercooling and reheat. The examples rebuild
the efficiency from the four state temperatures and follow the staged,
regenerated cycle towards Ericsson's and Carnot's.

*References:* G. B. Brayton, "Improvement in gas-engines," U.S. Patent
125,166 (1872).

.. minigallery:: ../../examples/thermodynamics/cycles/plot_03_brayton_cycle.py ../../examples/thermodynamics/cycles/plot_09_brayton_regeneration_intercooling.py

1873 -- Van der Waals's Equation of State
-----------------------------------------

Andrews had shown in 1869 that carbon dioxide has a critical temperature
above which it cannot be liquefied. In his Leiden thesis Johannes van der
Waals corrected the ideal gas for the volume of the molecules and their
mutual attraction, :math:`p = RT/(v - b) - a/v^2`. One equation then covers
gas and liquid and predicts the critical point. Maxwell's equal-area rule
(1875) gives the vapour pressure. Van der Waals received the 1910 Nobel Prize
in Physics.

*Implementation:*
:func:`engineeringkit.thermodynamics.van_der_waals_pressure` and
:func:`engineeringkit.thermodynamics.van_der_waals_critical_point`. The
example draws carbon dioxide's isotherms with Maxwell's construction.

*References:* J. D. van der Waals, *Over de continuïteit van den gas- en
vloeistoftoestand*, thesis, Leiden (1873).

.. minigallery:: ../../examples/thermodynamics/properties/plot_04_van_der_waals.py

1876 -- Otto's Four-Stroke Cycle
--------------------------------

Nicolaus Otto's engine compressed its charge of gas and air before igniting
it, a principle Beau de Rochas had published in 1862. Its four strokes of
intake, compression, power and exhaust are idealized as the air-standard
Otto cycle, with heat added and rejected at constant volume. The efficiency
:math:`1 - r^{1-\gamma}` depends only on the compression ratio, and
knock sets the limit on that ratio.

*Implementation:* :func:`engineeringkit.thermodynamics.otto_efficiency` and
:class:`engineeringkit.thermodynamics.OttoCycle`. The example draws the
cycle on the p-v diagram and checks the efficiency from the four states.

*References:* N. A. Otto, "Gas motor engine," U.S. Patent 194,047 (1877).

.. minigallery:: ../../examples/thermodynamics/cycles/plot_04_otto_cycle.py

1888 -- Antoine's Vapour-Pressure Equation
------------------------------------------

Measured vapour pressures bend away from the straight Clausius-Clapeyron
line because the latent heat falls with temperature. Charles Antoine added a
third constant, :math:`\log_{10} p = A - B/(C + T)`, which fits most liquids
within a fraction of a percent over a wide range. Tables of Antoine constants
remain a standard tool of chemical engineering.

*Implementation:*
:func:`engineeringkit.thermodynamics.antoine_vapor_pressure`. The example
fits Antoine and Clausius-Clapeyron to IF97 water and compares their errors.

*References:* C. Antoine, "Tensions des vapeurs; nouvelle relation entre
les tensions et les températures," *C. R. Acad. Sci.* 107 (1888) 681-684,
778-780, 836-837.

.. minigallery:: ../../examples/thermodynamics/properties/plot_05_antoine_equation.py

1889-1905 -- Gouy and Stodola: Lost Work
----------------------------------------

Georges Gouy defined the useful energy of a system in surroundings at
:math:`T_0` and showed that every irreversibility destroys work equal to
:math:`T_0` times the entropy it generates. Aurel Stodola found the same
result for steam turbines. The Gouy-Stodola theorem lets an engineer rank
the losses of a plant, from heat crossing temperature differences to
throttling and friction, in one currency.

*Implementation:*
:func:`engineeringkit.thermodynamics.exergy_destruction`. The example prices
heat transfer, a throttling valve and an imperfect turbine.

*References:* G. Gouy, "Sur l'énergie utilisable," *J. Phys. Théor. Appl.*
8 (1889) 501-518; A. Stodola, *Die Dampfturbinen*, 3rd ed. (Berlin:
Springer, 1905).

.. minigallery:: ../../examples/thermodynamics/exergy/plot_01_gouy_stodola.py

1893 -- Diesel's Compression-Ignition Cycle
-------------------------------------------

Rudolf Diesel set out to approach Carnot's efficiency by compressing air
alone until injected fuel ignites by itself. With no fuel present during
compression the engine cannot knock, so it can use compression ratios of 15
to 22. The air-standard Diesel cycle adds heat at constant pressure up to a
cutoff ratio. At the same compression ratio it is less efficient than Otto's
cycle, but its higher ratios make it more efficient in practice.

*Implementation:* :func:`engineeringkit.thermodynamics.diesel_efficiency`
and :func:`engineeringkit.thermodynamics.dual_efficiency`, the
limited-pressure cycle of fast diesels with heat added first at constant
volume. The examples compare Diesel and Otto cycles at their practical
compression ratios and split the heat addition of the dual cycle.

*References:* R. Diesel, *Theorie und Konstruktion eines rationellen
Wärmemotors* (Berlin: Springer, 1893).

.. minigallery:: ../../examples/thermodynamics/cycles/plot_05_diesel_cycle.py ../../examples/thermodynamics/cycles/plot_10_diesel_dual_cycle.py

1911 -- Carrier's Psychrometric Chart
-------------------------------------

Willis Carrier controlled the humidity of a Brooklyn printing plant in 1902
by cooling the air below its dew point and reheating it. His paper of 1911
treated moist air as an ideal mixture of dry air and water vapour. It gave
the formulae and the chart of humidity ratio against temperature that
air-conditioning engineers still use.

*Implementation:* :func:`engineeringkit.thermodynamics.humidity_ratio`,
:func:`engineeringkit.thermodynamics.relative_humidity`,
:func:`engineeringkit.thermodynamics.dew_point`, and for the processes
:func:`engineeringkit.thermodynamics.moist_air_enthalpy`,
:func:`engineeringkit.thermodynamics.adiabatic_saturation_temperature`,
:func:`engineeringkit.thermodynamics.adiabatic_mixing`,
:func:`engineeringkit.thermodynamics.cooling_coil` and
:func:`engineeringkit.thermodynamics.cooling_tower`. The examples draw the
chart, follow air through a cooling coil, condition a room and size a cooling
tower.

*References:* W. H. Carrier, "Rational psychrometric formulae," *Trans.
ASME* 33 (1911) 1005-1053.

.. minigallery:: ../../examples/thermodynamics/psychrometrics/plot_01_carrier_psychrometrics.py ../../examples/thermodynamics/psychrometrics/plot_02_carrier_air_conditioning.py

1955 -- Pitzer's Acentric Factor
--------------------------------

Van der Waals's law of corresponding states puts every fluid on one chart of
compressibility against reduced temperature and pressure. That holds for the
noble gases but not for larger molecules. Kenneth Pitzer and his Berkeley
group measured how far a fluid's reduced vapour pressure at
:math:`T_r = 0.7` falls below the simple fluids' 0.1. They defined the
acentric factor :math:`\omega = -\log_{10}p_r^{sat} - 1` and wrote
:math:`Z = Z^0 + \omega Z^1`. Three parameters, :math:`T_c`, :math:`p_c` and
:math:`\omega`, now describe most fluids to engineering accuracy.

*Implementation:* :func:`engineeringkit.thermodynamics.pitzer_acentric_factor`,
:func:`engineeringkit.thermodynamics.pitzer_virial_compressibility` and
:func:`engineeringkit.thermodynamics.compressibility_factor`. The example
finds water's acentric factor from IF97 and the compressibility of n-butane.

*References:* K. S. Pitzer, "The volumetric and thermodynamic properties of
fluids. I. Theoretical basis and virial coefficients," *J. Am. Chem. Soc.*
77 (1955) 3427-3433, :doi:`10.1021/ja01618a001`; K. S. Pitzer, D. Z. Lippmann,
R. F. Curl, C. M. Huggins and D. E. Petersen, "... II. Compressibility
factor, vapor pressure and entropy of vaporization," *J. Am. Chem. Soc.* 77
(1955) 3433-3440, :doi:`10.1021/ja01618a002`.

.. minigallery:: ../../examples/thermodynamics/properties/plot_08_pitzer_acentric_factor.py

1956 -- Rant Names Exergy
-------------------------

Gibbs's available energy, Gouy's useful energy and Keenan's availability
were the same quantity under different names. Zoran Rant proposed one short
word, exergy, for the work that energy could deliver when brought into
equilibrium with the environment. For heat at :math:`T_s` it is
:math:`Q(1 - T_0/T_s)`. Exergy balances show where work potential is lost,
which energy balances cannot.

*Implementation:* :func:`engineeringkit.thermodynamics.exergy_heat`. The
example plots the quality of heat against temperature and draws up the
exergy balance of a gas boiler.

*References:* Z. Rant, "Exergie, ein neues Wort für 'technische
Arbeitsfähigkeit'," *Forsch. Ingenieurwes.* 22 (1956) 36-37.

.. minigallery:: ../../examples/thermodynamics/exergy/plot_02_rant_exergy.py

1975 -- Curzon and Ahlborn: Efficiency at Maximum Power
-------------------------------------------------------

A Carnot engine reaches its efficiency only by running infinitely slowly.
Frank Curzon and Boye Ahlborn let a reversible engine exchange heat with
its reservoirs through finite conductances and maximized the power. The
efficiency at maximum power is :math:`1 - \sqrt{T_c/T_h}`, whatever the
conductances. Novikov and Chambadal had found it in 1957. It matches real
power stations better than Carnot's bound and founded finite-time
thermodynamics.

*Implementation:*
:func:`engineeringkit.thermodynamics.curzon_ahlborn_efficiency`. The example
maximizes the power numerically and compares three power stations.

*References:* F. L. Curzon and B. Ahlborn, "Efficiency of a Carnot engine
at maximum power output," *Am. J. Phys.* 43 (1975) 22-24.

.. minigallery:: ../../examples/thermodynamics/cycles/plot_06_curzon_ahlborn.py

1976 -- Peng and Robinson's Equation of State
---------------------------------------------

Ding-Yu Peng and Donald Robinson, at the University of Alberta, needed one
equation for the vapour-liquid equilibria of natural-gas mixtures. They kept
van der Waals's repulsion :math:`RT/(v - b)` and changed the attraction to
:math:`a\alpha(T)/(v^2 + 2bv - b^2)`, with :math:`\alpha` tied to Pitzer's
acentric factor. The critical compressibility became 0.307 and liquid
densities improved, while the equation stayed a cubic that is solved in closed
form. It remains the default equation of state of process simulators. With an
ideal-gas heat capacity it also describes a refrigerant well enough to design
a cycle with.

*Implementation:* :func:`engineeringkit.thermodynamics.peng_robinson_pressure`,
:func:`engineeringkit.thermodynamics.peng_robinson_compressibility`,
:func:`engineeringkit.thermodynamics.peng_robinson_saturation_pressure`,
:func:`engineeringkit.thermodynamics.refrigerant_properties`,
:class:`engineeringkit.thermodynamics.CubicFluid` and
:data:`engineeringkit.thermodynamics.R134A`. The example draws the
isotherms and dome of R-134a and compares them with its tables.

*References:* D.-Y. Peng and D. B. Robinson, "A new two-constant equation of
state," *Ind. Eng. Chem. Fundam.* 15 (1976) 59-64, :doi:`10.1021/i160057a011`.

.. minigallery:: ../../examples/thermodynamics/properties/plot_09_peng_robinson.py

1997 -- IAPWS-IF97 Steam Properties
-----------------------------------

The International Association for the Properties of Water and Steam
replaced the 1967 industrial formulation with IAPWS-IF97, developed under
Wolfgang Wagner. Separate fast equations cover five regions. In the
compressed liquid and the vapour each is a dimensionless Gibbs free energy
from which every property follows by differentiation. The release gives
nine-digit verification values, and every power-plant calculation since uses
it.

*Implementation:* :func:`engineeringkit.thermodynamics.steam_properties`,
:func:`engineeringkit.thermodynamics.saturation_pressure`,
:func:`engineeringkit.thermodynamics.saturation_temperature` and
:class:`engineeringkit.thermodynamics.SteamState` (all five regions,
with array inputs). The examples check the verification values, draw the T-s
diagram and follow superheated steam into region 5 above 800 °C.

*References:* W. Wagner et al., "The IAPWS Industrial Formulation 1997 for
the thermodynamic properties of water and steam," *J. Eng. Gas Turbines
Power* 122 (2000) 150-182.

.. minigallery:: ../../examples/thermodynamics/properties/plot_06_iapws_if97.py ../../examples/thermodynamics/properties/plot_10_iapws_if97_region5.py

See Also
--------

- :doc:`/api/thermodynamics`
- :doc:`/guides/thermodynamics`
- :doc:`/history/heat_transfer_breakthroughs`
- :doc:`/history/fluid_engineering_breakthroughs`

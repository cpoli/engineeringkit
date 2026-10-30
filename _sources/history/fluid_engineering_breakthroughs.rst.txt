Breakthroughs in Fluid Engineering
==================================


.. include:: /_generated/nav/fluid_engineering.rst

.. epigraph::

   "When you can measure what you are speaking about, and express it in
   numbers, you know something about it."
   -- William Thomson (Lord Kelvin), lecture to the Institution of Civil
   Engineers, 1883

Engineering hydraulics is older than the calculus. Archimedes and Stevin
settled fluids at rest, Torricelli timed the jet from an orifice, and the
eighteenth century learned to measure moving water with Pitot's and
Venturi's instruments while Euler wrote the equation of every turbine.
Viscous flow in pipes then took a century to understand, from Hagen and Poiseuille's capillaries
through Reynolds's dye streak to the friction charts of Nikuradse,
Colebrook and Moody. Alongside them came the scaling laws of Froude and
Buckingham, Stokes's drag and Couette's viscometer, the weirs, jumps and
specific energy of open channels, Pelton's wheel, the pressure waves of
water hammer, cavitation in turbines and the hand methods for city water
networks. The entries trace the ideas behind
:mod:`engineeringkit.fluid_engineering`. Each links to the implementing
function and to a gallery example devoted to it.

.. contents:: Timeline
   :local:
   :depth: 1

c. 250 BC -- Archimedes's Principle
-----------------------------------

In the first book of *On Floating Bodies* Archimedes proved, from the
postulate that fluid parts pressed less are pushed by those pressed more,
that a solid lighter than a fluid sinks until it displaces its own weight,
and that a heavier solid is lightened by the weight of the fluid it
displaces. The buoyant force :math:`F_B = \rho g V` underlies every float,
hull and hydrometer, and the weighing of bodies in water to find their
density.

*Implementation:* :func:`engineeringkit.fluid_engineering.buoyancy_force`.
The example weighs King Hiero's crown in water, reads its silver content
back from the apparent weight and floats an iceberg.

*References:* Archimedes, *On Floating Bodies*, Book I, Props. 5-7, in
T. L. Heath (ed.), *The Works of Archimedes* (Cambridge: Cambridge
University Press, 1897).

.. minigallery:: ../../examples/fluid_engineering/hydrostatics/plot_01_archimedes_buoyancy.py

1586 -- Stevin's Hydrostatic Paradox
------------------------------------

Simon Stevin, the Flemish engineer of the Dutch dikes and sluices, showed
that the force of water on the bottom of a vessel equals the weight of a
column standing on the bottom, whatever the shape of the vessel and however
little water it holds. He went on to find the force on a vertical or
sloping wall by slicing it into strips, the first calculation of
engineering hydrostatics. In modern form, :math:`p = \rho g h`, and the
force on a plane is the pressure at its centroid times its area.

*Implementation:*
:func:`engineeringkit.fluid_engineering.hydrostatic_pressure`,
:func:`engineeringkit.fluid_engineering.pressure_head`,
:func:`engineeringkit.fluid_engineering.hydrostatic_force_plane` and
:func:`engineeringkit.fluid_engineering.center_of_pressure_depth`. The
example fills three vessels of different shape to one depth and finds the
thrust on a dam face and where it acts.

*References:* S. Stevin, *De Beghinselen des Waterwichts* (Leiden:
Plantijn, 1586).

.. minigallery:: ../../examples/fluid_engineering/hydrostatics/plot_02_stevin_hydrostatic_paradox.py

1644 -- Torricelli's Law of Efflux
----------------------------------

Evangelista Torricelli, who succeeded Galileo at Florence, stated at the
end of his treatise on the motion of heavy bodies that water spouting from
a small hole in a vessel has the speed a body would gain falling from the
free surface, so that a jet pointed upward would rise back almost to the
surface. In modern form :math:`v = \sqrt{2gh}`. The jet contracts just
past a sharp-edged hole, the vena contracta, and the discharge coefficient
of about 0.61 that measures this contraction became a staple of hydraulic
tables. Daniel Bernoulli later derived Torricelli's law from the
conservation of energy.

*Implementation:* :func:`engineeringkit.fluid_engineering.orifice_flow`
and :func:`engineeringkit.fluid_engineering.tank_drain_time`. The example
drains a tank through an orifice and checks the square-root law of its
level against a step-by-step integration.

*References:* E. Torricelli, "De motu aquarum," in *De motu gravium
naturaliter descendentium et proiectorum*, *Opera geometrica* (Florence:
Masse and de Landis, 1644) 191-203. The date and page range were checked
against secondary sources (Whewell; S. R. Bistafa's translation), not the
1644 printing.

.. minigallery:: ../../examples/fluid_engineering/measurement/plot_03_torricelli_efflux.py

1732 -- Pitot's Tube
--------------------

Henri Pitot faced the open end of a bent glass tube into the current of
the Seine and saw the water inside rise above the river surface by an
amount that grew with the speed. His readings showed that the river flows
fastest near the surface, against the belief that speed grows with depth.
The rise is the velocity head :math:`v^2/(2g)`, a relation Bernoulli's
*Hydrodynamica* explained six years later, and Pitot's tube, combined with
a static port, became the standard probe for water and air speed.

*Implementation:* :func:`engineeringkit.fluid_engineering.pitot_velocity`
and :func:`engineeringkit.fluid_engineering.dynamic_pressure`. The example
reads river speeds from the rise in the tube and finds the Mach number at
which the incompressible formula overestimates airspeed by 1%.

*References:* H. Pitot, "Description d'une machine pour mesurer la vitesse
des eaux courantes et le sillage des vaisseaux," *Hist. Acad. R. Sci.*
(1732) 363-376.

.. minigallery:: ../../examples/fluid_engineering/measurement/plot_01_pitot_tube.py

1754 -- Euler's Turbine Equation
--------------------------------

Leonhard Euler, at the Berlin academy, analysed the reaction water wheel
that Segner had built, and in a long memoir written in 1753 he worked out
the theory of turbines in general. Its central result is that the torque
on a rotor equals the change in the angular momentum of the water passing
through it, so the head a rotor gives or takes is
:math:`H = (u_2V_{t2} - u_1V_{t1})/g`. The equation still opens the
design of every pump, fan, compressor and turbine; Jakob Ackeret built and
tested a turbine to Euler's proportions in 1944.

*Implementation:*
:func:`engineeringkit.fluid_engineering.euler_turbomachine_head`, with
:func:`engineeringkit.fluid_engineering.affinity_laws`. The example draws
the ideal head of an impeller with backward, radial and forward blades and
recovers the head affinity law from the velocity triangle.

*References:* L. Euler, "Théorie plus complette des machines qui sont
mises en mouvement par la réaction de l'eau," *Mém. Acad. Sci. Berlin* 10
(for 1754, publ. 1756) 227-295 (Eneström E222); J. Ackeret,
"Untersuchung einer nach den Euler'schen Vorschlägen (1754) gebauten
Wasserturbine" (1944), English translation by S. R. Bistafa,
arXiv:2108.12048 (2021). The volume and pages are from the Euler Archive;
the memoir itself was not consulted.

.. minigallery:: ../../examples/fluid_engineering/pumps/plot_02_euler_turbine_equation.py

1797 -- Venturi's Tube
----------------------

Giovanni Battista Venturi, professor of physics at Modena, studied water
flowing through conical contractions and expansions. He found that the
pressure falls in the throat and is mostly recovered when the cone widens
gently. Ninety years later Clemens Herschel built the Venturi water meter
on the effect: with continuity and Bernoulli's equation, the pressure drop
from inlet to throat gives the volume flow, and the recovery makes the
meter cheap to run.

*Implementation:* :func:`engineeringkit.fluid_engineering.venturi_flow`,
:func:`engineeringkit.fluid_engineering.continuity_velocity` and
:func:`engineeringkit.fluid_engineering.bernoulli_pressure`. The example
draws the pressure along a meter and reads the flow back from it.

*References:* G. B. Venturi, *Recherches expérimentales sur le principe de
la communication latérale du mouvement dans les fluides* (Paris, 1797);
C. Herschel, "The Venturi water meter," *Trans. ASCE* 17 (1887) 228-258.

.. minigallery:: ../../examples/fluid_engineering/measurement/plot_02_venturi_meter.py

1828 -- Bélanger's Hydraulic Jump
---------------------------------

Jean-Baptiste Bélanger, an engineer of the Ponts et Chaussées, published
in 1828 an essay on the numerical solution of steady flow in open
channels. It gave the backwater equation for gradually varied flow with a
step method to integrate it, identified the critical flow condition where
that equation is singular, and treated the abrupt rise of a fast, shallow
stream to a deep, slow one: the hydraulic jump. In 1828 he applied
Bernoulli's energy balance across the jump, which is wrong, since the
roller dissipates energy. His later lectures (1838, printed 1841) used the
momentum principle and gave the relation still known as the Bélanger
equation, :math:`y_2/y_1 = (\sqrt{1 + 8\mathrm{Fr}_1^2} - 1)/2`.

*Implementation:* :func:`engineeringkit.fluid_engineering.sequent_depth`
and :func:`engineeringkit.fluid_engineering.hydraulic_jump_head_loss`.
The example finds the jump below a sluice gate, checks that momentum is
conserved and energy is not, and draws the depth ratio and loss against
the upstream Froude number.

*References:* J.-B. Bélanger, *Essai sur la solution numérique de
quelques problèmes relatifs au mouvement permanent des eaux courantes*
(Paris: Carilian-Goeury, 1828); H. Chanson, "Development of the Bélanger
equation and backwater equation by Jean-Baptiste Bélanger (1828)," *J.
Hydraul. Eng.* 135 (2009) 159-163,
:doi:`10.1061/(ASCE)0733-9429(2009)135:3(159)`. The account of the 1828 and
1838 derivations follows Chanson; the essay itself was not consulted.

.. minigallery:: ../../examples/fluid_engineering/open_channel/plot_02_belanger_hydraulic_jump.py

1839-1840 -- Hagen and Poiseuille: Laminar Pipe Flow
----------------------------------------------------

Gotthilf Hagen, a Prussian hydraulic engineer, timed water through brass
tubes, and Jean Léonard Marie Poiseuille, a physician studying the
circulation of the blood, timed liquids through glass capillaries. Both
found that slow flow through a tube is proportional to the pressure drop
and to the fourth power of the radius. Stokes and others later derived the
parabolic velocity profile behind it. In Darcy's form the same law is the
laminar friction factor :math:`f = 64/Re`.

*Implementation:* :func:`engineeringkit.fluid_engineering.poiseuille_flow`,
:func:`engineeringkit.fluid_engineering.poiseuille_velocity` and
:func:`engineeringkit.fluid_engineering.laminar_friction_factor`. The
example integrates the parabolic profile, checks the fourth-power law and
recovers the pressure drop from :math:`64/Re`.

*References:* G. Hagen, "Über die Bewegung des Wassers in engen
cylindrischen Röhren," *Ann. Phys. Chem.* 46 (1839) 423-442;
J. L. M. Poiseuille, "Recherches expérimentales sur le mouvement des
liquides dans les tubes de très petits diamètres," *C. R. Acad. Sci.* 11
(1840) 961-967, 1041-1048.

.. minigallery:: ../../examples/fluid_engineering/pipe_flow/plot_01_hagen_poiseuille.py

1845-1857 -- Weisbach and Darcy: Head Loss in Pipes
---------------------------------------------------

Julius Weisbach, at the Freiberg mining academy, wrote the friction loss in
a pipe as a coefficient times :math:`(L/D)\,v^2/(2g)` and tabulated loss
coefficients for bends, valves and contractions. Henry Darcy, engineer of
the Dijon water supply, measured losses in new and incrusted pipes up to
0.5 m in diameter and showed that the coefficient depends on the diameter
and the state of the wall. Their equation, with Darcy's name on the
friction factor, is how every pipeline is sized.

*Implementation:* :func:`engineeringkit.fluid_engineering.darcy_head_loss`,
:func:`engineeringkit.fluid_engineering.minor_head_loss` and
:func:`engineeringkit.fluid_engineering.velocity_head`, with the
:func:`~engineeringkit.fluid_engineering.visualizers.plot_grade_lines`
visualizer. The example follows the energy and hydraulic grade lines from a
reservoir through two pipes to a free jet.

*References:* J. Weisbach, *Lehrbuch der Ingenieur- und
Maschinen-Mechanik*, vol. 1 (Braunschweig: Vieweg, 1845); H. Darcy,
*Recherches expérimentales relatives au mouvement de l'eau dans les
tuyaux* (Paris: Mallet-Bachelier, 1857).

.. minigallery:: ../../examples/fluid_engineering/pipe_flow/plot_02_darcy_weisbach.py

1851 -- Stokes's Drag on a Sphere
---------------------------------

Pendulum clocks and the pendulums used to measure gravity were corrected
for the air they carried with them, and the corrections disagreed with
experiment. George Gabriel Stokes, Lucasian professor at Cambridge, solved
the slow viscous flow around a sphere and a cylinder in a paper read in
December 1850. For the sphere he found the drag :math:`3\pi\mu dV`, valid
while inertia is negligible, and from it the speed at which a fine
particle settles, which he used to explain why clouds float. At larger
Reynolds numbers the measured drag curve takes over, levelling off near
:math:`C_D = 0.4`.

*Implementation:* :func:`engineeringkit.fluid_engineering.stokes_drag`,
:func:`engineeringkit.fluid_engineering.stokes_settling_velocity`,
:func:`engineeringkit.fluid_engineering.sphere_drag_coefficient` and
:func:`engineeringkit.fluid_engineering.terminal_velocity`. The example
draws the sphere drag curve and finds how fast sand grains of four sizes
fall through water, with and without Stokes's law.

*References:* G. G. Stokes, "On the effect of the internal friction of
fluids on the motion of pendulums," *Trans. Camb. Phil. Soc.* 9, Part II
(1851) 8-106 (read 9 December 1850); F. M. White, *Viscous Fluid Flow*,
3rd ed. (New York: McGraw-Hill, 2006), for the drag-curve fit.

.. minigallery:: ../../examples/fluid_engineering/drag/plot_01_stokes_drag.py

1858 -- James Thomson's Triangular Notch
----------------------------------------

James Thomson, professor of engineering at Belfast and elder brother of
William Thomson (Lord Kelvin), proposed to the British Association in
1858 that water be measured over a right-angled V-notch cut in a weir
board. Because the wetted section stays the same shape as the head
changes, the discharge coefficient is nearly constant, and because the
flow grows as :math:`H^{5/2}`, small flows still give a readable head. He
reported experiments on 90° and wider notches in 1861. Rectangular
sharp-crested weirs, with coefficients refined by Francis (1855) and
Rehbock (1929), remain the companion device for larger flows.

*Implementation:*
:func:`engineeringkit.fluid_engineering.v_notch_weir_flow`,
:func:`engineeringkit.fluid_engineering.rectangular_weir_flow` and
:func:`engineeringkit.fluid_engineering.rehbock_discharge_coefficient`.
The example compares the rating curves of a V-notch and a rectangular weir
and the precision with which each must be read.

*References:* J. Thomson, "On experiments on the measurement of water by
triangular notches in weir boards," *Rep. Br. Assoc. Adv. Sci.* (1858),
Trans. Sect., 181; (1861), Trans. Sect., 151. Later sources disagree on
the 1858 page (133 or 181); the report itself was not consulted.

.. minigallery:: ../../examples/fluid_engineering/measurement/plot_04_thomson_v_notch_weir.py

1868 -- Froude's Ship Models
----------------------------

William Froude proposed to the Admiralty that a ship's resistance could be
found by towing a model at the "corresponding speed", the ship's speed
times the square root of the length ratio. At that speed model and ship
make geometrically similar waves: their Froude numbers
:math:`v/\sqrt{gL}` agree. Because skin friction scales differently,
Froude subtracted it, measured on towed planks, and scaled only the rest.
The Admiralty paid for his tank at Torquay, and full-scale towing of HMS
*Greyhound* in 1871 confirmed the method.

*Implementation:* :func:`engineeringkit.fluid_engineering.froude_number`
and :func:`engineeringkit.fluid_engineering.reynolds_number`. The example
scales a model test to a ship with Froude's split, using the ITTC 1957
friction line.

*References:* W. Froude, "Observations and suggestions on the subject of
determining by experiment the resistance of ships," memorandum to the
Admiralty (1868), in *The Papers of William Froude* (London: Institution
of Naval Architects, 1955).

.. minigallery:: ../../examples/fluid_engineering/similarity/plot_01_froude_ship_models.py

1880 -- Pelton's Water Wheel
----------------------------

The hydraulic mines of the California gold country drove their stamps and
pumps with hurdy-gurdy wheels: jets of water from high heads striking flat
paddles. Lester Allan Pelton, a millwright at Camptonville, built a wheel
whose buckets split the jet with a central ridge and turned each half
almost back on itself, and patented it in October 1880. The momentum
equation explains the gain: a jet turned through :math:`\theta` delivers
a fraction :math:`(1 - \cos\theta)/2` of its power at the best bucket
speed, half the jet speed, which is one half for flat paddles and nearly
one for Pelton's buckets. Pelton wheels still run the highest-head
hydroelectric plants.

*Implementation:* :func:`engineeringkit.fluid_engineering.jet_force` and
:func:`engineeringkit.fluid_engineering.pelton_wheel_power`, with
:func:`engineeringkit.fluid_engineering.bend_force` for the same momentum
balance on a pipe bend. The example compares the efficiency of flat
paddles and Pelton buckets against bucket speed.

*References:* L. A. Pelton, "Water-wheel," US Patent 233,692 (granted 26
October 1880; filed 3 July 1880).

.. minigallery:: ../../examples/fluid_engineering/momentum/plot_01_pelton_wheel.py

1883 -- Reynolds's Dye Experiment
---------------------------------

Osborne Reynolds let water flow from a tank into a glass tube with a bell
mouth and fed a thin streak of dye into the entrance. Below a critical
speed the streak stayed a straight thread; above it the streak broke up
into eddies that filled the tube. The critical speed changed with the
diameter and the temperature, but always at the same value of
:math:`\rho v D/\mu`. Reynolds also showed that the resistance changes
from proportional to the speed to nearly its square at the same point.

*Implementation:* :func:`engineeringkit.fluid_engineering.reynolds_number`.
The example computes the critical speed against tube diameter and water
temperature and shows the jump in friction factor across the transition.

*References:* O. Reynolds, "An experimental investigation of the
circumstances which determine whether the motion of water shall be direct
or sinuous, and of the law of resistance in parallel channels," *Phil.
Trans. R. Soc.* 174 (1883) 935-982.

.. minigallery:: ../../examples/fluid_engineering/pipe_flow/plot_03_reynolds_transition.py

1889 -- Manning's Open-Channel Formula
--------------------------------------

Robert Manning, chief engineer of the Irish Board of Works, tested the
many open-channel formulas of his day against the best measurements of
canals and rivers. He found that Chézy's coefficient grows roughly as the
sixth root of the hydraulic radius and proposed
:math:`V = R^{2/3}S^{1/2}/n` in metric units, with a roughness
coefficient :math:`n`. He read the paper to the Institution of Civil
Engineers of Ireland in December 1889; it appeared in print in 1891.
Tabulated values of :math:`n` made it the standard for channels, culverts
and sewers.

*Implementation:* :func:`engineeringkit.fluid_engineering.manning_velocity`,
:func:`engineeringkit.fluid_engineering.manning_flow` and
:func:`engineeringkit.fluid_engineering.hydraulic_radius`. The example
draws the rating curve of a concrete channel, finds the normal depth for a
design flow and checks that the flow is subcritical.

*References:* R. Manning, "On the flow of water in open channels and
pipes," *Trans. Inst. Civil Eng. Ireland* 20 (1891) 161-207.

.. minigallery:: ../../examples/fluid_engineering/open_channel/plot_01_manning_formula.py

1890 -- Couette's Rotating-Cylinder Viscometer
----------------------------------------------

Maurice Couette, in his doctoral work under Gabriel Lippmann at Paris,
measured the viscosity of water between two coaxial cylinders, the outer
one rotating and the inner one hung from a torsion wire. In a narrow gap
the flow is plain shear between a fixed and a sliding wall, with a
straight velocity profile and a uniform stress :math:`\mu U/h`, so the
torque on the inner cylinder gives the viscosity directly. He found the
torque proportional to the speed until the flow became unstable, and his
instrument is the ancestor of the rotational viscometers and rheometers
used today. Arnulph Mallock built a similar apparatus at about the same
time.

*Implementation:*
:func:`engineeringkit.fluid_engineering.couette_viscometer_torque`,
:func:`engineeringkit.fluid_engineering.couette_velocity`,
:func:`engineeringkit.fluid_engineering.couette_flow_rate` and
:func:`engineeringkit.fluid_engineering.couette_shear_stress`. The example
reads a viscosity back from a measured torque and adds a pressure
gradient to plane Couette flow.

*References:* M. Couette, "Études sur le frottement des liquides," *Ann.
Chim. Phys.*, 6th ser., 21 (1890) 433-510. The volume and pages are from
secondary sources; the paper itself was not consulted.

.. minigallery:: ../../examples/fluid_engineering/viscous/plot_01_couette_viscometer.py

1898 -- Joukowsky's Water Hammer
--------------------------------

Bursts in the new Moscow water mains led Nikolai Joukowsky to experiment
on pipes up to 760 m long. He showed that a sudden valve closure raises the
pressure by :math:`\rho a\,\Delta v`, where :math:`a` is the speed of the
pressure wave, slowed below the speed of sound in water by the
elasticity of the wall as Korteweg had calculated in 1878. The wave
travels to the reservoir and back in :math:`2L/a`, and any closure faster
than that produces the full rise. He reported the work in 1898, and the
memoir appeared in 1900.

*Implementation:*
:func:`engineeringkit.fluid_engineering.water_hammer_pressure`,
:func:`engineeringkit.fluid_engineering.pipe_wave_speed` and
:func:`engineeringkit.fluid_engineering.critical_closure_time`. The
example compares the surge in steel, iron and PVC pipe and draws the
pressure at the valve after a sudden closure.

*References:* N. Joukowsky, "Über den hydraulischen Stoss in
Wasserleitungsröhren," *Mém. Acad. Imp. Sci. St.-Pétersbourg*, 8th ser.,
9(5) (1900) 1-71; English translation by O. Simin, *Proc. AWWA* 24 (1904)
341-424; D. J. Korteweg, "Ueber die Fortpflanzungsgeschwindigkeit des
Schalles in elastischen Röhren," *Ann. Phys.* 241 (1878) 525-542.

.. minigallery:: ../../examples/fluid_engineering/transients/plot_01_joukowsky_water_hammer.py

1912 -- Bakhmeteff's Specific Energy
------------------------------------

Boris Bakhmeteff, then teaching hydraulics at St Petersburg, measured the
energy of an open-channel stream from its own bed rather than from a fixed
datum. Plotted against depth, this specific energy
:math:`E = y + q^2/(2gy^2)` has a least value at the critical depth, where
the Froude number is 1, and any larger energy can be carried at two
alternate depths, one tranquil and one rapid. Bélanger had found the
critical condition in 1828 as a singularity of the backwater equation;
Bakhmeteff tied it to minimum energy, and the specific-energy diagram of
his 1932 *Hydraulics of Open Channels* became the standard way to reason
about gates, steps and contractions. Salas Edwards and Domínguez reached
the same idea independently at about the same time.

*Implementation:* :func:`engineeringkit.fluid_engineering.specific_energy`,
:func:`engineeringkit.fluid_engineering.critical_depth` and
:func:`engineeringkit.fluid_engineering.alternate_depths`. The example
draws the specific-energy diagram and finds the water surface over a
raised step in the bed.

*References:* B. A. Bakhmeteff, *O neravnomernom dvizhenii zhidkosti v
otkrytom rusle* [Varied flow of liquids in open channels] (St Petersburg,
1912); B. A. Bakhmeteff, *Hydraulics of Open Channels* (New York:
McGraw-Hill, 1932); O. Castro-Orgaz, H. Chanson et al., "Boris A.
Bakhmeteff and the development of specific energy and momentum concepts,"
*J. Hydraul. Eng.* (2018). The 1912 book was not consulted; its priority
follows Castro-Orgaz and Chanson.

.. minigallery:: ../../examples/fluid_engineering/open_channel/plot_03_bakhmeteff_specific_energy.py

1914 -- Buckingham's Pi Theorem
-------------------------------

Edgar Buckingham, a physicist at the US Bureau of Standards, proved that
any physically meaningful relation among :math:`n` quantities whose
dimensions involve :math:`r` independent base units can be written as a
relation among :math:`n - r` dimensionless products. Rayleigh and Reynolds
had used dimensional reasoning case by case; Buckingham made it a method,
and his paper named the groups :math:`\Pi`. The exponents of the groups
form the null space of the dimension matrix.

*Implementation:* :func:`engineeringkit.fluid_engineering.buckingham_pi`
(exact rational row reduction in
:func:`engineeringkit.fluid_engineering.utils.integer_null_space`),
:func:`engineeringkit.fluid_engineering.euler_number` and
:func:`engineeringkit.fluid_engineering.weber_number`. The example derives
the four groups of pipe friction and the groups of droplet break-up, and
shows three dissimilar pipes with equal groups giving one Euler number.

*References:* E. Buckingham, "On physically similar systems; illustrations
of the use of dimensional equations," *Phys. Rev.* 4 (1914) 345-376.

.. minigallery:: ../../examples/fluid_engineering/similarity/plot_02_buckingham_pi.py

1924 -- Thoma's Cavitation Parameter
------------------------------------

Faster, lower-head water turbines of the early twentieth century had their
runners pitted by cavitation: vapour bubbles that form where the pressure
drops to the vapour pressure and collapse against the blades. Dieter
Thoma, at the Technische Hochschule München, showed that turbines of one
design begin to cavitate at the same ratio :math:`\sigma` of the suction
head margin to the machine head. The Thoma number, and the net positive
suction head behind it, set the setting level of every turbine and the
suction lift of every pump.

*Implementation:* :func:`engineeringkit.fluid_engineering.thoma_sigma`,
:func:`engineeringkit.fluid_engineering.npsh_available` and
:func:`engineeringkit.fluid_engineering.cavitation_number`. The example
finds the largest suction lift of a pump at three water temperatures.

*References:* D. Thoma, "Die Kavitation bei Wasserturbinen,"
*Wasserkraft-Jahrbuch* (Munich, 1924).

.. minigallery:: ../../examples/fluid_engineering/pumps/plot_01_thoma_cavitation.py

1933 -- Nikuradse's Sand-Roughened Pipes
----------------------------------------

Johann Nikuradse, in Ludwig Prandtl's institute at Göttingen, coated brass
pipes with sieved sand of six grain sizes and measured their friction over
a wide range of Reynolds numbers. After his smooth-pipe measurements of
1932, which fixed Prandtl and von Kármán's logarithmic law, he showed that
rough pipes follow the smooth law until the grains reach through the
viscous sublayer and then settle on a friction factor that depends on the
relative roughness alone.

*Implementation:*
:func:`engineeringkit.fluid_engineering.nikuradse_friction_factor` and
:func:`engineeringkit.fluid_engineering.smooth_pipe_friction_factor`. The
example draws the smooth and fully rough limits for Nikuradse's six
roughnesses and finds where each becomes fully rough.

*References:* J. Nikuradse, "Strömungsgesetze in rauhen Rohren,"
*VDI-Forschungsheft* 361 (1933); J. Nikuradse, "Gesetzmäßigkeiten der
turbulenten Strömung in glatten Rohren," *VDI-Forschungsheft* 356 (1932).

.. minigallery:: ../../examples/fluid_engineering/pipe_flow/plot_04_nikuradse_roughness.py

1936 -- Hardy Cross's Pipe-Network Method
-----------------------------------------

Hardy Cross, professor of structural engineering at the University of
Illinois, had already relaxed the moments in frames one joint at a time.
He applied the same idea to water distribution systems. Starting from
flows that satisfy continuity at every junction, he corrected one loop at
a time by a Newton step that makes the loop's head losses sum to zero.
The corrections preserve continuity, so the method converges on both sets
of conditions with only hand arithmetic. Network analysis in water
utilities relied on it for thirty years.

*Implementation:* :func:`engineeringkit.fluid_engineering.hardy_cross`
(returning a :class:`~engineeringkit.fluid_engineering.NetworkResult`) and
:func:`engineeringkit.fluid_engineering.pipe_resistance`. The example
balances a two-loop network and checks it against a direct Newton solve.

*References:* H. Cross, "Analysis of flow in networks of conduits or
conductors," *Univ. Illinois Eng. Exp. Stn. Bull.* 286 (1936).

.. minigallery:: ../../examples/fluid_engineering/networks/plot_01_hardy_cross.py

1939 -- Colebrook's Turbulent Friction Factor
---------------------------------------------

Commercial pipes do not have Nikuradse's uniform sand, and their friction
leaves the smooth law gradually instead of with his dip. Cyril Colebrook,
after experiments with Cedric White on pipes coated with mixed grains,
added the smooth and rough terms inside one logarithm. The result fits
commercial pipes from smooth to fully rough with an "equivalent sand
roughness" for each material. It is implicit in the friction factor, so
it is solved by iteration.

*Implementation:*
:func:`engineeringkit.fluid_engineering.colebrook_friction_factor`. The
example solves the equation by fixed-point substitution and shows how it
joins the smooth and fully rough limits.

*References:* C. F. Colebrook, "Turbulent flow in pipes, with particular
reference to the transition region between the smooth and rough pipe
laws," *J. Inst. Civil Eng.* 11 (1939) 133-156; C. F. Colebrook and
C. M. White, "Experiments with fluid friction in roughened pipes," *Proc.
R. Soc. A* 161 (1937) 367-381.

.. minigallery:: ../../examples/fluid_engineering/pipe_flow/plot_05_colebrook_equation.py

1944 -- Moody's Diagram
-----------------------

Lewis Moody, professor of hydraulics at Princeton, plotted Colebrook's
equation for a range of relative roughnesses on one log-log chart with the
laminar line, and tabulated equivalent roughnesses for commercial
materials. He showed how to use the chart for the three classic pipe
problems: finding the loss, the flow or the diameter. Hunter Rouse had
drawn a similar chart in 1943, but Moody's became the one in every
textbook.

*Implementation:*
:func:`~engineeringkit.fluid_engineering.visualizers.plot_moody_diagram`,
with :func:`engineeringkit.fluid_engineering.colebrook_friction_factor`
and :func:`engineeringkit.fluid_engineering.laminar_friction_factor`. The
example draws the chart, reads one point from it, and solves Moody's
diameter problem.

*References:* L. F. Moody, "Friction factors for pipe flow," *Trans. ASME*
66 (1944) 671-684.

.. minigallery:: ../../examples/fluid_engineering/pipe_flow/plot_06_moody_diagram.py

1976 -- Swamee and Jain's Explicit Friction Factor
--------------------------------------------------

Iterating Colebrook's equation was tedious on the programmable calculators
of the 1970s. Prabhata Swamee and Akalank Jain fitted explicit formulas for
the friction factor, the flow and the diameter of a pipe, each within a few
percent of the implicit solution over the range of engineering pipes.
Their friction factor became the standard explicit substitute for
Colebrook in hand calculations and spreadsheets.

*Implementation:*
:func:`engineeringkit.fluid_engineering.swamee_jain_friction_factor`. The
example maps its error against Colebrook and checks their explicit
diameter formula against an iterative solution.

*References:* P. K. Swamee and A. K. Jain, "Explicit equations for
pipe-flow problems," *J. Hydraul. Div. ASCE* 102 (1976) 657-664.

.. minigallery:: ../../examples/fluid_engineering/pipe_flow/plot_07_swamee_jain.py

See Also
--------

- :doc:`/api/fluid_engineering`
- :doc:`/guides/fluid_engineering`
- :doc:`/history/aerodynamics_breakthroughs`
- :doc:`/history/heat_transfer_breakthroughs`

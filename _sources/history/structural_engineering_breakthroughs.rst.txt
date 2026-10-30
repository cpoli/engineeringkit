Breakthroughs in Structural Engineering
=======================================


.. include:: /_generated/nav/structural_engineering.rst

.. epigraph::

   "Strength is essential but otherwise unimportant."
   -- attributed to Hardy Cross

Builders raised cathedrals and bridges for centuries by rule and
precedent. Structural engineering became a science in three steps.
Seventeenth- and eighteenth-century natural philosophers asked how strong
a beam is and when a column buckles. Engineers of the nineteenth century,
faced with iron railway bridges, worked out the theory of beams, trusses
and statically indeterminate structures, and the energy theorems that tie
loads to displacements. The twentieth century added plastic collapse,
relaxation methods for frames and, finally, matrix methods run on a
computer. The entries trace the ideas behind
:mod:`engineeringkit.structural_engineering`. Each links to the
implementing function and to a gallery example devoted to it.

.. contents:: Timeline
   :local:
   :depth: 1

1638 -- Galileo's Cantilever
----------------------------

In the second day of the *Two New Sciences*, Galileo Galilei asked what
load breaks a beam built into a wall. He treated the beam as a lever
turning about the bottom edge of the wall section, with the whole section
pulling with the same stress. The breaking moment then grows with the
width and the *square* of the depth, so a plank set on edge is stronger
than the same plank laid flat, and large animals and buildings need
proportionally thicker members. The scaling was right. The constant was
not: it took Navier's neutral axis at the centroid to divide it by three.

*Implementation:* :func:`engineeringkit.structural_engineering.rectangular_section_modulus`
and :func:`engineeringkit.structural_engineering.rectangular_second_moment`.
The example compares Galileo's breaking moment with the elastic one and
checks the :math:`bh^2` scaling.

*References:* G. Galilei, *Discorsi e dimostrazioni matematiche intorno a
due nuove scienze* (Elsevier, Leiden, 1638), second day.

.. minigallery:: ../../examples/structural_engineering/beams/plot_01_galileo_cantilever.py

1675 -- Hooke's Inverted Chain
------------------------------

At the end of *A Description of Helioscopes* Robert Hooke printed a Latin
anagram for "the true mathematical and mechanical form of all manner of
arches". Decoded after his death, it read *ut pendet continuum flexile,
sic stabit contiguum rigidum inversum*: as hangs the flexible line, so but
inverted will stand the rigid arch. A hanging chain can only pull, so it
takes the shape in which its loads cause no bending. Turned upside down,
the same shape is an arch that only pushes. A three-hinged arch is
statically determinate, and its moment is the simple-beam moment less the
thrust times the rise, :math:`M = M_0 - Hy`, which vanishes on the
funicular shape.

*Implementation:* :func:`engineeringkit.structural_engineering.three_hinged_arch`
and :class:`engineeringkit.structural_engineering.ArchResult`. The example
checks that a parabolic arch under a load uniform along the span carries
no moment and the same force as the matching cable, and shows the bending
a point load causes.

*References:* R. Hooke, *A Description of Helioscopes, and Some Other
Instruments* (John Martyn, London, 1676 [1675]), postscript; the solution
was published in R. Waller (ed.), *The Posthumous Works of Robert Hooke*
(1705). Dates checked against secondary accounts only.

.. minigallery:: ../../examples/structural_engineering/cables/plot_01_hooke_inverted_chain.py


1678 -- Hooke's *Ut tensio, sic vis*
------------------------------------

Robert Hooke published his law as an anagram in 1676 and its solution two
years later: "as the extension, so the force". Springs, wires and wooden
rods all stretched in proportion to the load. Written for a bar of length
:math:`L`, area :math:`A` and modulus :math:`E` (Young's, 1807), the
extension is :math:`\delta = PL/(AE)`, and the energy stored is the area
under the load-extension line, :math:`U = P\delta/2`. The linear law is
the constitutive relation that turns every later structural theory from
statics into a calculation of displacements.

*Implementation:* :func:`engineeringkit.structural_engineering.axial_deformation`
and :func:`engineeringkit.structural_engineering.strain_energy_axial`. The
example repeats Hooke's experiment on a steel wire and checks the
extension and the stored energy.

*References:* R. Hooke, *Lectures de Potentia Restitutiva, or of Spring*
(John Martyn, London, 1678).

.. minigallery:: ../../examples/structural_engineering/members/plot_01_hooke_axial.py

1691 -- The Catenary
--------------------

Galileo had guessed that a hanging chain is a parabola. In 1690 Jacob
Bernoulli asked the readers of *Acta Eruditorum* for its true shape, and
in 1691 three answers appeared there, by his brother Johann, by Leibniz and
by Huygens: the *catenary* :math:`y = a(\cosh(x/a) - 1)`, where
:math:`a = H/w` is the horizontal tension divided by the weight per unit
length. The tension at a support is :math:`T = H + wf`. A cable that
carries a deck uniform along the span, rather than its own weight, hangs
in a parabola with :math:`H = wL^2/(8f)`, and for shallow cables the two
curves are nearly the same.

*Implementation:* :func:`engineeringkit.structural_engineering.catenary_cable`
and :func:`engineeringkit.structural_engineering.parabolic_cable`. The
example sizes a power-line conductor both ways and shows where the
parabola stops being a good approximation.

*References:* Joh. Bernoulli, G. W. Leibniz and C. Huygens, solutions of
the catenary problem, *Acta Eruditorum* (1691). Checked against secondary
accounts only.

.. minigallery:: ../../examples/structural_engineering/cables/plot_02_catenary.py


1744 -- Euler's Column Buckling
-------------------------------

In an appendix on elastic curves to his book on the calculus of
variations, Leonhard Euler found the load at which a straight elastic
column can stay in equilibrium bent as well as straight. For a pin-ended
column it is :math:`P_{cr} = \pi^2EI/L^2`. Other end conditions change the
length of the buckled half-wave, written today as an effective length
:math:`KL`. Euler's load does not depend on the strength of the material
at all, which is why slender columns fail by stability and stocky ones by
crushing. Engineers distrusted the formula for a century, because the
columns they tested were not slender enough for it to apply.

*Implementation:* :func:`engineeringkit.structural_engineering.euler_buckling_load`,
:func:`engineeringkit.structural_engineering.radius_of_gyration` and
:func:`engineeringkit.structural_engineering.slenderness_ratio`. The example
finds Euler's load as the lowest eigenvalue of the finite-difference column
equation and plots the Euler hyperbola against yield.

*References:* L. Euler, *Methodus inveniendi lineas curvas maximi minimive
proprietate gaudentes* (Bousquet, Lausanne, 1744), Additamentum I, "De
curvis elasticis".

.. minigallery:: ../../examples/structural_engineering/stability/plot_01_euler_buckling.py

1784 -- Coulomb's Torsion
-------------------------

Charles-Augustin Coulomb hung a weight from a thin wire, twisted it and
timed the oscillations. The restoring torque was proportional to the
angle of twist, to the fourth power of the diameter and inversely to the
length, the law :math:`\phi = TL/(GJ)` with :math:`J = \pi d^4/32`. The
torsion balance built on it measured the electric force the next year and,
in Cavendish's hands, the density of the Earth. For engineers it gave the
shear stress :math:`\tau = Tr/J` and stiffness of shafts.

*Implementation:* :func:`engineeringkit.structural_engineering.torsion_angle`
and :func:`engineeringkit.structural_engineering.torsion_shear_stress`.
The example simulates Coulomb's torsion pendulum, recovers the shear
modulus from its period and checks the fourth-power law.

*References:* C.-A. Coulomb, "Recherches théoriques et expérimentales sur
la force de torsion et sur l'élasticité des fils de métal," *Hist. Acad.
Roy. Sci.* (1784, published 1787) 229-269.

.. minigallery:: ../../examples/structural_engineering/members/plot_02_coulomb_torsion.py

1807 -- Young's Eccentric Column
--------------------------------

In his *Lectures on Natural Philosophy* Thomas Young analysed a column
loaded off its axis. The load bends it, the bending increases the lever arm
of the load, and the peak stress follows the *secant formula*,
:math:`\sigma_\max = (P/A)[1 + (ec/r^2)\sec((KL/2r)\sqrt{P/EA})]`. Stress
grows faster than load and becomes unbounded at the Euler load, so a
column's factor of safety must be applied to the load. Solving the
formula for the load at first yield gives a family of column curves, one
for each eccentricity ratio :math:`ec/r^2`.

*Implementation:* :func:`engineeringkit.structural_engineering.secant_formula_stress`
and :func:`engineeringkit.structural_engineering.secant_formula_load`. The
example shows the nonlinear rise of stress with load and draws the
secant-formula column curves.

*References:* T. Young, *A Course of Lectures on Natural Philosophy and the
Mechanical Arts* (Joseph Johnson, London, 1807), vol. 2. The attribution
follows the secondary literature; the lectures were not consulted.

.. minigallery:: ../../examples/structural_engineering/stability/plot_03_young_secant_formula.py


1826 -- Navier's Beam Theory
----------------------------

Claude-Louis Navier's lectures at the École des ponts et chaussées put
beam theory in the form still taught. Plane sections remain plane, the
neutral axis passes through the centroid, the stress varies linearly as
:math:`\sigma = My/I`, and the curvature of the beam is :math:`M/(EI)`.
With these, the deflection of any statically determinate beam follows
from two integrations of the moment diagram. The *Résumé des leçons* also
treated arches and suspension bridges, and made elastic theory a working
tool of French engineers.

*Implementation:* :func:`engineeringkit.structural_engineering.bending_stress`
and :class:`engineeringkit.structural_engineering.SimplySupportedBeam`. The
example checks the midspan deflection :math:`5wL^4/(384EI)` and that the
stress profile carries the applied moment with no net force.

*References:* C.-L. Navier, *Résumé des leçons données à l'École des ponts
et chaussées sur l'application de la mécanique à l'établissement des
constructions et des machines* (Paris, 1826).

.. minigallery:: ../../examples/structural_engineering/beams/plot_02_navier_beam.py

1829 -- Poisson's Ratio
-----------------------

A bar stretched along its axis gets thinner. Siméon Denis Poisson's
molecular theory of elasticity predicted that the lateral strain is
always one quarter of the axial strain. Experiments by Wertheim and others
soon showed that the ratio :math:`\nu` depends on the material, and the
argument over it led to the modern two-constant theory of isotropic
elasticity. :math:`\nu` lies between -1 and 1/2 for a stable material,
and :math:`\nu = 1/2` means no change of volume.

*Implementation:* :func:`engineeringkit.structural_engineering.poisson_lateral_strain`.
The example compares the contraction of steel, rubber and cork and checks
the first-order volume change :math:`(1 - 2\nu)\varepsilon`.

*References:* S.-D. Poisson, "Mémoire sur l'équilibre et le mouvement des
corps élastiques," *Mém. Acad. Roy. Sci. Inst. France* 8 (1829) 357-570.

.. minigallery:: ../../examples/structural_engineering/members/plot_03_poisson_ratio.py

1833 -- Lamé's Thick-Walled Cylinder
------------------------------------

Gabriel Lamé and Émile Clapeyron, engineers in the Russian service, sent the
Paris Academy a memoir on the internal equilibrium of elastic solids,
published in 1833. Among its applications is the hollow cylinder under
internal and external pressure: :math:`\sigma_r = A - B/r^2` and
:math:`\sigma_\theta = A + B/r^2`. The hoop stress is greatest at the bore.
Adding thickness helps less and less, which drove the shrink-fitted and
wire-wound gun barrels of the nineteenth century. For a thin wall the
solution reduces to the boiler formula :math:`pr/t`.

*Implementation:* :func:`engineeringkit.structural_engineering.thick_wall_cylinder_stresses`
and :func:`engineeringkit.structural_engineering.thin_wall_vessel_stresses`.
The example plots the stresses through a thick wall and finds the radius
ratio, about 10, at which the thin-wall formula is 5% low.

*References:* G. Lamé and É. Clapeyron, "Mémoire sur l'équilibre intérieur
des corps solides homogènes," *Mémoires présentés par divers savans à
l'Académie royale des sciences de l'Institut de France* 4 (1833) 465-562.
The volume and page are from secondary bibliographies.

.. minigallery:: ../../examples/structural_engineering/vessels/plot_01_lame_thick_cylinder.py


1847 -- Whipple's Method of Joints
----------------------------------

American railway bridges of the 1840s were built from rules of thumb, and
some of them fell down. Squire Whipple, a canal engineer from Utica, showed
that the force in every member of a pin-jointed truss follows from statics
alone. Each joint is in equilibrium, so its two equations give the forces
in two members. Work from a support across the truss, two unknowns at a
time, and every force is found. His book was the first to analyse trusses
rather than test them, and Herman Haupt reached the same results
independently in 1851.

*Implementation:* :func:`engineeringkit.structural_engineering.method_of_joints`.
The example solves a six-panel Pratt truss and checks every chord force
against the method of sections, :math:`M/h`.

*References:* S. Whipple, *A Work on Bridge Building* (H. H. Curtiss,
Utica, 1847).

.. minigallery:: ../../examples/structural_engineering/trusses/plot_01_whipple_method_of_joints.py

1856 -- Jourawski's Shear Stress
--------------------------------

Building timber bridges for the St Petersburg-Moscow railway, Dmitrii
Ivanovich Jourawski saw large beams split along their length at
mid-depth, where the bending stress is zero. He balanced the change in
bending force on the part of the section beyond a cut against the shear
on the cut. That gave the shear flow :math:`q = VQ/I`, which the keys of a
built-up beam must carry, and the shear stress :math:`\tau = VQ/(It)`. The
analysis was done in 1855 and published in Paris the next year.
Saint-Venant later found the exact solution, but Jourawski's formula is
the one engineers use.

*Implementation:* :func:`engineeringkit.structural_engineering.transverse_shear_stress`,
:func:`engineeringkit.structural_engineering.shear_flow`,
:func:`engineeringkit.structural_engineering.composite_section` and
:meth:`engineeringkit.structural_engineering.SectionProperties.first_moment`.
The example checks that the shear stress in a rectangle and an I-beam adds
up to the shear force, and spaces the screws of a flitch beam from its
transformed section.

*References:* D. J. Jourawski, "Sur la résistance d'un corps prismatique et
d'une pièce composée en bois ou en tôle de fer à une force perpendiculaire
à leur longueur," *Annales des ponts et chaussées*, 3e série, 12 (1856)
328-351.

.. minigallery:: ../../examples/structural_engineering/shear/plot_01_jourawski_shear_stress.py


1857 -- Clapeyron's Three-Moment Equation
-----------------------------------------

A beam continuous over several supports is statically indeterminate: the
support moments cannot be found from equilibrium alone. Émile Clapeyron
met the problem while rebuilding the Asnières railway bridge over the
Seine. He required the slope of the beam to be continuous over each
interior support, which gives one linear equation per support linking its
moment to those of its two neighbours. The system is tridiagonal and
small, and it made continuous girders, lighter than a row of simple spans,
practical to design.

*Implementation:* :func:`engineeringkit.structural_engineering.three_moment_equation`.
The example checks the textbook coefficients for equal spans, then
integrates the curvature of an unequal three-span beam to show it does not
move at the interior supports.

*References:* É. Clapeyron, "Calcul d'une poutre élastique reposant
librement sur des appuis inégalement espacés," *C. R. Acad. Sci. Paris* 45
(1857) 1076-1080.

.. minigallery:: ../../examples/structural_engineering/beams/plot_03_clapeyron_three_moment.py

1864 -- Maxwell's Rigidity Count and Reciprocal Theorem
-------------------------------------------------------

James Clerk Maxwell's short paper on frames made two contributions. A
plane frame of :math:`j` joints needs :math:`2j - 3` bars to be stiff, or,
counting the :math:`r` support reactions, :math:`m + r = 2j` for a
statically determinate truss. The count is necessary but not sufficient,
since bars can be put where they are not needed. The paper also stated the
reciprocal theorem: in a linear elastic structure, the deflection at
:math:`i` caused by a unit load at :math:`j` equals the deflection at
:math:`j` caused by a unit load at :math:`i`. Betti generalised it in 1872,
and Mohr built the force method for indeterminate trusses on it.

*Implementation:* :func:`engineeringkit.structural_engineering.truss_determinacy`
and :func:`engineeringkit.structural_engineering.maxwell_betti_reciprocity`.
The example counts four trusses, one of which passes the count and is
still a mechanism, and checks that the flexibility matrix of a redundant
truss is symmetric.

*References:* J. C. Maxwell, "On the calculation of the equilibrium and
stiffness of frames," *Phil. Mag.* 27 (1864) 294-299; E. Betti, "Teoria
dell'elasticità," *Il Nuovo Cimento* 7-8 (1872) 5-21.

.. minigallery:: ../../examples/structural_engineering/trusses/plot_02_maxwell_rigidity.py

1867 -- Winkler's Influence Lines
---------------------------------

A railway bridge must be checked for a train standing in every position.
Emil Winkler, and Otto Mohr independently the next year, plotted the
value of one response, such as the moment at a chosen section, against
the position of a single unit load. With this *influence line* the
response to any train is the sum of the axle loads times the ordinates
under them, and the worst position of the train can be read off the plot.
Influence lines remain the standard tool for moving loads on bridges and
crane girders.

*Implementation:* :func:`engineeringkit.structural_engineering.influence_line`.
The example finds the worst position of a four-axle locomotive on a 20 m
span and checks it against superposed point-load solutions.

*References:* E. Winkler, *Die Lehre von der Elasticität und Festigkeit*
(Dominicus, Prague, 1867); O. Mohr, "Beitrag zur Theorie der Holz- und
Eisen-Constructionen," *Z. Architekten- und Ingenieur-Vereins zu
Hannover* 14 (1868) 19-51.

.. minigallery:: ../../examples/structural_engineering/beams/plot_04_winkler_influence_lines.py

1874 -- Mohr's Unit-Load Method
-------------------------------

Otto Mohr applied the principle of virtual work to trusses. To find how
far one joint moves, put an imaginary unit load there. Its member forces
:math:`n`, working through the real extensions :math:`NL/(AE)`, equal the
unit load times the displacement: :math:`\delta = \sum nNL/(AE)`. Both
force systems come from statics, so a determinate truss needs no
stiffness matrix. The same argument gives :math:`\int mM/(EI)\,dx` for
beams, and, with the redundant forces as unknowns, the force method for
indeterminate structures that Maxwell had begun ten years earlier.

*Implementation:* :func:`engineeringkit.structural_engineering.unit_load_truss_deflection`
and :func:`engineeringkit.structural_engineering.unit_load_beam_deflection`.
The example finds the midspan deflection of a Warren truss member by
member and checks it against the direct stiffness method.

*References:* O. Mohr, "Beitrag zur Theorie des Fachwerks," *Zeitschrift
des Architekten- und Ingenieur-Vereins zu Hannover* 20 (1874) 509-526.

.. minigallery:: ../../examples/structural_engineering/energy/plot_02_mohr_unit_load.py


1879 -- Castigliano's Theorem
-----------------------------

In his 1873 engineering thesis at Turin, extended into a book in 1879,
Alberto Castigliano showed that the partial derivative of the strain
energy with respect to a load is the displacement at that load, in its
direction. For a beam this is the integral of
:math:`(M/EI)\,\partial M/\partial P` along the span. A *dummy* load,
set to zero after differentiating, gives the displacement where no load
acts. The companion theorem of least work found the redundant forces of an
indeterminate structure by minimising the energy.

*Implementation:* :func:`engineeringkit.structural_engineering.castigliano_deflection`
(integrating :math:`M\,\partial M/\partial P/(EI)` with Simpson's rule) and
:func:`engineeringkit.structural_engineering.strain_energy_bending_constant_moment`.
The example finds the tip deflection of a tapered cantilever, checks it by
double integration and by the energy balance, and uses a dummy load at
midspan.

*References:* A. Castigliano, *Théorie de l'équilibre des systèmes
élastiques et ses applications* (A. F. Negro, Turin, 1879).

.. minigallery:: ../../examples/structural_engineering/energy/plot_01_castigliano.py

1882 -- Mohr's Circle
---------------------

Otto Mohr showed that the normal and shear stresses on every plane
through a point, plotted against each other, fall on one circle. Its
centre is the mean normal stress, its radius the largest in-plane shear
stress, and its intersections with the axis are the principal stresses.
Turning the element by :math:`\theta` moves the point by :math:`2\theta`
around the circle. Culmann had drawn the stresses in a beam graphically in
1866; Mohr's construction made the transformation one picture, and he
used it to state his theory of failure. For a three-dimensional state he
drew three circles, one through each pair of principal stresses; the
largest gives the absolute maximum shear :math:`(\sigma_1 - \sigma_3)/2`.

*Implementation:* :func:`engineeringkit.structural_engineering.mohr_circle`,
:func:`engineeringkit.structural_engineering.stress_transformation` and
:func:`engineeringkit.structural_engineering.principal_stresses_3d`.
The example traces the circle by rotating an element, checks the principal
stresses against the eigenvalues of the stress tensor, finds the plane of
maximum shear, and draws the three circles for a tank wall, where the
largest shear acts out of the plane.

*References:* O. Mohr, "Über die Darstellung des Spannungszustandes und
des Deformationszustandes eines Körperelementes," *Zivilingenieur* 28
(1882) 113-156.

.. minigallery:: ../../examples/structural_engineering/stress/plot_01_mohr_circle.py

1889 -- Engesser's Tangent-Modulus Load
---------------------------------------

Columns of intermediate length failed well below Euler's load. Friedrich
Engesser argued that once the stress passes the proportional limit, a
column bends against the material's tangent modulus
:math:`E_t = d\sigma/d\varepsilon`, so the critical stress solves
:math:`\sigma = \pi^2E_t(\sigma)/(KL/r)^2`. Engesser later replaced it by
a "double-modulus" load, which was thought correct for half a century,
until F. R. Shanley showed in 1947 that real columns start to bend at the
tangent-modulus load. Column curves in design codes descend from it.

*Implementation:* :func:`engineeringkit.structural_engineering.tangent_modulus_load`
(solving for the buckling stress with Brent's method). The example builds
the column curve of a Ramberg-Osgood aluminium alloy and checks each
stress against Engesser's equation.

*References:* F. Engesser, "Ueber die Knickfestigkeit gerader Stäbe,"
*Z. Architekten- und Ingenieur-Vereins zu Hannover* 35 (1889) 455-462;
F. R. Shanley, "Inelastic column theory," *J. Aeronaut. Sci.* 14 (1947)
261-268.

.. minigallery:: ../../examples/structural_engineering/stability/plot_02_engesser_tangent_modulus.py

1893 -- Johnson's Parabola
--------------------------

Tests on real columns fell well below Euler's hyperbola at low
slenderness, where columns yield or buckle inelastically. J. B. Johnson,
C. W. Bryan and F. E. Turneaure fitted a parabola that starts at the yield
stress and meets Euler's curve tangentially at half of it, at the
transition slenderness :math:`C_c = \sqrt{2\pi^2E/\sigma_y}`. End conditions
enter through the effective length :math:`KL`. With a factor of safety,
the same curve was the American steel column formula for most of the
twentieth century.

*Implementation:* :func:`engineeringkit.structural_engineering.johnson_critical_stress`
and :func:`engineeringkit.structural_engineering.effective_length_factor`.
The example draws the design curve, checks that the two branches meet with
the same slope, and compares four end conditions.

*References:* J. B. Johnson, C. W. Bryan and F. E. Turneaure, *The Theory
and Practice of Modern Framed Structures* (Wiley, New York, 1893). The
date follows the secondary literature; the first edition was not
consulted.

.. minigallery:: ../../examples/structural_engineering/stability/plot_04_johnson_parabola.py


1896 -- Bredt's Thin-Walled Torsion
-----------------------------------

Rudolf Bredt, an engineer with a crane-building firm, showed that in a
closed thin-walled tube of any shape the shear flow :math:`q = \tau t` is
the same all round the wall, :math:`q = T/(2A_m)`, where :math:`A_m` is
the area enclosed by the wall's centreline. The twist follows as
:math:`\phi = (TL/4A_m^2G)\oint ds/t`. A closed section is far stiffer in
torsion than the same wall slit open: the ratio is :math:`3(r/t)^2` for
a circular tube. Box girders, wing boxes and monocoque fuselages are
designed with Bredt's two formulas.

*Implementation:* :func:`engineeringkit.structural_engineering.thin_walled_torsion_shear_stress`
and :func:`engineeringkit.structural_engineering.thin_walled_torsion_angle`.
The example checks Bredt against Coulomb for a circular tube and compares
it with the same tube slit along its length.

*References:* R. Bredt, "Kritische Bemerkungen zur Drehungselastizität,"
*Zeitschrift des Vereines deutscher Ingenieure* 40 (1896) 785ff. The page
is from secondary bibliographies.

.. minigallery:: ../../examples/structural_engineering/members/plot_04_bredt_thin_walled_torsion.py


1914 -- Kazinczy's Plastic Hinge
--------------------------------

Gábor Kazinczy loaded steel beams built into concrete walls in Budapest
and found that they carried far more than elastic theory allowed. When
the most stressed section yielded through its whole depth, it went on
turning at a constant plastic moment :math:`M_p`, a *plastic hinge*, and
the load moved to other sections until enough hinges had formed to make a
mechanism. Plastic design of steel frames grew from his tests, through
Maier-Leibnitz in Germany and Baker's group at Cambridge, into the
limit-state codes.

*Implementation:* :func:`engineeringkit.structural_engineering.plastic_section_modulus`
and :func:`engineeringkit.structural_engineering.plastic_collapse_load`.
The example follows a built-in beam under a uniform load from first yield
through hinges at the walls to collapse at :math:`16M_p/L^2`.

*References:* G. Kazinczy, "Kísérletek befalazott tartókkal" (Experiments
with clamped girders), *Betonszemle* 2 (1914); J. F. Baker, M. R. Horne
and J. Heyman, *The Steel Skeleton*, vol. 2 (Cambridge University Press,
1956).

.. minigallery:: ../../examples/structural_engineering/plasticity/plot_01_kazinczy_plastic_hinge.py

1915 -- Maney's Slope-Deflection Method
---------------------------------------

George A. Maney at the University of Minnesota wrote the moment at each
end of a member in terms of the rotations of its ends and of its chord:
:math:`M_{ij} = (2EI/L)(2\theta_i + \theta_j - 3\psi) + M^F_{ij}`.
Moment equilibrium at each free joint then gives one equation per unknown
rotation. Support settlements enter through the chord rotation
:math:`\psi`. Hardy Cross's moment distribution (1930) solves the same
equations by relaxation, and the stiffness method is their matrix form.
Johnson's *Modern Framed Structures* credits Otto Mohr with the idea before
Maney found it independently.

*Implementation:* :func:`engineeringkit.structural_engineering.slope_deflection`
and :class:`engineeringkit.structural_engineering.SlopeDeflectionResult`.
The example analyses a three-span girder, checks it against moment
distribution and lets a pier settle.

*References:* G. A. Maney, "Secondary stresses in rigid frames," *Studies
in Engineering* (University of Minnesota, Minneapolis, 1915). Checked
against secondary accounts only.

.. minigallery:: ../../examples/structural_engineering/frames/plot_03_maney_slope_deflection.py


1919 -- Macaulay's Singularity Functions
----------------------------------------

Integrating :math:`EIy'' = M` piece by piece between loads leaves two
constants of integration per piece, tied together by continuity
conditions. W. H. Macaulay, a Cambridge mathematician, wrote the moment once for
the whole beam with brackets :math:`\langle x - a\rangle^n` that vanish to
the left of :math:`a`. Integrated as ordinary powers, they leave only two
constants. Clebsch had used the idea in 1862, but Macaulay's two-page note
brought it into English textbooks. With the reactions as unknowns too,
statically indeterminate beams are solved the same way as determinate
ones. Superposition of standard cases, the older method, gives the same
answers.

*Implementation:* :class:`engineeringkit.structural_engineering.Beam` and
:func:`engineeringkit.structural_engineering.superpose_beam_results`. The
example solves a propped cantilever under a point load, a couple and a
triangular load, and repeats it by superposition.

*References:* W. H. Macaulay, "A note on the deflection of beams,"
*Messenger of Mathematics* 48 (1919) 129-130; A. Clebsch, *Theorie der
Elasticität fester Körper* (Teubner, Leipzig, 1862).

.. minigallery:: ../../examples/structural_engineering/beams/plot_05_macaulay_singularity_functions.py


1921 -- Maillart's Shear Centre
-------------------------------

Carl von Bach's tests at Stuttgart had shown channel beams twisting under
loads through the web. The Swiss engineer Robert Maillart explained why.
The shear flows in the flanges form a couple, so a load bends the beam
without twisting it only if it passes through a point behind the web, the
*shear centre*, at :math:`e = 3b^2t_f/(6bt_f + ht_w)`. In a section such as
an angle the principal axes are inclined, so a vertical load also bends the
beam sideways. The neutral axis then follows the generalized flexure
formula, not the plane of the load.

*Implementation:* :func:`engineeringkit.structural_engineering.channel_shear_center`,
:func:`engineeringkit.structural_engineering.unsymmetrical_bending_stress`,
:func:`engineeringkit.structural_engineering.neutral_axis_angle`,
:func:`engineeringkit.structural_engineering.channel_section` and
:func:`engineeringkit.structural_engineering.angle_section`. The example
locates the shear centre of a channel from the flange couple and the
neutral axis of an angle under a vertical load.

*References:* R. Maillart, "Zur Frage der Biegung," *Schweizerische
Bauzeitung* 77 (1921) 195-197.

.. minigallery:: ../../examples/structural_engineering/shear/plot_02_maillart_shear_center.py


1930 -- Hardy Cross's Moment Distribution
-----------------------------------------

Continuous beams and building frames led to large sets of simultaneous
equations that no engineer wanted to solve by hand. Hardy Cross proposed
relaxation instead. Lock every joint and apply the fixed-end moments, then
release one joint at a time: its unbalanced moment is shared among the
members in proportion to their stiffness, and half of each share is
carried over to the far end. Repeat until the joints balance. The method
fitted on a sheet of paper, converged in a few cycles and was the main
tool of frame analysis until computers replaced it.

*Implementation:* :func:`engineeringkit.structural_engineering.moment_distribution`.
The example balances a three-span beam built in at one end, checks the
end moments against the stiffness method and plots the decay of the
unbalance.

*References:* H. Cross, "Analysis of continuous frames by distributing
fixed-end moments," *Proc. ASCE* 56 (1930) 919-928.

.. minigallery:: ../../examples/structural_engineering/frames/plot_01_hardy_cross_moment_distribution.py

1938 -- Simmons and Ruge's Bonded Strain Gauge
----------------------------------------------

Stress cannot be measured, only strain. Edward E. Simmons at Caltech (from
1936) and Arthur C. Ruge at MIT (April 1938) independently cemented fine
resistance wire to a surface, so that the wire's resistance changed with
the strain beneath it. The device was sold as the SR-4 gauge. One gauge
gives the strain in one direction. A *rosette* of three at known angles
gives the full plane strain state, and Hooke's law turns it into stress.
Strain gauges made experimental stress analysis a routine laboratory and
field method.

*Implementation:* :func:`engineeringkit.structural_engineering.strain_rosette`
and :func:`engineeringkit.structural_engineering.plane_stress_from_strains`.
The example simulates a rosette on a shaft under bending and torsion and
recovers the torque from three noisy readings.

*References:* E. E. Simmons, US Patent 2,292,549 (11 August 1942). The
dates of Simmons's and Ruge's work are from secondary accounts.

.. minigallery:: ../../examples/structural_engineering/stress/plot_02_strain_gauge_rosette.py

1956 -- Turner, Clough, Martin and Topp: the Direct Stiffness Method
--------------------------------------------------------------------

At Boeing, M. J. Turner, R. W. Clough, H. C. Martin and L. J. Topp needed
the stiffness of swept jet wings, which hand methods could not give. They
wrote a stiffness matrix for each element, added the matrices into one
global matrix according to the joints each element connects, and solved
:math:`K_{ff}u_f = F_f` on a computer. Indeterminate structures cost no
more than determinate ones. The paper is the origin of the finite element
method; Clough gave it that name in 1960. A near-mechanism gives a matrix
that is invertible in floating point but meaningless, so the condition
number of :math:`K_{ff}` is part of the answer.

*Implementation:* :class:`engineeringkit.structural_engineering.PlanarTruss`
(with :math:`\operatorname{cond}(K_{ff})` and a warning near a mechanism).
The example checks a redundant truss against Clapeyron's energy theorem,
then flattens a joint until the warning fires.

*References:* M. J. Turner, R. W. Clough, H. C. Martin and L. J. Topp,
"Stiffness and deflection analysis of complex structures," *J. Aeronaut.
Sci.* 23 (1956) 805-823.

.. minigallery:: ../../examples/structural_engineering/trusses/plot_03_direct_stiffness.py

1956 -- Livesley's Plane Frames by Computer
-------------------------------------------

Moment distribution assumes that the joints do not move sideways. A
building frame under wind does sway, and hand methods needed an extra
correction for every storey. R. K. Livesley programmed the stiffness
method for rigid-jointed frames on the Manchester computer. Each member
got a matrix with axial and bending terms, three degrees of freedom per
joint, and the sway came out of the same solve as the rotations. Matrix
frame analysis became the everyday tool of structural design offices.

*Implementation:* :class:`engineeringkit.structural_engineering.PlanarFrame`.
The example sways a fixed-base portal frame, checks the result against
the rigid-beam and pinned-beam limits and checks the equilibrium of the
reactions under wind and gravity.

*References:* R. K. Livesley, "The application of an electronic digital
computer to some problems of structural analysis," *Struct. Eng.* 34
(1956) 1-12.

.. minigallery:: ../../examples/structural_engineering/frames/plot_02_planar_frame.py

See Also
--------

- :doc:`/api/structural_engineering`
- :doc:`/guides/structural_engineering`
- :doc:`/history/materials_breakthroughs`
- :doc:`/history/mechanical_engineering_breakthroughs`

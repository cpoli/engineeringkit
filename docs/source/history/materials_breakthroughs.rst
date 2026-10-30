Breakthroughs in Materials
==========================


.. include:: /_generated/nav/materials.rst

.. epigraph::

   "Ut tensio, sic vis" (as the extension, so the force).
   -- Robert Hooke, *Lectures de Potentia Restitutiva*, 1678

The strength of materials began as a list of breaking loads and became a
science in three stages, while X-ray diffraction, diffusion and
transformation kinetics opened up the structure inside. The nineteenth century gave elasticity its
constants and its thermal stresses, and gave plasticity its first yield
criterion. Around 1900 engineers learned that parts fail by fatigue and
creep long before they reach their static strength, and wrote empirical
laws for both. The twentieth century explained those laws from the inside:
cracks explained why real solids are weak, dislocations why metals are
soft, and fracture mechanics turned crack size into a design variable. The
entries trace the ideas behind :mod:`engineeringkit.materials`. Each links
to the implementing function and to a gallery example devoted to it.

.. contents:: Timeline
   :local:
   :depth: 1

1807 -- Young's Modulus
-----------------------

Hooke's law said that extension is proportional to force, but the constant
belonged to a particular spring or wire. Thomas Young's *Course of Lectures
on Natural Philosophy and the Mechanical Arts* gave a constant for the
material itself. Young described his "modulus of elasticity" as the height
of a column of the material whose weight would compress it to nothing if
the law stayed linear, a length of thousands of kilometres for steel.
Navier (1826) restated it as stress over strain, :math:`E = \sigma/\varepsilon`,
the form every tensile test now reports.

*Implementation:* :func:`engineeringkit.materials.hooke_stress`,
:func:`engineeringkit.materials.normal_stress`,
:func:`engineeringkit.materials.normal_strain`,
:func:`engineeringkit.materials.resilience_modulus`,
:func:`engineeringkit.materials.toughness` and
:class:`engineeringkit.materials.StressStrainCurve`. The example reduces a
simulated tensile test to a stress-strain curve, fits the modulus and
measures resilience and toughness.

*References:* T. Young, *A Course of Lectures on Natural Philosophy and the
Mechanical Arts* (Johnson, London, 1807), vol. 1, lecture XIII; C.-L.
Navier, *Résumé des leçons données à l'École des Ponts et Chaussées sur
l'application de la mécanique* (Didot, Paris, 1826).

.. minigallery:: ../../examples/materials/elasticity/plot_01_young_modulus.py

1838 -- Duhamel's Thermal Stress
--------------------------------

Jean-Marie Duhamel extended the equations of elasticity to bodies heated
unevenly. A temperature change :math:`\Delta T` produces a free strain
:math:`\alpha\,\Delta T`; where the body cannot expand, an elastic strain
cancels it and a stress :math:`-E\alpha\,\Delta T` appears. The stress does
not depend on the length of the bar, so a long weld is no worse than a short
one, but it is large: 40 K in restrained steel gives about 100 MPa. Franz
Neumann reached the same equations independently in 1841.

*Implementation:* :func:`engineeringkit.materials.thermal_strain` and
:func:`engineeringkit.materials.constrained_thermal_stress`. The example
computes the stress in continuously welded and jointed rail.

*References:* J.-M.-C. Duhamel, "Mémoire sur le calcul des actions
moléculaires développées par les changements de température dans les corps
solides," *Mém. Divers Savants Acad. R. Sci. Inst. France* 5 (1838)
440-498.

.. minigallery:: ../../examples/materials/elasticity/plot_02_duhamel_thermal_stress.py

1852 -- Lamé's Elastic Constants
--------------------------------

Navier and Poisson built elasticity from forces between molecules and
concluded that an isotropic solid has one elastic constant, with Poisson's
ratio fixed at 1/4. Cauchy and Green, working from the continuum, found two,
and measurements on rubber and metals sided with them. Gabriel Lamé's
*Leçons*, the first book on the mathematical theory of elasticity, wrote the
isotropic law with the two constants :math:`\lambda` and :math:`\mu`,
:math:`\sigma_{ij} = \lambda\,\varepsilon_{kk}\delta_{ij} + 2\mu\,\varepsilon_{ij}`.
Any two of :math:`E`, :math:`\nu`, :math:`G`, :math:`K`, :math:`\lambda`
give the rest.

*Implementation:* :func:`engineeringkit.materials.lame_parameters`,
:func:`engineeringkit.materials.shear_modulus` and
:func:`engineeringkit.materials.bulk_modulus`. The example recovers
uniaxial tension from Lamé's law and plots the moduli against Poisson's
ratio.

*References:* G. Lamé, *Leçons sur la théorie mathématique de l'élasticité
des corps solides* (Bachelier, Paris, 1852).

.. minigallery:: ../../examples/materials/elasticity/plot_03_lame_constants.py

1855 -- Fick's Laws of Diffusion
--------------------------------

Adolf Fick, prosector in anatomy at Zürich, measured salt diffusing up
tall columns of water and described it with Fourier's law of heat
conduction: the flux is proportional to the concentration gradient,
:math:`J = -D\,\partial C/\partial x`, and conservation gives
:math:`\partial C/\partial t = D\,\partial^2 C/\partial x^2`. Roberts-Austen
measured diffusion of gold in solid lead in 1896, and Arrhenius-type
diffusivities, :math:`D = D_0 e^{-Q_d/RT}`, made Fick's second law the tool
for timing heat treatments. Its error-function solution for a surface held
at fixed concentration sets the furnace time for case carburizing.

*Implementation:* :func:`engineeringkit.materials.fick_flux`,
:func:`engineeringkit.materials.diffusion_coefficient`,
:func:`engineeringkit.materials.semi_infinite_concentration` and
:func:`engineeringkit.materials.diffusion_time`. The example computes a
steady flux and the carburizing time of a gear at three temperatures.

*References:* A. Fick, "Ueber Diffusion," *Ann. Phys.* 170 (1855) 59-86;
W. C. Roberts-Austen, "On the diffusion of metals," *Phil. Trans. R. Soc.
Lond. A* 187 (1896) 383-415.

.. minigallery:: ../../examples/materials/diffusion/plot_01_fick_carburizing.py

1864 -- Tresca's Maximum-Shear Yield Criterion
----------------------------------------------

Henri Tresca, at the Conservatoire des Arts et Métiers in Paris, punched and
extruded lead, copper and iron and found that solid metal flows like a
liquid once the largest shear stress reaches a constant. With principal
stresses :math:`\sigma_1 \ge \sigma_2 \ge \sigma_3`, yield is
:math:`\sigma_1 - \sigma_3 = \sigma_y`. Saint-Venant (1870) built the first
theory of plastic flow on it. It is the first criterion for yielding under
combined stresses, and it remains the conservative choice in pressure
vessel codes.

*Implementation:* :func:`engineeringkit.materials.tresca_equivalent` and
:func:`engineeringkit.materials.factor_of_safety`. The example checks a
shaft in bending and torsion and plots Tresca's hexagon.

*References:* H. Tresca, "Mémoire sur l'écoulement des corps solides soumis
à de fortes pressions," *C. R. Acad. Sci. Paris* 59 (1864) 754-758.

.. minigallery:: ../../examples/materials/yield/plot_01_tresca_criterion.py

1889-1929 -- Voigt and Reuss: Bounds on a Composite
---------------------------------------------------

Woldemar Voigt (1889) estimated the elastic constants of a polycrystal by
assuming that all its crystals strain alike, which averages their
stiffnesses. Adolf Reuss (1929) assumed instead that they all carry the
same stress, which averages their compliances. For a fibre composite the
two assumptions are exact along and across the fibres: the rules of
mixtures. Rodney Hill (1952) proved that they bound the modulus of every
arrangement of the phases, and Hashin and Shtrikman (1963) narrowed the
bounds for isotropic mixtures.

*Implementation:* :func:`engineeringkit.materials.rule_of_mixtures_E_longitudinal`,
:func:`engineeringkit.materials.rule_of_mixtures_E_transverse` and
:func:`engineeringkit.materials.density_rule_of_mixtures`. The example
builds both bounds from springs and plots the band between them for glass
fibre in epoxy.

*References:* W. Voigt, "Ueber die Beziehung zwischen den beiden
Elasticitätsconstanten isotroper Körper," *Ann. Phys.* 274 (1889) 573-587;
A. Reuss, "Berechnung der Fließgrenze von Mischkristallen auf Grund der
Plastizitätsbedingung für Einkristalle," *Z. Angew. Math. Mech.* 9 (1929)
49-58; R. Hill, "The elastic behaviour of a crystalline aggregate," *Proc.
Phys. Soc. A* 65 (1952) 349-354.

.. minigallery:: ../../examples/materials/composites/plot_01_voigt_reuss_bounds.py

1899 -- Goodman's Mean-Stress Diagram
-------------------------------------

Wöhler's tests showed that a steady stress lowers the alternating stress a
part can survive. Gerber (1874) fitted a parabola through Wöhler's data.
John Goodman, professor of engineering at Leeds, replaced it in his
textbook *Mechanics Applied to Engineering* with a straight line from the
endurance limit at zero mean stress to the ultimate strength at zero
amplitude, :math:`\sigma_a/S_e + \sigma_m/S_{ut} = 1/n`. The line is slightly
conservative, easy to use, and still the default in machine design.

*Implementation:* :func:`engineeringkit.materials.goodman_safety_factor`,
with :func:`engineeringkit.materials.gerber_safety_factor` and
:func:`engineeringkit.materials.soderberg_safety_factor` for the rival
lines. The example finds the safety factor of a preloaded bolt and compares the
Gerber, Goodman and Soderberg lines.

*References:* J. Goodman, *Mechanics Applied to Engineering* (Longmans,
Green, London, 1899).

.. minigallery:: ../../examples/materials/fatigue/plot_01_goodman_diagram.py

1900 -- Brinell's Hardness Test
-------------------------------

Johan August Brinell, chief engineer at the Fagersta ironworks, needed a
quick check of the steels he was heat-treating. He pressed a 10 mm
hardened steel ball into the surface under a fixed load and divided the
load by the curved area of the dent. He presented the test at the
International Congress for Testing Materials at the 1900 Paris
Exposition, and Axel Wahlberg's English account (1901) spread it. Because
an indentation is constrained plastic flow, the number tracks the tensile
strength, :math:`S_{ut} \approx 3.45\,\mathrm{HB}` MPa for steels, so a
small mark on a finished part estimates its strength. Smith and
Sandland's Vickers diamond (1922) extended the scale to hardened steel.

*Implementation:* :func:`engineeringkit.materials.brinell_hardness`,
:func:`engineeringkit.materials.brinell_tensile_strength` and
:func:`engineeringkit.materials.vickers_hardness`. The example converts a
dent to HB and HB to strength across the range of steels.

*References:* J. A. Brinell, "Mémoire sur les épreuves à bille en acier,"
*Congrès International des Méthodes d'Essai des Matériaux de
Construction* (Paris, 1900), vol. 2, 83-94; A. Wahlberg, "Brinell's method
of determining hardness and other properties of iron and steel," *J. Iron
Steel Inst.* 59 (1901) 243-298. The page ranges were taken from secondary
sources.

.. minigallery:: ../../examples/materials/hardness/plot_01_brinell_hardness.py

1910 -- Basquin's S-N Law
-------------------------

August Wöhler's railway-axle tests of the 1860s established the stress-life
(S-N) curve and the endurance limit of steel. Olin Basquin, of Northwestern
University, noticed that on log-log axes the finite-life part of the
curve is a straight line, :math:`\sigma_a = \sigma_f'(2N)^b`, with
:math:`b` around -0.1. A fixed fraction of strength is lost for each decade
of life, and two parameters describe a material's whole high-cycle fatigue
behaviour.

*Implementation:* :func:`engineeringkit.materials.basquin_cycles` and
:func:`engineeringkit.materials.utils.power_law_fit`. The example fits the
law to scattered endurance tests and recovers its constants.

*References:* O. H. Basquin, "The exponential law of endurance tests,"
*Proc. Am. Soc. Test. Mater.* 10 (1910) 625-630.

.. minigallery:: ../../examples/materials/fatigue/plot_02_basquin_sn_curve.py

1913 -- The Braggs' Law of X-ray Diffraction
--------------------------------------------

Max von Laue, Walter Friedrich and Paul Knipping showed in 1912 that
crystals diffract X-rays. William Lawrence Bragg, a research student at
Cambridge, explained the spots as reflections from planes of atoms that
reinforce when :math:`n\lambda = 2d\sin\theta`; his paper was read in
November 1912 and printed in 1913. With the X-ray spectrometer of his
father, William Henry Bragg, the two solved the structures of rock salt,
diamond and the metals, and showed that iron is body-centred and copper
face-centred cubic. The packing of spheres in those cells fixes the
density of a perfect crystal. Father and son shared the 1915 Nobel Prize
in Physics.

*Implementation:* :func:`engineeringkit.materials.bragg_angle`,
:func:`engineeringkit.materials.interplanar_spacing`,
:func:`engineeringkit.materials.diffraction_peaks`,
:func:`engineeringkit.materials.cubic_lattice_parameter`,
:func:`engineeringkit.materials.atomic_packing_factor` and
:func:`engineeringkit.materials.theoretical_density`. The example builds
the powder patterns of iron and copper and reads copper's density back
from its first peak.

*References:* W. L. Bragg, "The diffraction of short electromagnetic waves
by a crystal," *Proc. Camb. Phil. Soc.* 17 (1913) 43-57; W. H. Bragg and
W. L. Bragg, "The reflection of X-rays by crystals," *Proc. R. Soc. Lond.
A* 88 (1913) 428-438.

.. minigallery:: ../../examples/materials/crystals/plot_01_bragg_diffraction.py

1913 -- Von Mises's Yield Criterion
-----------------------------------

Richard von Mises replaced Tresca's hexagon with the smooth surface through
its corners, yield when the second invariant of the deviatoric stress
reaches :math:`\sigma_y^2/3`. Maksymilian Huber had proposed the same rule in
1904, and Heinrich Hencky (1924) explained it as a critical elastic energy
of distortion. The criterion predicts shear yield at
:math:`\sigma_y/\sqrt3`, 15% above Tresca's :math:`\sigma_y/2`, and Taylor and
Quinney's tension-torsion tests on tubes of copper, aluminium and steel
(1931) fell closer to it.

*Implementation:* :func:`engineeringkit.materials.von_mises_plane_stress`,
:func:`engineeringkit.materials.von_mises_stress`,
:func:`engineeringkit.materials.principal_stresses` and
:func:`engineeringkit.materials.visualizers.plot_yield_surfaces`. The
example checks the distortion-energy form, shows that a hydrostatic stress
leaves a three-dimensional state's equivalent stress unchanged, and
compares the two criteria for a tube in tension and torsion.

*References:* R. von Mises, "Mechanik der festen Körper im
plastisch-deformablen Zustand," *Nachr. Ges. Wiss. Göttingen,
Math.-Phys. Kl.* (1913) 582-592; G. I. Taylor and H. Quinney, "The plastic
distortion of metals," *Phil. Trans. R. Soc. Lond. A* 230 (1931) 323-362.

.. minigallery:: ../../examples/materials/yield/plot_02_von_mises_criterion.py

1920 -- Griffith's Theory of Brittle Fracture
---------------------------------------------

Alan Arnold Griffith, at the Royal Aircraft Establishment, asked why glass
breaks at a hundredth of the strength its atomic bonds imply. His answer,
read to the Royal Society in 1920 and published in 1921, was that small
cracks concentrate the stress, and that a crack grows when the elastic
energy it releases exceeds the surface energy it creates. Strength then
falls as the inverse square root of crack size,
:math:`\sigma_f = \sqrt{2E\gamma_s/\pi a}`. Griffith confirmed it on scratched
glass bulbs and tubes, and freshly drawn glass fibres came close to the
theoretical strength.

*Implementation:* :func:`engineeringkit.materials.griffith_fracture_stress`.
The example finds the maximum of the energy balance and plots strength
against flaw size for glass.

*References:* A. A. Griffith, "The phenomena of rupture and flow in
solids," *Phil. Trans. R. Soc. Lond. A* 221 (1921) 163-198.

.. minigallery:: ../../examples/materials/fracture/plot_01_griffith_fracture.py

1929 -- Norton's Power-Law Creep
--------------------------------

Steam turbines and cracking stills of the 1920s ran hot enough for their
steel to stretch slowly under constant load. Frederick Norton, at MIT,
measured creep rates of steels at temperature and found the steady rate
rising as a power of stress, :math:`\dot\varepsilon = A\sigma^n`, with
:math:`n` around 4 to 6. In the 1950s Dorn and co-workers showed that the
temperature dependence is an Arrhenius factor :math:`e^{-Q/RT}` with an
activation energy close to that of self-diffusion, which tied creep to the
climb of dislocations.

*Implementation:* :func:`engineeringkit.materials.arrhenius_creep_rate`.
The example recovers the stress exponent and the activation energy from
simulated creep tests.

*References:* F. H. Norton, *The Creep of Steel at High Temperatures*
(McGraw-Hill, New York, 1929).

.. minigallery:: ../../examples/materials/creep/plot_01_norton_creep.py

1934 -- Taylor, Orowan and Polanyi: Dislocations
------------------------------------------------

Jacob Frenkel (1926) estimated that shearing a perfect crystal needs about
:math:`G/2\pi`, thousands of times more than metal single crystals need. In
1934 Geoffrey Taylor, Egon Orowan and Michael Polanyi independently
proposed the dislocation, a line defect whose glide moves one row of atoms
at a time, so that slip needs only a small stress. Rudolf Peierls (1940) and
Frank Nabarro (1947) estimated the lattice's resistance to it, which falls
exponentially as the slip planes get farther apart. Dislocations were seen
directly in the electron microscope in 1956.

*Implementation:* :func:`engineeringkit.materials.ideal_shear_strength` and
:func:`engineeringkit.materials.peierls_stress`. The example compares
Frenkel's estimate, the Peierls stress and measured strengths of copper and
aluminium.

*References:* J. Frenkel, "Zur Theorie der Elastizitätsgrenze und der
Festigkeit kristallinischer Körper," *Z. Phys.* 37 (1926) 572-609; G. I.
Taylor, "The mechanism of plastic deformation of crystals. Part I.
Theoretical," *Proc. R. Soc. Lond. A* 145 (1934) 362-387; E. Orowan, "Zur
Kristallplastizität. III," *Z. Phys.* 89 (1934) 634-659; M. Polanyi, "Über
eine Art Gitterstörung, die einen Kristall plastisch machen könnte," *Z.
Phys.* 89 (1934) 660-664; R. Peierls, "The size of a dislocation," *Proc.
Phys. Soc.* 52 (1940) 34-37; F. R. N. Nabarro, "Dislocations in a simple
cubic lattice," *Proc. Phys. Soc.* 59 (1947) 256-272.

.. minigallery:: ../../examples/materials/plasticity/plot_01_dislocations.py

1939 -- Avrami's Transformation Kinetics
----------------------------------------

Recrystallization, pearlite growing from austenite and the crystallization
of polymers all follow an S-shaped curve in time. Melvin Avrami took
nuclei that appear at random and grow until they impinge, counted the
volume they would fill if they could overlap, and removed the overlap
statistically: the fraction transformed is :math:`y = 1 - e^{-kt^n}`, with
:math:`n` from 1 to 4 set by the modes of nucleation and growth.
Kolmogorov (1937) and Johnson and Mehl (1939) derived the same law, so it
is also called JMAK. Its isothermal curves, read at 1% and 99%, are the
time-temperature-transformation diagrams of heat treatment, and the lever
rule on the phase diagram gives the amount they approach.

*Implementation:* :func:`engineeringkit.materials.avrami_fraction`,
:func:`engineeringkit.materials.avrami_time`,
:func:`engineeringkit.materials.avrami_fit` and
:func:`engineeringkit.materials.lever_rule`. The example simulates random
nuclei on a grid, recovers :math:`n \approx 3`, and finds the pearlite
fraction of a 0.40% carbon steel.

*References:* M. Avrami, "Kinetics of phase change. I. General theory,"
*J. Chem. Phys.* 7 (1939) 1103-1112; A. N. Kolmogorov, "On the statistical
theory of the crystallization of metals," *Izv. Akad. Nauk SSSR, Ser.
Mat.* 1 (1937) 355-359; W. A. Johnson and R. F. Mehl, "Reaction kinetics
in processes of nucleation and growth," *Trans. AIME* 135 (1939) 416-442.

.. minigallery:: ../../examples/materials/phases/plot_01_avrami_kinetics.py

1939 -- Weibull's Strength of Brittle Materials
-----------------------------------------------

Waloddi Weibull, professor at the Royal Institute of Technology in
Stockholm, explained why brittle solids scatter in strength and why
large ones are weaker. A part fails at its weakest flaw, so its survival
probability is the product of those of its volume elements,
:math:`1 - P_f = \exp[-(V/V_0)(\sigma/\sigma_0)^m]`. The Weibull modulus
:math:`m` measures the scatter, and a part 1000 times the size of a test
bar is weaker by :math:`1000^{1/m}`. Weibull generalised the distribution
in 1951, and it became the standard lifetime model of reliability
engineering (:doc:`/history/reliability_breakthroughs`).

*Implementation:* :func:`engineeringkit.materials.weibull_failure_probability`,
:func:`engineeringkit.materials.weibull_size_effect` and
:func:`engineeringkit.materials.weibull_mean_strength`, which evaluate
:func:`engineeringkit.reliability.weibull_cdf` and
:func:`engineeringkit.reliability.weibull_mttf`. The example fits the
modulus to bend-bar strengths and scales the strength to a larger part.

*References:* W. Weibull, "A statistical theory of the strength of
materials," *Ingeniörsvetenskapsakademiens Handlingar* 151 (Stockholm,
1939).

.. minigallery:: ../../examples/materials/fracture/plot_04_weibull_strength.py

1943 -- Ramberg and Osgood's Stress-Strain Curve
------------------------------------------------

Aluminium alloys and stainless steels yield gradually, with no sharp yield
point, and wartime aircraft designers needed their curves in a form they
could compute with. Walter Ramberg and William Osgood, at the National
Bureau of Standards for NACA, described the whole curve by three
parameters, :math:`\varepsilon = \sigma/E + K(\sigma/E)^n`. Written with the
0.2% offset yield strength, :math:`\varepsilon = \sigma/E + 0.002(\sigma/\sigma_y)^n`,
it remains the standard description of hardening in buckling and fatigue
analysis.

*Implementation:* :func:`engineeringkit.materials.ramberg_osgood_strain`
and :meth:`engineeringkit.materials.StressStrainCurve.offset_yield_strength`.
The example plots the family of curves and recovers the offset yield
strength from sampled data.

*References:* W. Ramberg and W. R. Osgood, "Description of stress-strain
curves by three parameters," NACA Technical Note 902 (1943).

.. minigallery:: ../../examples/materials/plasticity/plot_02_ramberg_osgood.py

1945 -- Palmgren and Miner's Cumulative Damage
----------------------------------------------

Arvid Palmgren (1924), estimating the lives of ball bearings under varying
loads, assumed that each cycle at amplitude :math:`\sigma_i` uses up
:math:`1/N_i` of the life. Milton Miner, at Douglas Aircraft, tested riveted
aluminium joints under mixed amplitudes and found failure near
:math:`\sum n_i/N_i = 1`. The rule ignores the order of the loads, which
tests show matters, but it is simple and turned load spectra into fatigue
lives for every aircraft since.

*Implementation:* :func:`engineeringkit.materials.miner_damage`. The example
sums the damage of a block spectrum and reduces it to an equivalent
constant amplitude.

*References:* A. Palmgren, "Die Lebensdauer von Kugellagern," *Z. Ver.
Dtsch. Ing.* 68 (1924) 339-341; M. A. Miner, "Cumulative damage in
fatigue," *J. Appl. Mech.* 12 (1945) A159-A164.

.. minigallery:: ../../examples/materials/fatigue/plot_03_palmgren_miner.py

1951-1953 -- Hall and Petch: Grain-Size Strengthening
-----------------------------------------------------

Eric Hall measured the yield stress of mild steels with different grain
sizes, and Norman Petch the cleavage strength of iron. Both found strength
rising linearly with the inverse square root of grain size,
:math:`\sigma_y = \sigma_0 + k\,d^{-1/2}`. Dislocations pile up at grain
boundaries, and the stress at the head of a pile-up grows with its length.
Grain refinement, by controlled rolling, became the way to make steels
both stronger and tougher.

*Implementation:* :func:`engineeringkit.materials.hall_petch`. The example
fits the Hall-Petch line to scattered yield strengths.

*References:* E. O. Hall, "The deformation and ageing of mild steel: III
Discussion of results," *Proc. Phys. Soc. B* 64 (1951) 747-753; N. J.
Petch, "The cleavage strength of polycrystals," *J. Iron Steel Inst.* 174
(1953) 25-28.

.. minigallery:: ../../examples/materials/plasticity/plot_03_hall_petch.py

1952 -- The Larson-Miller Parameter
-----------------------------------

Power plants are designed for 100,000 hours of creep, longer than anyone
can test. Frank Larson and James Miller, of General Electric, showed that
rupture times at different temperatures collapse onto one curve of stress
against :math:`P = T(C + \log_{10} t_r)`, with :math:`C` about 20 for hours.
The parameter follows from an Arrhenius rupture time, and lets a short test
at a higher temperature stand in for a long life at a lower one. Sherby
and Dorn (1954) and Manson and Haferd (1953) proposed rival parameters, but
Larson-Miller remains the one in design codes.

*Implementation:* :func:`engineeringkit.materials.larson_miller_parameter`
and :func:`engineeringkit.materials.larson_miller_rupture_time`. The example
collapses rupture tests onto a master curve and finds the short test
equivalent to a 100,000-hour life.

*References:* F. R. Larson and J. Miller, "A time-temperature relationship
for rupture and creep stresses," *Trans. ASME* 74 (1952) 765-775.

.. minigallery:: ../../examples/materials/creep/plot_02_larson_miller.py

1954 -- Coffin and Manson: Low-Cycle Fatigue
--------------------------------------------

Parts of nuclear plant and jet engines heated and cooled a few thousand
times fail by cycles of plastic strain, which Basquin's stress-life law
cannot describe. Louis Coffin, at General Electric, and Samson Manson, at
NACA's Lewis laboratory, independently found that the plastic strain
amplitude follows a power law of life,
:math:`\Delta\varepsilon_p/2 = \varepsilon_f'(2N)^c`, with :math:`c` near
-0.6. Added to Basquin's elastic line it gives the strain-life curve used
for notched parts and thermal fatigue.

*Implementation:* :func:`engineeringkit.materials.coffin_manson_cycles` and
:func:`engineeringkit.materials.strain_life_amplitude`. The example draws
the strain-life curve, finds its transition life and solves for the life
at a given strain.

*References:* L. F. Coffin Jr., "A study of the effects of cyclic thermal
stresses on a ductile metal," *Trans. ASME* 76 (1954) 931-950; S. S.
Manson, "Behavior of materials under conditions of thermal stress," NACA
Report 1170 (1954).

.. minigallery:: ../../examples/materials/fatigue/plot_04_coffin_manson.py

1957 -- Irwin's Stress Intensity Factor
---------------------------------------

Griffith's energy balance predicted absurdly low strengths for steel,
whose cracks blunt by plastic flow, and the brittle fractures of welded
Liberty ships and tankers needed explaining. George Irwin, at the US Naval
Research Laboratory, added the plastic work to Griffith's surface energy
(1948) and then showed that the stresses near every crack tip have one
shape, scaled by the stress intensity factor :math:`K_I = Y\sigma\sqrt{\pi a}`.
Fracture occurs when :math:`K_I` reaches the toughness :math:`K_{Ic}`, a
material property that can be measured. Linear-elastic fracture mechanics
was born.

*Implementation:* :func:`engineeringkit.materials.stress_intensity` and
:func:`engineeringkit.materials.critical_crack_length`. The example shows
that Irwin's criterion reproduces Griffith's for glass and compares the
crack tolerance of two steels.

*References:* G. R. Irwin, "Analysis of stresses and strains near the end
of a crack traversing a plate," *J. Appl. Mech.* 24 (1957) 361-364.

.. minigallery:: ../../examples/materials/fracture/plot_02_irwin_stress_intensity.py

1961 -- Neuber's Rule
---------------------

Kirsch (1898) found the factor of 3 at a circular hole and Inglis (1913)
the factor :math:`1 + 2a/b` at an elliptical one, but in fatigue sharp
notches hurt less than these elastic factors suggest, and in ductile
parts the notch root yields. Peterson and Neuber described the first
effect with a notch sensitivity, :math:`K_f = 1 + q(K_t - 1)`. Heinz
Neuber, at the Technical University of Munich, solved the second for
shear-loaded notches of any nonlinear material and found that the stress
and strain concentrations keep :math:`K_\sigma K_\varepsilon = K_t^2`.
With the cyclic stress-strain curve, :math:`\sigma\varepsilon = (K_t S)^2/E`
gives the local strain that strain-life fatigue analysis needs, and it is
used for notches under any loading.

*Implementation:* :func:`engineeringkit.materials.neuber_notch_stress`,
:func:`engineeringkit.materials.elliptical_hole_concentration`,
:func:`engineeringkit.materials.peterson_notch_sensitivity`,
:func:`engineeringkit.materials.neuber_notch_sensitivity` and
:func:`engineeringkit.materials.fatigue_notch_factor`. The example finds
the fatigue notch factor of a hole and intersects the Neuber hyperbola
with a Ramberg-Osgood curve.

*References:* H. Neuber, "Theory of stress concentration for
shear-strained prismatical bodies with arbitrary nonlinear stress-strain
law," *J. Appl. Mech.* 28 (1961) 544-550; C. E. Inglis, "Stresses in a
plate due to the presence of cracks and sharp corners," *Trans. Inst.
Naval Archit.* 55 (1913) 219-241.

.. minigallery:: ../../examples/materials/fatigue/plot_05_neuber_rule.py

1961 -- Reissner and Stavsky: Classical Laminate Theory
-------------------------------------------------------

Fibre-reinforced plastics are stiff along the fibres and soft across, so
parts are laid up from plies at several angles. Eric Reissner and Yehuda
Stavsky, at MIT, extended Kirchhoff's plate theory to such stacks and
showed that a laminate that is not symmetric about its midplane couples
stretching to bending: the :math:`\mathbf B` matrix of
:math:`[\mathbf N; \mathbf M] = [\mathbf A, \mathbf B; \mathbf B, \mathbf D][\boldsymbol\varepsilon^0; \boldsymbol\kappa]`.
Pister and Dong (1959) had treated layered plates earlier. Summing each
ply's rotated stiffness through the thickness gives the laminate's
response and every ply's stress, the starting point of composite design.

*Implementation:* :func:`engineeringkit.materials.lamina_stiffness`,
:func:`engineeringkit.materials.transformed_stiffness`,
:func:`engineeringkit.materials.classical_laminate` and
:class:`engineeringkit.materials.Laminate`. The example checks a
quasi-isotropic carbon-epoxy laminate and shows a cross-ply one curling
under tension.

*References:* E. Reissner and Y. Stavsky, "Bending and stretching of
certain types of heterogeneous aeolotropic elastic plates," *J. Appl.
Mech.* 28 (1961) 402-408; K. S. Pister and S. B. Dong, "Elastic bending of
layered plates," *J. Eng. Mech. Div. ASCE* 85 (1959) 1-10.

.. minigallery:: ../../examples/materials/composites/plot_03_classical_laminate.py

1963 -- Paris's Crack-Growth Law
--------------------------------

The Comet airliner disasters of 1954 showed that fatigue cracks grow in
service. Paul Paris, with Gomez and Anderson (1961), proposed that the
growth per cycle depends on the range of the stress intensity factor alone,
:math:`da/dN = C\,\Delta K^m`. Journals rejected the idea at first; Paris and
Erdogan's 1963 paper, comparing it with the laws then in use on aluminium
test data, made it standard. Integrated from the largest crack an
inspection can miss to the critical size, it sets the inspection intervals
of damage-tolerant aircraft design.

*Implementation:* :func:`engineeringkit.materials.paris_crack_growth`,
:func:`engineeringkit.materials.crack_growth_life` and
:class:`engineeringkit.materials.CrackGrowthResult`. The example fits the
law to growth-rate data and integrates the life of a cracked panel with a
finite-width correction.

*References:* P. C. Paris, M. P. Gomez and W. E. Anderson, "A rational
analytic theory of fatigue," *The Trend in Engineering* 13 (1961) 9-14;
P. Paris and F. Erdogan, "A critical analysis of crack propagation laws,"
*J. Basic Eng.* 85 (1963) 528-533.

.. minigallery:: ../../examples/materials/fracture/plot_03_paris_law.py

1969 -- The Halpin-Tsai Equations
---------------------------------

The Voigt and Reuss bounds lie far apart for fibre composites, and the
elasticity solutions for regular fibre arrays were too complicated for
design. John Halpin and Stephen Tsai, at the US Air Force Materials
Laboratory, interpolated those solutions with one formula,
:math:`P/P_m = (1 + \xi\eta V_f)/(1 - \eta V_f)`, whose parameter
:math:`\xi` runs from the Reuss bound at zero to the Voigt bound at
infinity: about 2 for the transverse modulus of round fibres, 1 for the
shear modulus and :math:`2\ell/d` for short fibres. The report was never
published in a journal; Halpin and Kardos reviewed the equations in 1976.

*Implementation:* :func:`engineeringkit.materials.halpin_tsai`. The
example estimates the transverse and shear moduli of glass-epoxy and the
stiffness of aligned short fibres.

*References:* J. C. Halpin and S. W. Tsai, "Effects of environmental
factors on composite materials," AFML-TR-67-423, Air Force Materials
Laboratory (1969); J. C. Halpin and J. L. Kardos, "The Halpin-Tsai
equations: a review," *Polym. Eng. Sci.* 16 (1976) 344-352.

.. minigallery:: ../../examples/materials/composites/plot_02_halpin_tsai.py

See Also
--------

- :doc:`/api/materials`
- :doc:`/guides/materials`
- :doc:`/history/structural_engineering_breakthroughs`
- :doc:`/history/reliability_breakthroughs`

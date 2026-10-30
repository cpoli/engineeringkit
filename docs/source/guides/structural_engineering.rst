Structural engineering: loads, signs and linear elasticity
==========================================================

Elastic formulas use modulus E and G in Pa, lengths in m, areas in m^2,
second/polar moments in m^4, force in N and moment/torque in N m. Strain is
dimensionless. Isotropic Poisson ratio must satisfy -1 < nu < 0.5. Euler
buckling assumes a slender ideal elastic column, with effective length K*L;
it does not model yielding, imperfections or post-buckling.

Beam convention
---------------

Positive point loads and uniform loads act downward. Reactions act upward,
sagging moment is positive, and downward deflection is negative. A load at a
support has zero span deflection. Shear at a point load uses the right-side
value; plotted endpoints do not include the support reaction jumps. The
Euler-Bernoulli relation is E*I*y''=M with this sign convention. Small slopes,
constant E*I, and negligible shear deformation are assumed.

.. code-block:: python

    from engineeringkit.structural_engineering import SimplySupportedBeam
    result = SimplySupportedBeam(6, 200e9, 1e-5).solve_point_load(12000, 2)
    assert result.reactions == (8000, 4000)
    fig, axes = result.plot()

The tests check force/moment equilibrium, zero support deflections, numerical
curvature versus M/(EI), and the uniform-load midspan result
-5*w*L^4/(384*E*I). Axial strain energy U=P^2*L/(2*A*E) is differentiated
independently to recover axial displacement.

Beams on any supports
---------------------

``Beam(length, E, I, supports)`` takes supports as ``(position, kind)`` with
kind ``"pin"``, ``"roller"`` (both restrain deflection only) or ``"fixed"``
(deflection and slope), anywhere on the beam, so cantilevers, overhangs,
propped cantilevers and continuous beams are all one call. ``solve`` accepts
point loads ``(x, P)`` downward positive, couples ``(x, C)`` clockwise
positive, and distributed loads ``(start, end, w)`` or
``(start, end, w_start, w_end)`` varying linearly, downward positive. The
moment is written with Macaulay brackets and integrated exactly; the
reactions and the two constants of integration come from one linear system
of equilibrium and support conditions, so indeterminate beams need no extra
input. A support layout that allows rigid-body motion raises a ValueError.
The result is a ``BeamResult`` with one reaction per support, the reaction
moments of the fixed supports (clockwise positive) and the slope. Shear at a
load takes the value just to its right, except at the right end.
``superpose_beam_results`` adds the responses of one beam to several load
cases; it requires the same stations and supports.

.. code-block:: python

    from engineeringkit.structural_engineering import Beam
    beam = Beam(6.0, 200e9, 8e-6, ((0.0, "fixed"), (6.0, "roller")))
    result = beam.solve(distributed_loads=[(0.0, 6.0, 1e3)])
    result.reactions        # (3750.0, 2250.0): 5wL/8 and 3wL/8
    result.support_moments  # (-4500.0,): -wL^2/8

Sections, unsymmetrical bending and shear
-----------------------------------------

Cross-sections use Hibbeler's and Gere's axes: ``y`` vertical (up), ``z``
horizontal, beam axis ``x``. ``Iz`` is the second moment about the
horizontal centroidal axis, ``Iy`` about the vertical one and ``Iyz`` the
product of inertia. ``composite_section`` builds a ``SectionProperties`` from
rectangles ``(b, h, z, y[, weight])``, with ``(z, y)`` the rectangle's
centre; weight -1 removes a hole and a modular ratio ``n = E/E_ref`` gives a
transformed section, whose properties are in the reference material. The
shapes ``i_section``, ``t_section``, ``channel_section``, ``angle_section``
and ``box_section`` use overall dimensions, ignore fillets and corner radii,
and put the bottom of the section at ``y = 0``.
``SectionProperties.first_moment(y)`` is Q of the part above a cut, and
``width_at(y)`` the true width there (holes subtracted, modular ratios
ignored).

``unsymmetrical_bending_stress`` is the generalized flexure formula. Positive
``Mz`` compresses the fibres at positive ``y``, so for a horizontal beam with
downward loads ``Mz`` is the sagging moment; ``My`` follows the right-hand
rule about ``y``. Loads are assumed to pass through the shear centre.
``neutral_axis_angle`` is measured from ``z``, counterclockwise.
``transverse_shear_stress`` is Jourawski's VQ/(It), averaged across the
cut; it is accurate for narrow cuts and thin walls, not across wide
flanges. ``channel_shear_center`` uses centreline dimensions of a
thin-walled channel and returns the distance behind the web centreline.

Torsion of thin-walled tubes, pressure vessels and stress at a point
--------------------------------------------------------------------

Bredt's formulas apply to closed thin-walled sections free to warp:
``thin_walled_torsion_shear_stress`` takes the area enclosed by the wall's
centreline and the local thickness; ``thin_walled_torsion_angle`` takes the
integral of ds/t around the centreline. Open sections are not covered.
``thin_wall_vessel_stresses`` uses the inner radius and neglects the radial
stress (valid for r/t of about 10 or more); ``thick_wall_cylinder_stresses``
is Lamé's elastic solution away from the ends, with the axial stress of a
closed cylinder. ``principal_stresses_3d`` returns the eigenvalues of the
stress tensor in decreasing order; the absolute maximum shear is
(sigma_1 - sigma_3)/2. ``strain_rosette`` takes three gauge readings at
angles counterclockwise from ``x`` (rectangular 0-45-90 by default) and
returns the engineering shear strain; ``plane_stress_from_strains`` is
Hooke's law for plane stress, valid at a free surface.

Column design
-------------

``effective_length_factor`` returns theoretical K values; codes recommend
larger ones because no connection is perfectly rigid. The secant formula
assumes a linear elastic, initially straight column with an eccentric load
and returns the peak compressive stress as a positive number; it is
undefined at and beyond the Euler load. ``secant_formula_load`` solves for
first yield by Brent's method. ``johnson_critical_stress`` is the parabola
for KL/r below Cc = sqrt(2 pi^2 E / sigma_y) and Euler's hyperbola above,
with no safety factor.

Trusses
-------

Nodes have shape (n,2); member indices start at zero. DOFs are x0,y0,x1,y1,...
and fixed DOFs have zero prescribed displacement. Members are pin-jointed,
axial-only, and use a shared area and modulus. Member force is positive in
tension. ``truss_determinacy`` returns Maxwell's count m + r - 2j; it is
necessary but not sufficient for stability. ``method_of_joints`` requires a
count of zero and solves joint by joint from statics, raising a ValueError if
every remaining joint has more than two independent unknowns.
``PlanarTruss.solve`` refuses a negative count, raises a ValueError when the
reduced stiffness K_ff is singular to working precision, and otherwise
returns cond(K_ff) in the result and issues a RuntimeWarning above
``cond_limit`` (1e12 by default): a near-mechanism gives a matrix that is
invertible in floating point but physically meaningless.

Continuous beams and frames
---------------------------

``three_moment_equation`` assumes simply supported ends, one modulus and a
uniform load per span; support moments are sagging positive. ``influence_line``
is for a simply supported span and a unit downward load. ``moment_distribution``
assumes no joint translation (no sway) and uses Hardy Cross's sign
convention: member-end moments are clockwise positive. ``PlanarFrame`` has
three DOFs per node (u, v, theta), counterclockwise rotations and moments,
Euler-Bernoulli members with axial stiffness, and uniform member loads along
the local -y axis converted to fixed-end nodal loads. Its end moments are
counterclockwise positive, the negative of Cross's. ``slope_deflection``
takes the same member tuples as ``moment_distribution`` and solves the
joint equations directly, with clockwise rotations and moments; known chord
rotations (support settlements, psi = Delta/L clockwise) can be given per
member, but sway is not solved for. The tests check the
three-moment equation, moment distribution and the frame solver against one
another, the frame against cantilever and fixed-beam formulas, and the truss
solvers against each other.

Stress, energy, stability and plasticity
----------------------------------------

Mohr circle uses (center, radius, sigma_max, sigma_min), all Pa; rotations in
``stress_transformation`` are counterclockwise. ``castigliano_deflection``
either returns a derivative dU/dP computed by the caller or integrates
M (dM/dP)/(EI) along a beam with Simpson's rule; only bending energy is
included. Maxwell-Betti compares two supplied displacements; it does not
establish linearity of an arbitrary model. ``tangent_modulus_load`` accepts a
tangent modulus Et(sigma) that is positive at zero stress and nonincreasing,
and solves sigma = pi^2 Et(sigma)/(KL/r)^2 by Brent's method.
``plastic_collapse_load`` is rigid-plastic: no strain hardening, instability,
or reduction of Mp by shear or axial force. ``unit_load_truss_deflection``
sums n N L/(A E) over the members; in an indeterminate truss the real forces
must come from a compatible solution. ``unit_load_beam_deflection`` is the
keyword form of ``castigliano_deflection`` under its virtual-work name.

Cables and arches
-----------------

Cables are perfectly flexible and inextensible, with supports at the same
level and the sag measured at midspan. ``parabolic_cable`` takes a load per
unit horizontal length (a suspended deck), ``catenary_cable`` a weight per
unit length of cable; both return the horizontal tension, the tension at the
supports and the cable length. ``three_hinged_arch`` is a parabolic arch
with pinned springings at the same level and a crown hinge at midspan,
under a load uniform along the span and point loads, downward positive. Its
moment is sagging positive (intrados in tension), its axial force tension
positive, and the thrust is reported as a positive number. The shape is not
updated under load.

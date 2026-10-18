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
counterclockwise positive, the negative of Cross's. The tests check the
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
or reduction of Mp by shear or axial force.

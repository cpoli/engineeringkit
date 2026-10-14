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
tension. Singular reduced stiffness produces a ValueError identifying an
unstable/insufficiently constrained model. Near-mechanisms can still be poorly
conditioned; a solvable matrix does not guarantee a useful physical model.
The single-bar benchmark checks displacement, force and reaction balance.

Mohr circle uses (center, radius, sigma_max, sigma_min), all Pa. Castigliano's
helper simply returns a derivative computed by the caller; it is not a
symbolic or numerical strain-energy solver. Maxwell-Betti compares two supplied
displacements; it does not establish linearity of an arbitrary model.

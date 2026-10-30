Heat transfer: absolute temperatures and resistance accounting
==============================================================

Absolute temperatures are K; differences are K. Conductivity is W/(m K),
convection coefficient W/(m^2 K), heat rate W, and resistance K/W. Radiation
uses absolute fourth powers; passing Celsius silently changes the physics.
Geometry and properties are constant, and the equations are ideal baselines.

Conduction and convection
-------------------------

Full wall resistances include area: L/(k*A), ln(r2/r1)/(2*pi*k*L), and
(1/r1-1/r2)/(4*pi*k). Cylinder/sphere radii require r2>r1>0. composite_wall
returns each temperature drop from Q*R and includes supplied film resistances.
A missing film is represented by None, not zero; zero h is rejected. At least
one physical layer is required.

The legacy overall_u name returns conductance 1/sum(R) [W/K] for full
resistances. It does not return area-normalized U [W/(m^2 K)]. To obtain U,
divide the conductance by area. The documentation now makes this distinction
explicit without changing existing numeric behavior.

fin_heat_rate_convective_tip adds the tip face, with h_tip = h unless
given; h_tip = 0 is the adiabatic fin, and L = 0 the bare base area.

Lumped cooling requires Bi=h*Lc/k typically below 0.1, where Lc=V/A. The
function cannot check Bi without conductivity; that assessment remains the
caller's responsibility. Semi-infinite conduction requires alpha>0 and t>0;
t=0 is excluded rather than resolving the discontinuous boundary condition.

Exchangers and limits
---------------------

NTU is nonnegative and Cr=Cmin/Cmax lies in [0,1]. Counterflow at Cr=1 uses
NTU/(1+NTU). Stable expm1 expressions preserve tiny NTU, and log1p preserves
near-equal terminal temperature differences in LMTD. Temperature crossing is
excluded by requiring both terminal differences positive.

exchanger_effectiveness_crossflow takes mixed="none" (both streams
unmixed; the usual approximate formula, within 4% of the exact series for
NTU <= 5 and about 1% for 2 <= NTU <= 5), "cmax" or "cmin" (exact). Every
arrangement gives 1 - exp(-NTU) at Cr = 0. exchanger_effectiveness_shell_tube
takes the total NTU and the number of shell passes, each with an even number
of tube passes; many shells approach counterflow.

lmtd_correction_factor(P, R, shells) uses the tube-side effectiveness
P = (t_o - t_i)/(T_i - t_i) and R = (T_i - T_o)/(t_o - t_i) = C_tube/C_shell,
and multiplies the *counterflow* LMTD. It is symmetric, F(P, R) = F(PR, 1/R).
A P that the arrangement cannot reach for that R (where F would fall to
zero, a temperature cross) raises ValueError instead of returning NaN; add
shells. The F-LMTD and effectiveness-NTU methods agree to round-off, which
the tests check.

Transient conduction and numerics
---------------------------------

plane_wall_transient takes dimensionless inputs: Fo = alpha*t/L^2 and
Bi = h*L/k on the half-thickness L, and x/L measured from the midplane.
The default 100 series terms are accurate for Fo above about 1e-3; near
Fo = 0 the truncated series oscillates, so pass more terms. terms=1 is
Heisler's one-term chart, within 2% at the midplane for Fo > 0.2. Bi is a
scalar; the eigenvalues of zeta*tan(zeta) = Bi are bracketed exactly
(utils.plane_wall_eigenvalues).

conduction_1d uses evenly spaced nodes, fixed surface temperatures or
adiabatic faces (None), and constant properties. The explicit scheme is
stable only for r = alpha*dt/dx^2 <= 1/2; it is not refused above that so
the instability can be studied, and result.mesh_fourier reports r.
Crank-Nicolson is unconditionally stable but rings for a few steps after a
discontinuous initial or boundary condition when r is large.

Convection correlations
-----------------------

flat_plate_nusselt is the isothermal flat plate in parallel flow: laminar
0.332 Re^(1/2) Pr^(1/3) locally and 0.664 on average (Pohlhausen),
turbulent 0.0296 Re^(4/5) Pr^(1/3), and the mixed average with transition at
Re_c = 5e5. Re_c = 0 makes it turbulent from the leading edge and Re_c = inf
keeps it laminar. The Pr^(1/3) laws are meant for 0.6 < Pr < 60.
churchill_bernstein_nusselt covers a cylinder in crossflow for Re Pr > 0.2,
properties at the film temperature. zukauskas_nusselt is the average over a
tube bank: Re is based on tube_bank_max_velocity, properties on the mean
fluid temperature except Pr_s, Re must lie in [10, 2e6] (enforced), and banks
of fewer than 20 rows are corrected with C2 interpolated from Zukauskas's
table. In the range 100-1000 it uses Zukauskas's constants (0.52 and 0.71,
m = 0.5) where Incropera et al. suggest the single-cylinder correlation.

laminar_tube_nusselt returns 48/11 (uniform flux) or 3.657 (uniform wall
temperature, Graetz's limit) for fully developed laminar flow in a circular
tube. entry_length gives 0.05 Re Pr D below Re = 2300 (Pr = 1 for the
hydrodynamic length) and 10 D above it. hausen_nusselt is the mean over a
thermal entry length with a developed velocity profile and uniform wall
temperature; it tends to 3.66.

chilton_colburn_mass_transfer gives h_m = h/(rho cp Le^(2/3)) [m/s], with
lewis_number = alpha/D_AB. It is exact between correlations with Pr^(1/3)
and Sc^(1/3) (the laminar flat plate) and approximate for turbulent flow;
low mass-transfer rates are assumed.

dittus_boelter_nusselt is for fully developed turbulent tube flow,
Re above about 1e4 and 0.6 < Pr < 160, with typical error 25%; it is not
checked. churchill_chu_nusselt is for a vertical plate with properties at
the film temperature; laminar=True selects the more accurate laminar form
for Ra below about 1e9. rayleigh_number takes the magnitude of the
temperature difference.

Radiation and view factors
--------------------------

All radiation temperatures are absolute. planck_spectral_emissive_power
returns W/(m^2 m) per metre of wavelength (multiply by 1e-6 for per
micrometre) and is written to underflow, not overflow, at short wavelength.
radiation_heat_rate treats a small gray surface in large isothermal
surroundings. The view factors assume diffuse surfaces: aligned, directly
opposed rectangles and coaxial parallel disks.

radiation_enclosure solves Oppenheim's network for n diffuse gray surfaces:
T, epsilon and A have one entry per surface and F[i, j] is F_ij. The rows of F
must sum to one and A_i F_ij = A_j F_ji must hold (both to 1e-6), otherwise
ValueError: a missing surface (an opening) is a surface at the surroundings'
temperature with epsilon = 1. Surfaces flagged reradiating are insulated;
their T is ignored and solved for. The net heat rates sum to zero.
two_surface_radiation is the two-surface closed form (emissivities in
(0, 1]); radiation_shield_heat_flux the large parallel plates with thin
shields, given one emissivity or a pair per shield. view_factor_reciprocity
returns F_ji and refuses a result above one.

Two-dimensional conduction and shape factors
--------------------------------------------

conduction_2d solves steady conduction in a rectangle with the four edge
temperatures fixed, constant k and no generation, on an nx by ny grid that
includes the edges (corners take the mean of their two edges and enter no
equation). heat_rates gives W per metre of depth entering through each edge;
they sum to zero to round-off. The error falls as dx^2. The shape factors
return S [m] with Q = S k dT: buried sphere and cylinder (isothermal ground
surface, cylinder length much larger than D), eccentric cylinders and two
parallel cylinders. Geometries that overlap raise ValueError.

Transient cylinders and spheres
-------------------------------

cylinder_transient and sphere_transient take Fo = alpha t/r_o^2 and
Bi = h r_o/k on the *outer radius* (not V/A) and r/r_o. Their eigenvalues
(utils.cylinder_eigenvalues, utils.sphere_eigenvalues) are bracketed
between known zeros. For small Bi they tend to the lumped exp(-2 Bi Fo) and
exp(-3 Bi Fo); terms=1 is Heisler's one-term chart, good for Fo > 0.2.

Phase change
------------

neumann_stefan_front is the one-phase Stefan problem: the original phase
stays at the melting point and equal densities are assumed. Ste =
cp*dT/latent must be positive. nusselt_film_condensation is laminar film
condensation on a vertical wall (film Reynolds number below about 1800;
ripples above about 30 raise h by 10-20%); pass Rohsenow's
h_fg + 0.68*cp*(Tsat - Ts) to include film subcooling.
rohsenow_boiling_heat_flux is saturated nucleate pool boiling, with C_sf =
0.013 (water on polished copper) and n = 1 (water; 1.7 for other liquids)
by default. It grows as dT_excess^3 and is valid only up to the critical
heat flux, which it does not check; errors of 100% in flux (about 25% in
superheat) are common. zuber_critical_heat_flux uses Zuber's C = pi/24 by default; 0.149 fits
large flat heaters. Both require rho_l > rho_v.

.. code-block:: python

    from engineeringkit.heat_transfer import composite_wall
    result = composite_wall(400, 300, [(.1, 1), (.2, 2)], A=2)
    # Q=1000 W, interfaces [400, 350, 300] K.
    fig, ax = result.plot()

Tests check interface energy balance, cylindrical thin-wall limit, equal-
temperature radiation, cooling endpoints and exchanger limiting cases, the
plane-wall eigenvalue equation and fixed-surface limit, energy conservation
of insulated finite-difference conduction, Planck's integral, view-factor
limits and the Neumann front equation. The curriculum models are checked
against the Blasius-Pohlhausen similarity solution, Graetz's eigenvalue,
the exact crossflow series, the effectiveness-NTU/F-LMTD identity, the
cylinder and sphere lumped, fixed-surface and energy limits, the 2-D Fourier
series and its second-order convergence, shape-factor identities and
the closed-form two- and three-surface networks. Fin
formulas assume a straight constant-section fin with an adiabatic tip.

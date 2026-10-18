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

Phase change
------------

neumann_stefan_front is the one-phase Stefan problem: the original phase
stays at the melting point and equal densities are assumed. Ste =
cp*dT/latent must be positive. nusselt_film_condensation is laminar film
condensation on a vertical wall (film Reynolds number below about 1800;
ripples above about 30 raise h by 10-20%); pass Rohsenow's
h_fg + 0.68*cp*(Tsat - Ts) to include film subcooling.
zuber_critical_heat_flux uses Zuber's C = pi/24 by default; 0.149 fits
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
limits and the Neumann front equation. Fin
formulas assume a straight constant-section fin with an adiabatic tip.

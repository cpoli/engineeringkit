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

.. code-block:: python

    from engineeringkit.heat_transfer import composite_wall
    result = composite_wall(400, 300, [(.1, 1), (.2, 2)], A=2)
    # Q=1000 W, interfaces [400, 350, 300] K.
    fig, ax = result.plot()

Tests check interface energy balance, cylindrical thin-wall limit, equal-
temperature radiation, cooling endpoints and exchanger limiting cases. Fin
formulas assume a straight constant-section fin with an adiabatic tip.

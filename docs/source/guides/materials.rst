Materials: structure, elasticity, yield, fatigue, fracture and creep conventions
================================================================================

Stresses and moduli are Pa, strains dimensionless, lengths m, stress
intensity Pa sqrt(m), energies per volume J/m^3, temperatures K, angles
rad, diffusivities m^2/s and molar masses kg/mol. Hardness numbers are the
exception: they keep their conventional kgf/mm^2 values.
Engineering (nominal) stress and strain and small-strain elasticity are
assumed throughout. Every function accepts NumPy arrays and broadcasts its
arguments, except where a vector is stated.

Crystals
--------

cubic_lattice_parameter gives the cell edge of touching hard spheres
(sc 2R, bcc 4R/sqrt(3), fcc 2*sqrt(2)*R); atomic_packing_factor also
accepts "hcp" with the ideal c/a = sqrt(8/3). theoretical_density is
n*A/(Vc*N_A) for a perfect crystal; vacancies lower the measured value.
interplanar_spacing is for cubic cells only and needs integer Miller
indices. bragg_angle returns theta, not 2*theta, and raises when
n*lambda > 2d. diffraction_peaks lists first-order peaks with the BCC
(h+k+l even) and FCC (all odd or all even) selection rules; planes with the
same h^2+k^2+l^2 share a peak and are listed once. Peak intensities are
not modelled.

Diffusion and phases
--------------------

diffusion_coefficient is D0*exp(-Qd/RT) with R = 8.314 J/(mol K).
fick_flux is the steady first law. semi_infinite_concentration and
diffusion_time are the error-function solution for a constant D, a
uniform initial concentration and a surface held at Cs from t = 0; the
part must be thicker than a few sqrt(D*t). Concentrations may be in any
unit, as long as C, C0 and Cs share it.

lever_rule takes the tie-line compositions read from a phase diagram (no
diagrams are built in) and returns (W_alpha, W_beta); weight percent gives
mass fractions. avrami_fraction is 1 - exp(-k*t^n) for an isothermal
transformation, with k in 1/s^n; avrami_time inverts it and avrami_fit
fits k and n on Avrami axes, so trim points near 0 and 1 first.

.. code-block:: python

    from engineeringkit.materials import diffusion_coefficient, diffusion_time
    D = diffusion_coefficient(2.3e-5, 148e3, 1223.0)    # C in austenite
    t = diffusion_time(0.5e-3, 0.80, D, 0.25, 1.20)     # s to 0.8 wt% at 0.5 mm

Elasticity
----------

For an isotropic stable elastic material E>0 and -1<nu<0.5. The
incompressible limit nu=0.5 makes the bulk modulus and Lamé's lambda
infinite and is excluded. lame_parameters returns (lambda, mu) with mu the
shear modulus. constrained_thermal_stress is the uniaxial restraint
-E*alpha*dT; a plate restrained in two directions carries
-E*alpha*dT/(1-nu), which is not computed.

StressStrainCurve requires matching finite vectors with increasing strain.
toughness integrates stress d(strain) by the trapezoid rule; young_modulus
fits the first n samples with 2<=n<=sample count, and the caller must
choose an elastic segment. offset_yield_strength finds the first crossing
of the curve with the elastic line shifted by the offset strain (0.002 by
default) by linear interpolation, and raises if the curve never reaches it.

Hardness
--------

brinell_hardness and vickers_hardness take the load in N and lengths in m
and return the conventional HB and HV numbers (kgf/mm^2). Brinell is
valid for dents between 0.24 and 0.6 ball diameters. brinell_tensile_strength
is the empirical 3.45 MPa per HB for steels, good to about 10%; pass
another factor for other alloys.

Yield
-----

tresca_equivalent takes the principal stresses along the last axis and
returns their spread; include the zero out-of-plane stress for plane
stress, or equal biaxial tension will be wrongly reported as safe.
von_mises_plane_stress assumes sigma_z = tau_xz = tau_yz = 0. Both are for
isotropic, pressure-insensitive ductile metals; Tresca's equivalent is never
smaller than von Mises's and at most 2/sqrt(3) times it.
principal_stresses takes the six components of a three-dimensional state
(shears default to zero) and returns sigma_1 >= sigma_2 >= sigma_3 along
the last axis, ready for tresca_equivalent; von_mises_stress is the
three-dimensional equivalent and equals von_mises_plane_stress when
sigma_z = tau_yz = tau_zx = 0.

Composites
----------

The Voigt (longitudinal) and Reuss (transverse) rules of mixtures require
Vf in [0,1]. They are the exact moduli of continuous aligned fibres loaded
along and across them, and bounds for any other arrangement; real
transverse moduli lie above the Reuss value. Voids and fibre waviness are
not modelled. halpin_tsai interpolates between them with xi >= 0: xi = 0 is
Reuss, xi -> infinity is Voigt, xi = 2 for E2 and 1 for G12 of round fibres,
2*l/d for E1 of aligned short fibres.

classical_laminate builds the ABD matrix from the ply properties (E1, E2,
nu12, G12), the fibre angles in radians listed from the bottom ply up, and
one thickness or one per ply. The returned Laminate solves for midplane
strains and curvatures, gives each ply's stresses in its material axes at
mid-height, and the in-plane engineering constants from A. Kirchhoff
kinematics and perfect bonding are assumed; free-edge interlaminar
stresses, hygrothermal strains and failure criteria (Tsai-Wu, Hashin) are
not modelled.

Fatigue
-------

Basquin uses sigma_a = sigma_f*(2N)^b with b<0, Coffin-Manson
eps_pa = eps_f*(2N)^c with c<0, and both return N in cycles, not reversals.
strain_life_amplitude returns (total, elastic, plastic) amplitudes; no mean
stress correction (Morrow, Smith-Watson-Topper) is applied. Goodman
restricts mean stress to nonnegative values. Miner sums n_i/N_i over
matching vectors and ignores load order; it is not a probabilistic fatigue
predictor. No endurance limit is built into the S-N law.

gerber_safety_factor solves n*sa/Se + (n*sm/Sut)^2 = 1 in closed form and
soderberg_safety_factor uses the yield strength in place of Sut; for the
same stresses Soderberg <= Goodman <= Gerber. endurance_limit_estimate is
0.5*Sut capped at 700 MPa, for steels only. The Marin factors follow
Shigley: marin_surface_factor (ground, machined or cold-drawn, hot-rolled,
forged), marin_size_factor for round bars in bending or torsion between
2.79 and 254 mm (use 1 for axial load), marin_load_factor (1, 0.85, 0.59)
and marin_reliability_factor 1 - 0.08*z for reliability in [0.5, 1).
marin_endurance_limit multiplies them; the temperature factor kd is left
to the caller.

Notches
-------

elliptical_hole_concentration is Inglis's 1 + 2a/b for a wide plate.
peterson_notch_sensitivity and neuber_notch_sensitivity take the root
radius and a material length in m, and fatigue_notch_factor is
1 + q*(Kt - 1). neuber_notch_stress solves sigma*eps = (Kt*S)^2/E with the
Ramberg-Osgood curve of ramberg_osgood_strain; for fatigue use the cyclic
curve and Kf, and pass a stress range with the doubled curve for the
Massing hysteresis branch. Nominal stresses should be elastic.

Fracture
--------

K = Y*sigma*sqrt(pi*a) assumes tensile mode-I loading with a supplied
geometry factor. Critical length inverts it. These linear-elastic models do
not check small-scale yielding: the plastic zone (K/sigma_y)^2/(2*pi) must
be small compared with the crack and the thickness. Griffith uses surface
energy J/m^2, a central crack of half-length a and plane stress.

Paris uses da/dN = C*(delta_K)^m; C must be calibrated to Pa sqrt(m), not
MPa sqrt(m) (multiply by 1e-6**m), and only the stable-growth regime is
modelled. crack_growth_life integrates the law by the trapezoid rule on a
geometric grid of crack lengths; Y may be a function of crack length.

The Weibull strength functions assume a uniformly stressed volume:
weibull_failure_probability is 1 - exp(-(V/V0)*(sigma/sigma_0)^m),
evaluated by engineeringkit.reliability.weibull_cdf with the scale
sigma_0*(V0/V)^(1/m). For bending or other non-uniform stress, replace V
by the effective volume. weibull_size_effect scales a strength between
volumes at equal failure probability.

.. code-block:: python

    from engineeringkit.materials import crack_growth_life, critical_crack_length
    a_c = critical_crack_length(35e6, 100e6)
    life = crack_growth_life(5e-32, 3.5, 100e6, 1e-3, a_c).life   # cycles

Creep and plasticity
--------------------

arrhenius_creep_rate is the steady (secondary) rate only; primary and
tertiary creep are not modelled. larson_miller_parameter takes the rupture
time in seconds and converts it to hours internally, so that the usual
C = 20 applies; the parameter is returned in K. ideal_shear_strength is
Frenkel's sinusoidal estimate G*b/(2*pi*a); peierls_stress is the
Peierls-Nabarro estimate, a guide to trends rather than an accurate value.
hall_petch holds down to grain sizes of tens of nanometres.
ramberg_osgood_strain uses the offset form eps = sigma/E +
offset*(sigma/sigma_y)^n with n >= 1 and is odd in stress.

Tests check the textbook values of copper's density, the (220) reflection
of iron, Fick's flux, a carburizing time, the lever rule, a Brinell number
and Shigley's Marin factors; the error-function profile against Fick's
second law, Avrami fit round trips, Gerber-Goodman-Soderberg ordering,
Neuber's product, principal-stress invariants, Halpin-Tsai's limits at the
Reuss and Voigt bounds, the ABD matrix of single, quasi-isotropic and
cross-ply laminates, Weibull's size effect against the reliability
distribution, and also the isotropic modulus identities, Lamé's law in uniaxial
tension, elastic triangle energy, fitted modulus and offset yield, the
Tresca and von Mises ratio in shear, fracture inversion and Griffith-Irwin
agreement, Paris integration against the closed form and quadrature,
composite bounds, the strain-life transition, the Larson-Miller round trip
and Miner accumulation.

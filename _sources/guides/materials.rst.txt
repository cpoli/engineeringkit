Materials: elasticity, yield, fatigue, fracture and creep conventions
=====================================================================

Stresses and moduli are Pa, strains dimensionless, lengths m, stress
intensity Pa sqrt(m), energies per volume J/m^3 and temperatures K.
Engineering (nominal) stress and strain and small-strain elasticity are
assumed throughout. Every function accepts NumPy arrays and broadcasts its
arguments, except where a vector is stated.

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

Yield
-----

tresca_equivalent takes the principal stresses along the last axis and
returns their spread; include the zero out-of-plane stress for plane
stress, or equal biaxial tension will be wrongly reported as safe.
von_mises_plane_stress assumes sigma_z = tau_xz = tau_yz = 0. Both are for
isotropic, pressure-insensitive ductile metals; Tresca's equivalent is never
smaller than von Mises's and at most 2/sqrt(3) times it.

Composites
----------

The Voigt (longitudinal) and Reuss (transverse) rules of mixtures require
Vf in [0,1]. They are the exact moduli of continuous aligned fibres loaded
along and across them, and bounds for any other arrangement; real
transverse moduli lie above the Reuss value. Voids and fibre waviness are
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

Tests check the isotropic modulus identities, Lamé's law in uniaxial
tension, elastic triangle energy, fitted modulus and offset yield, the
Tresca and von Mises ratio in shear, fracture inversion and Griffith-Irwin
agreement, Paris integration against the closed form and quadrature,
composite bounds, the strain-life transition, the Larson-Miller round trip
and Miner accumulation.

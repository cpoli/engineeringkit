Materials: elasticity, fracture and fatigue conventions
=======================================================

Stress/moduli are Pa, strain dimensionless, crack length m and stress intensity
Pa sqrt(m). Engineering strains and small-strain elasticity are assumed.
For an isotropic stable elastic material E>0 and -1<nu<0.5. The incompressible
limit nu=0.5 makes bulk modulus infinite and is excluded.

StressStrainCurve requires matching finite vectors with increasing strain.
Toughness integrates stress d(strain), producing J/m^3; young_modulus fits the
first n samples with 2<=n<=sample count. It is the caller's responsibility to
select a genuinely elastic segment rather than fitting yield/plastic response.
SciPy trapezoid is used to retain compatibility with the stated NumPy minimum.

Fracture and fatigue
--------------------

K=Y*sigma*sqrt(pi*a) assumes tensile mode-I loading with a supplied dimensionless
geometry factor. Critical length inverts this relation. These linear-elastic
models do not check small-scale yielding or supply crack geometry corrections.
Griffith uses surface energy J/m^2 and an ideal brittle central crack convention;
real fracture toughness depends on constraint and plasticity.

Basquin uses sigma_a=sigma_f*(2*N)^b with b<0 and returns N in cycles, not
reversals. Paris uses da/dN=C*(delta_K)^m; C must be calibrated to Pa sqrt(m),
not MPa sqrt(m), and the stable-growth regime is the caller's responsibility.
Miner sums n_i/N_i and assumes sequence-independent damage; it is not a
probabilistic fatigue predictor. Goodman restricts mean stress to nonnegative
values; compressive-mean-stress corrections require another model.

.. code-block:: python

    from engineeringkit.materials import critical_crack_length, stress_intensity
    a = critical_crack_length(50e6, 100e6, 1.12)
    assert abs(stress_intensity(1.12, 100e6, a) / 50e6 - 1) < 1e-12

Tests check the isotropic modulus identity, elastic triangle energy, fitted
modulus, fracture inversion, composite Voigt/Reuss bounds, yield invariants
and Miner accumulation. Rule-of-mixtures inputs require Vf in [0,1].

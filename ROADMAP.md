# Roadmap

What engineeringkit does not do yet, in priority order, and the plan of
every domain's history page. engineeringkit is a teaching toolkit, so the
measure is not feature parity with python-control, CoolProp or a finite
element code. It is whether each domain covers the canonical topics of an
undergraduate course in that field, tells their history, and checks that
every example still shows the right engineering.

All twelve domains exist and are tested (0.1.0, 465 tests, 100% statement
coverage), but only two have been rebuilt to the current standard:
`aerodynamics` and `communications` have the `core/systems/utils/visualizers`
layout, a dated history page of 15+ breakthroughs and one gallery example
per breakthrough. `control_systems` has the layout and 13 examples, but its
history page has not been rewritten. The other nine still have three undated
history sections and two generic examples each. Rebuilding them is the
phase checklist below. Items 1 to 4 are the work that is not tied to a
single domain.

The project conventions apply to everything here: SI inputs validated in
`_validation.py`, NumPy-style docstrings with units and the guide Notes
block, an analytical or benchmark test for every new model, 100% production
statement coverage, `filterwarnings = ["error"]`, and, for each breakthrough,
a history entry with its own gallery example, titled after it, that ends in
a check of the result it illustrates and is linked from no other entry.

## Phase checklist

Each domain is rebuilt in its own session, from this file. A rebuild means
the `core/ systems/ utils/ visualizers/ tests/` layout (as in
`aerodynamics`), the public API kept backward compatible through the
domain's `__init__.py`, a history page in the dated format of
`communications_breakthroughs.rst`, and the gallery examples named below.
Tick the box when ruff, mypy, pytest, coverage and the `-W` docs build are green.

- [x] 1. Twelve domains, validation, benchmarks, guides, 100% coverage (0.1.0)
- [ ] 2. Domain rebuilds, one per session
  - [x] `aerodynamics` (21 breakthroughs, 21 examples)
  - [x] `communications` (17 breakthroughs, 17 examples)
  - [ ] `control_systems`: layout and 13 examples done; history page and
        the remaining examples to write
  - [ ] `electrical_engineering`
  - [ ] `fluid_engineering`
  - [ ] `heat_transfer`
  - [ ] `materials`
  - [ ] `mechanical_engineering`
  - [ ] `reliability`
  - [ ] `signal_processing`
  - [ ] `structural_engineering`
  - [ ] `thermodynamics`
- [ ] 3. Cleanup once every domain is rebuilt
  - [ ] Delete the per-domain `plot_<domain>_models.py` and
        `plot_<domain>_benchmark.py` examples, `examples/all_domains.py`,
        `examples/phase1_gallery.py`, `gallery/integrators.png` and
        `engineeringkit/tests/test_phase1.py` if nothing else uses them
  - [ ] Move the remaining sibling `test_*.py` files into `tests/` packages
  - [ ] Update the Documentation and domain-layout paragraphs of CLAUDE.md
        and AGENTS.md (they still describe single-module domains and
        `plot_<domain>_models.py`)
  - [ ] Rewrite `docs/source/limitations.rst` around what is still open
- [ ] 4. Course material (item 4)

## 1. Release infrastructure

**Why.** engineeringkit cannot be cited or installed from PyPI, and the
hosted docs do not exist. The sibling kits all have this already.

**Now.** The directory is not a git repository. There is only `ci.yml`: no
`CITATION.cff`, no release or docs workflow. The CHANGELOG says no external
publication has been made.

**Approach.** Port from physicskit, renaming as needed:
- `git init` with author cpoli374 <cpoli374@gmail.com>, then create the
  GitHub repository `cpoli/engineeringkit`.
- `CITATION.cff`, plus `release.yml`: on a `vX.Y.Z` tag, check the tag
  against `engineeringkit.__version__` and `CITATION.cff`, build, publish
  to PyPI with trusted publishing, create the GitHub Release from the
  CHANGELOG section and rebuild `gh-pages`.
- `docs.yml`, which deploys `gh-pages` on every push to `main`. Then
  enable Pages from `gh-pages`.
- CI matrix to 3.10-3.14, as in physicskit.
- After the first release, enable Zenodo and put that release's own
  **version DOI** (never the concept DOI) in `CITATION.cff` and the README
  badge (`https://zenodo.org/badge/DOI/<version-doi>.svg`).
- A conda-forge `recipe.yaml` later, with only `version` in `context` and
  `name: engineeringkit` and the PyPI source URL written out.

## 2. A blocking type check

**Done** (2026-10-13). `python -m mypy` checks the whole package with
`check_untyped_defs`, `strict_equality` and `warn_unreachable` (tests at
the signature level only), as physicskit does, and runs as its own
blocking CI job. The package was already clean, so no code changed.
Rebuilt domains must keep it clean.

## 3. Numerical gaps from `limitations.rst`

**Why.** These are known places where a model can return a wrong or
misleading number. Each fits the domain rebuild named in brackets, and
should be done there unless it is needed sooner.

- **High-order IIR filters** [`signal_processing`]. `butterworth_lowpass`
  and `butterworth_highpass` return `(b, a)` polynomials, which lose
  precision above about order 8. Add an `output="sos"` option and make
  `filter_signal` accept second-order sections (the Jackson 1970 entry).
- **Stability margins with several crossovers** [`control_systems`].
  `stability_margins` reads the *first* gain and phase crossovers. Return
  the worst margin over all crossovers, and report which crossover it is,
  so conditionally stable loops come out right.
- **Truss conditioning** [`structural_engineering`]. `PlanarTruss.solve`
  catches `LinAlgError` only for exactly singular stiffness. Report
  `cond(K_ff)` in the result, warn above a threshold, and check Maxwell's
  count `m + r = 2j` before solving (the Maxwell 1864 entry).
- **Real fluids** [`thermodynamics`]. Every cycle uses ideal-gas or
  cold-air-standard relations. Add van der Waals, saturation via
  Clausius-Clapeyron and Antoine, and a small IAPWS-IF97 subset (regions
  1, 2 and 4), so the Rankine cycle can be computed with steam properties.
  CoolProp stays out of scope (no new dependencies).
- **Docstrings and vectorization**. Result classes and the
  multidimensional helpers do not all have full NumPy-style parameter
  docstrings. Scalar parameters such as `gamma` or Weibull `beta` are not
  vectorized the same way everywhere. Fix both in each domain as it is
  rebuilt.
- **Castigliano**. `castigliano_deflection` takes a strain-energy
  derivative the caller has already computed. The Castigliano entry
  should integrate `M dM/dP / EI` along a beam itself.

## 4. Course material

**Why.** The galleries show one idea at a time. A course also needs
worked problems that span several functions, and the other kits have
tutorials and a quickstart notebook.

**Approach.** After the rebuilds:
- Four or five tutorials in `docs/source/tutorials/`, each using several
  domains. Candidates: *a beam from load to fatigue life*
  (`structural_engineering`, `materials`, `reliability`); *a feedback loop
  in the time and frequency domains* (`control_systems`,
  `signal_processing`); *a power plant from cycle to cooling tower*
  (`thermodynamics`, `heat_transfer`, `fluid_engineering`); *a radio link
  from antenna to bit errors* (`communications`, `signal_processing`).
- `notebooks/quickstart.ipynb` with a Colab badge, served in JupyterLite.
- An exercise set per domain, built on the gallery examples, with
  solutions in `sphinx-design` dropdowns that run in the docs build.

## Scope decisions

- **Out of scope**, since the dependencies stay NumPy, SciPy and
  Matplotlib: CAD and meshing, general finite element or CFD solvers,
  wrappers around CoolProp, python-control or SPICE, and unit-carrying
  types (inputs stay raw SI numbers).
- **Candidate domains**, not scheduled, each needing 15+ breakthroughs:
  `geotechnical` (Terzaghi's effective stress and consolidation, Mohr-
  Coulomb, Rankine earth pressure, Casagrande), `power_systems` (load
  flow, symmetrical components and fault analysis, swing equation;
  otherwise these go in `electrical_engineering`), `manufacturing` (Taylor
  tool life, Merchant's cutting force, statistical process control).

## Where engineeringkit stands

465 tests, 100% production statement coverage. The other numbers are
from the history pages and `examples/`.

| Domain | Breakthroughs | Gallery examples | Target |
|---|---:|---:|---:|
| `aerodynamics` | 21 | 21 | done |
| `communications` | 17 | 17 | done |
| `control_systems` | 3 (old format) | 13 | 20 |
| `electrical_engineering` | 3 (old format) | 2 | 16 |
| `fluid_engineering` | 3 (old format) | 2 | 17 |
| `heat_transfer` | 3 (old format) | 2 | 17 |
| `materials` | 3 (old format) | 2 | 17 |
| `mechanical_engineering` | 3 (old format) | 2 | 16 |
| `reliability` | 3 (old format) | 2 | 17 |
| `signal_processing` | 3 (old format) | 2 | 17 |
| `structural_engineering` | 3 (old format) | 2 | 17 |
| `thermodynamics` | 3 (old format) | 2 | 18 |

## Planned breakthroughs

Each list is in chronological order, as it will appear on the history
page. It gives the function that carries the entry and the gallery example
in brackets as `<topic>/plot_NN_<slug>.py`. **new** marks a function that
does not exist yet. A list may shrink to 15 if an entry turns out to be
weak, and must not fall below 15. Dates and attributions are to be checked
against the primary reference while the entry is written.

### control_systems

Entries 1-13 already have examples; their history entries are to be written.

1. 1788 -- Watt's centrifugal governor: `feedback`, `step_response` [`feedback/plot_01_watt_centrifugal_governor.py`]
2. 1868 -- Maxwell's *On Governors*: `poles` [`stability/plot_01_maxwell_on_governors.py`]
3. 1877 -- Routh's stability criterion: `routh_first_column`, `routh_stable` [`stability/plot_02_routh_stability_criterion.py`]
4. 1892 -- Lyapunov's direct method: `lyapunov_equation`, `lyapunov_stable` [`stability/plot_03_lyapunov_direct_method.py`]
5. 1895 -- Hurwitz's determinants: `hurwitz_determinants` [`stability/plot_04_hurwitz_determinants.py`]
6. 1922 -- Minorsky's PID ship steering: `PID` [`feedback/plot_02_minorsky_pid_ship_steering.py`]
7. 1927 -- Black's negative-feedback amplifier: `feedback`, `dc_gain` [`feedback/plot_03_black_negative_feedback_amplifier.py`]
8. 1932 -- Nyquist's stability criterion: `nyquist`, `nyquist_encirclements` [`frequency/plot_01_nyquist_stability_criterion.py`]
9. 1934 -- Hazen's theory of servomechanisms: `second_order_metrics` [`time_response/plot_01_hazen_servomechanism.py`]
10. 1940 -- Bode's gain and phase margins: `bode`, `stability_margins` (all crossovers, item 3) [`frequency/plot_02_bode_stability_margins.py`]
11. 1942 -- Ziegler-Nichols tuning: `ziegler_nichols` [`feedback/plot_04_ziegler_nichols_tuning.py`]
12. 1948 -- Evans's root locus: `root_locus` [`root_locus/plot_01_evans_root_locus.py`]
13. 1952 -- Ragazzini and Zadeh: sampled-data systems and the zero-order hold: `discretize_zoh` [`sampled/plot_01_ragazzini_zadeh_zero_order_hold.py`]
14. 1957 -- Smith's predictor for dead time: **new** `pade_delay`, `smith_predictor` [`time_response/plot_02_smith_predictor.py`]
15. 1960 -- Kalman's controllability and observability: `is_controllable`, `is_observable` [`state_space/plot_01_kalman_controllability_observability.py`]
16. 1960 -- Kalman's linear-quadratic regulator: `lqr`, `dlqr` [`optimal/plot_01_kalman_lqr.py`]
17. 1961 -- Jury's stability test for sampled systems: **new** `jury_stable` [`sampled/plot_02_jury_stability_test.py`]
18. 1961 -- The Kalman-Bucy filter: `kalman_gain` [`optimal/plot_02_kalman_bucy_filter.py`]
19. 1964 -- Luenberger's observer: `observer_gain` [`state_space/plot_02_luenberger_observer.py`]
20. 1967 -- Wonham's pole assignment: `state_feedback_poles`, `feedforward_gain` [`state_space/plot_03_wonham_pole_assignment.py`]
21. 1978 -- Doyle: guaranteed margins for LQG, "there are none": `lqr`, `kalman_gain`, `stability_margins` [`optimal/plot_03_doyle_lqg_margins.py`]

### electrical_engineering

1. 1827 -- Ohm's law: `ohms_law`, `series_resistance`, `parallel_resistance` [`dc_circuits/plot_01_ohm_law.py`]
2. 1843 -- Wheatstone's bridge: `wheatstone_balanced` [`dc_circuits/plot_02_wheatstone_bridge.py`]
3. 1845 -- Kirchhoff's circuit laws and nodal analysis: **new** `nodal_analysis` [`dc_circuits/plot_03_kirchhoff_nodal_analysis.py`]
4. c. 1840 -- Jacobi's maximum power theorem: `maximum_power_transfer` [`networks/plot_01_jacobi_maximum_power.py`]
5. 1855 -- Kelvin's telegraph cable and the RC time constant: `rc_time_constant`, `rc_charging_voltage` [`transients/plot_01_kelvin_rc_cable.py`]
6. 1876-1887 -- Heaviside's telegrapher's equations and the distortionless line: `transmission_line_z0`, `transmission_line_velocity`, **new** `distortionless_condition` [`transmission_lines/plot_01_heaviside_telegrapher.py`]
7. 1883 -- Thévenin's equivalent circuit: `thevenin_from_source` [`networks/plot_02_thevenin_equivalent.py`]
8. 1886 -- Stanley's transformer: `transformer_secondary` [`ac_circuits/plot_01_stanley_transformer.py`]
9. 1888 -- Tesla's polyphase system: `three_phase_power` [`power/plot_01_tesla_three_phase.py`]
10. 1893 -- Steinmetz's phasors and complex impedance: `capacitor_impedance`, `inductor_impedance`, `series_rlc_impedance` [`ac_circuits/plot_02_steinmetz_phasors.py`]
11. 1897 -- Lodge's syntonic tuning: series resonance: `resonance_frequency`, `SeriesRLC` [`ac_circuits/plot_03_lodge_resonance.py`]
12. 1918 -- Fortescue's symmetrical components: **new** `symmetrical_components` [`power/plot_02_fortescue_symmetrical_components.py`]
13. 1926 -- Norton's equivalent circuit: `norton_from_thevenin` [`networks/plot_03_norton_equivalent.py`]
14. 1927 -- Budeanu's reactive power and the power factor: `real_power`, `reactive_power`, `power_factor` [`power/plot_03_budeanu_reactive_power.py`]
15. 1939 -- Smith's chart: **new** `reflection_coefficient`, `input_impedance` [`transmission_lines/plot_02_smith_chart.py`]
16. 1956 -- Ward and Hale's digital load flow: **new** `gauss_seidel_load_flow` [`power/plot_04_ward_hale_load_flow.py`]

### fluid_engineering

1. c. 250 BC -- Archimedes's principle: `buoyancy_force` [`hydrostatics/plot_01_archimedes_buoyancy.py`]
2. 1586 -- Stevin's hydrostatic paradox: `hydrostatic_pressure`, `pressure_head` [`hydrostatics/plot_02_stevin_hydrostatic_paradox.py`]
3. 1732 -- Pitot's tube: `pitot_velocity`, `dynamic_pressure` [`measurement/plot_01_pitot_tube.py`]
4. 1797 -- Venturi's tube: `venturi_flow`, `continuity_velocity` [`measurement/plot_02_venturi_meter.py`]
5. 1839-1840 -- Hagen and Poiseuille: laminar pipe flow: `poiseuille_flow`, `laminar_friction_factor` [`pipe_flow/plot_01_hagen_poiseuille.py`]
6. 1845-1857 -- Weisbach and Darcy: head loss in pipes: `darcy_head_loss`, `minor_head_loss` [`pipe_flow/plot_02_darcy_weisbach.py`]
7. 1868 -- Froude's ship models: `froude_number` [`similarity/plot_01_froude_ship_models.py`]
8. 1883 -- Reynolds's dye experiment: `reynolds_number` [`pipe_flow/plot_03_reynolds_transition.py`]
9. 1889 -- Manning's open-channel formula: **new** `manning_velocity` [`open_channel/plot_01_manning_formula.py`]
10. 1898 -- Joukowsky's water hammer: `water_hammer_pressure` [`transients/plot_01_joukowsky_water_hammer.py`]
11. 1914 -- Buckingham's Pi theorem: **new** `buckingham_pi` (null space of the dimensional matrix), `euler_number`, `weber_number` [`similarity/plot_02_buckingham_pi.py`]
12. 1924 -- Thoma's cavitation parameter: `cavitation_number` [`pumps/plot_01_thoma_cavitation.py`]
13. 1933 -- Nikuradse's sand-roughened pipes: **new** `nikuradse_friction_factor` [`pipe_flow/plot_04_nikuradse_roughness.py`]
14. 1936 -- Hardy Cross's pipe-network method: **new** `hardy_cross` [`networks/plot_01_hardy_cross.py`]
15. 1939 -- Colebrook's turbulent friction factor: `colebrook_friction_factor` [`pipe_flow/plot_05_colebrook_equation.py`]
16. 1944 -- Moody's diagram: `colebrook_friction_factor`, `laminar_friction_factor` [`pipe_flow/plot_06_moody_diagram.py`]
17. 1976 -- Swamee and Jain's explicit friction factor: `swamee_jain_friction_factor` [`pipe_flow/plot_07_swamee_jain.py`]

Pump affinity laws and specific speed (`pump_specific_speed`, `pump_hydraulic_power`) can join item 12 or replace a weak entry.

### heat_transfer

1. 1701 -- Newton's law of cooling: `convection_heat_rate`, `lumped_temperature` [`convection/plot_01_newton_cooling.py`]
2. 1804 -- Biot's heated bar and the Biot number: `biot_number` [`transient/plot_01_biot_number.py`]
3. 1822 -- Fourier's *Théorie analytique de la chaleur*: `fourier_heat_rate`, `plane_wall_resistance`, `composite_wall` [`conduction/plot_01_fourier_law.py`]
4. 1859 -- Kirchhoff's law of thermal radiation: `radiation_heat_rate`, `radiation_coefficient` [`radiation/plot_01_kirchhoff_emissivity.py`]
5. 1862 -- Kelvin's age of the Earth: `semi_infinite_temperature` [`transient/plot_02_kelvin_age_of_earth.py`]
6. 1879-1884 -- Stefan and Boltzmann: `stefan_boltzmann_emissive_power` [`radiation/plot_02_stefan_boltzmann.py`]
7. 1891 -- Stefan's problem: a moving melt front: **new** `neumann_stefan_front` [`phase_change/plot_01_stefan_problem.py`]
8. 1893 -- Wien's displacement law: `blackbody_peak_wavelength` [`radiation/plot_03_wien_displacement.py`]
9. 1916 -- Nusselt's film condensation: **new** `nusselt_film_condensation` [`phase_change/plot_02_nusselt_condensation.py`]
10. 1922 -- Harper and Brown's fin efficiency: `fin_parameter`, `fin_efficiency_adiabatic`, `fin_heat_rate_adiabatic` [`conduction/plot_02_fin_efficiency.py`]
11. 1928 -- Nusselt's unit-sphere view factor: **new** `view_factor_parallel_plates`, `view_factor_coaxial_disks` [`radiation/plot_04_nusselt_view_factor.py`]
12. 1930 -- Dittus and Boelter: turbulent pipe convection: **new** `dittus_boelter_nusselt` [`convection/plot_02_dittus_boelter.py`]
13. 1934 -- Nukiyama's boiling curve and Zuber's critical heat flux (1958): **new** `zuber_critical_heat_flux` [`phase_change/plot_03_nukiyama_boiling_curve.py`]
14. 1947 -- Heisler's transient conduction charts: **new** `plane_wall_transient` (one-term and series), `fourier_number` [`transient/plot_03_heisler_charts.py`]
15. 1947 -- Crank and Nicolson: numerical conduction: **new** `conduction_1d` (explicit and Crank-Nicolson) [`conduction/plot_03_crank_nicolson.py`]
16. 1955 -- Kays and London: effectiveness-NTU: `exchanger_effectiveness_parallel`, `exchanger_effectiveness_counterflow`, `log_mean_temperature_difference` [`exchangers/plot_01_kays_london_ntu.py`]
17. 1975 -- Churchill and Chu: natural convection over all Rayleigh numbers: **new** `churchill_chu_nusselt` [`convection/plot_03_churchill_chu.py`]

The critical insulation radius (`critical_insulation_radius_*`) and `overall_u` are used in examples 3 and 16 rather than given an entry.

### materials

1. 1807 -- Young's modulus: `normal_stress`, `normal_strain`, `hooke_stress`, `StressStrainCurve` [`elasticity/plot_01_young_modulus.py`]
2. 1838 -- Duhamel's thermal stress: `thermal_strain`, `constrained_thermal_stress` [`elasticity/plot_02_duhamel_thermal_stress.py`]
3. 1852 -- Lamé's elastic constants: `shear_modulus`, `bulk_modulus`, `lame_parameters` [`elasticity/plot_03_lame_constants.py`]
4. 1864 -- Tresca's maximum-shear yield criterion: `tresca_equivalent`, `factor_of_safety` [`yield/plot_01_tresca_criterion.py`]
5. 1889-1929 -- Voigt and Reuss: composite bounds: `rule_of_mixtures_E_longitudinal`, `rule_of_mixtures_E_transverse`, `density_rule_of_mixtures` [`composites/plot_01_voigt_reuss_bounds.py`]
6. 1899 -- Goodman's mean-stress diagram: `goodman_safety_factor` [`fatigue/plot_01_goodman_diagram.py`]
7. 1910 -- Basquin's S-N law: `basquin_cycles` [`fatigue/plot_02_basquin_sn_curve.py`]
8. 1913 -- Von Mises's yield criterion: `von_mises_plane_stress` [`yield/plot_02_von_mises_criterion.py`]
9. 1920 -- Griffith's theory of brittle fracture: `griffith_fracture_stress` [`fracture/plot_01_griffith_fracture.py`]
10. 1929 -- Norton's power-law creep: `arrhenius_creep_rate` [`creep/plot_01_norton_creep.py`]
11. 1934 -- Taylor, Orowan and Polanyi: dislocations and the shear strength gap: **new** `ideal_shear_strength`, `peierls_stress` [`plasticity/plot_01_dislocations.py`]
12. 1943 -- Ramberg and Osgood's stress-strain curve: **new** `ramberg_osgood_strain` [`plasticity/plot_02_ramberg_osgood.py`]
13. 1945 -- Palmgren and Miner's cumulative damage: `miner_damage` [`fatigue/plot_03_palmgren_miner.py`]
14. 1951-1953 -- Hall and Petch: grain-size strengthening: `hall_petch` [`plasticity/plot_03_hall_petch.py`]
15. 1952 -- The Larson-Miller parameter: **new** `larson_miller_parameter` [`creep/plot_02_larson_miller.py`]
16. 1954 -- Coffin and Manson: low-cycle fatigue: **new** `coffin_manson_cycles` [`fatigue/plot_04_coffin_manson.py`]
17. 1957 -- Irwin's stress intensity factor: `stress_intensity`, `critical_crack_length` [`fracture/plot_02_irwin_stress_intensity.py`]
18. 1963 -- Paris's crack-growth law: `paris_crack_growth` [`fracture/plot_03_paris_law.py`]

`resilience_modulus` and `toughness` are used in example 1.

### mechanical_engineering

1. 1673 -- Huygens and Steiner: moments of inertia and the parallel-axis theorem: `disk_inertia`, `solid_sphere_inertia`, `parallel_axis` [`dynamics/plot_01_huygens_moment_of_inertia.py`]
2. 1765 -- Euler's involute gear tooth: `gear_ratio`, `gear_output_speed`, `gear_output_torque` [`transmission/plot_01_euler_involute_gear.py`]
3. 1769 -- Euler and Eytelwein: belt friction: **new** `belt_tension_ratio`, `belt_speed` [`transmission/plot_02_capstan_equation.py`]
4. 1785 -- Coulomb's laws of friction and the power screw: `screw_torque_raise` [`machine_elements/plot_01_coulomb_power_screw.py`]
5. 1829 -- Coriolis defines work and kinetic energy: `work_constant_force`, `kinetic_energy`, `rotational_energy`, `power_torque` [`dynamics/plot_02_coriolis_work_energy.py`]
6. 1875 -- Reuleaux's kinematics of machinery: the slider-crank: `slider_crank_displacement` [`mechanisms/plot_01_reuleaux_slider_crank.py`]
7. 1877 -- Rayleigh's *Theory of Sound*: free vibration: `natural_frequency`, `damping_ratio`, `damped_frequency`, `simulate_sdof` [`vibration/plot_01_rayleigh_free_vibration.py`]
8. 1883 -- Grashof's condition for four-bar linkages: `fourbar_grashof` [`mechanisms/plot_02_grashof_condition.py`]
9. 1883 -- Petroff's journal-bearing friction: **new** `petroff_friction_torque` [`machine_elements/plot_02_petroff_bearing.py`]
10. 1909 -- Frahm's vibration absorber and Den Hartog's tuning (1928): **new** `tuned_mass_damper` [`vibration/plot_02_frahm_absorber.py`]
11. 1919 -- Jeffcott's rotor and critical speed: `rotating_unbalance_force`, **new** `jeffcott_response` [`rotordynamics/plot_01_jeffcott_rotor.py`]
12. 1921 -- Holzer's method for torsional vibration: **new** `holzer_frequencies` [`vibration/plot_03_holzer_torsional.py`]
13. 1930s -- Base excitation and vibration isolation: `transmissibility` [`vibration/plot_04_vibration_isolation.py`]
14. 1947 -- Lundberg and Palmgren's rolling-bearing life: `bearing_l10_life` [`machine_elements/plot_03_lundberg_palmgren_l10.py`]
15. 1954 -- Freudenstein's equation for four-bar synthesis: **new** `freudenstein_coefficients` [`mechanisms/plot_03_freudenstein_equation.py`]
16. 1955 -- Denavit and Hartenberg's link parameters: **new** `dh_transform`, `forward_kinematics` [`mechanisms/plot_04_denavit_hartenberg.py`]

Item 13 has no single originator; date it from a primary source (Soderberg's 1924 isolator work is a candidate) or drop it. `flywheel_energy`, `angular_momentum`, `centripetal_force` and `natural_frequency_hz` are used inside other examples.

### reliability

1. 1693 -- Halley's life table and the hazard rate: **new** `life_table_hazard` [`lifetimes/plot_01_halley_life_table.py`]
2. 1837 -- Poisson's distribution of rare failures: `poisson_failures_probability`, `expected_failures` [`lifetimes/plot_02_poisson_failures.py`]
3. 1947 -- Freudenthal's stress-strength interference: **new** `stress_strength_reliability` [`structural/plot_01_freudenthal_interference.py`]
4. 1951 -- Weibull's distribution: `weibull_reliability`, `weibull_cdf`, `weibull_hazard`, `weibull_mttf` [`lifetimes/plot_03_weibull_distribution.py`]
5. 1952 -- Epstein and Sobel: exponential life testing: `exponential_reliability`, `exponential_mttf`, `failure_rate_from_mtbf` [`lifetimes/plot_04_epstein_sobel_life_testing.py`]
6. 1956 -- Von Neumann's redundancy: `parallel_reliability` [`systems/plot_01_von_neumann_redundancy.py`]
7. 1957 -- The AGREE report: series systems and MTBF: `series_reliability`, `exponential_hazard` [`systems/plot_02_agree_series_systems.py`]
8. 1958 -- Kaplan and Meier's estimator: **new** `kaplan_meier` [`lifetimes/plot_05_kaplan_meier.py`]
9. 1961 -- Birnbaum, Esary and Saunders: k-out-of-n systems: `k_of_n_reliability` [`systems/plot_03_k_out_of_n.py`]
10. 1962 -- Watson's fault tree (Minuteman): **new** `fault_tree_probability` [`risk/plot_01_watson_fault_tree.py`]
11. 1964 -- Duane's reliability growth: `reliability_growth_duane` [`growth/plot_01_duane_growth.py`]
12. 1965 -- Barlow and Proschan: age replacement: `renewal_interval_cost_rate` [`maintenance/plot_01_barlow_proschan_replacement.py`]
13. 1965 -- Markov models of repairable systems: `availability`, `unavailability`, `repair_rate`, `maintainability`, **new** `markov_availability` [`maintenance/plot_02_markov_availability.py`]
14. 1969 -- Birnbaum's component importance: `birnbaum_importance_series` [`systems/plot_04_birnbaum_importance.py`]
15. 1974 -- Hasofer and Lind's reliability index: **new** `hasofer_lind_index` [`structural/plot_02_hasofer_lind.py`]
16. 1975 -- WASH-1400: event trees: **new** `event_tree_probability` [`risk/plot_02_wash1400_event_tree.py`]
17. 1998 -- IEC 61508 and safety integrity levels: `safety_integrity_pfd` [`risk/plot_03_iec61508_sil.py`]

`bathtub_hazard` illustrates entry 1, with Halley's infant and old-age mortality.

### signal_processing

1. 1898 -- Schuster's periodogram: `dft`, `idft`, `power_spectral_density` [`spectral/plot_01_schuster_periodogram.py`]
2. 1899 -- Gibbs's phenomenon: **new** `fourier_partial_sum` [`fourier/plot_01_gibbs_phenomenon.py`]
3. 1910 -- Haar's wavelet: **new** `haar_transform` [`wavelets/plot_01_haar_wavelet.py`]
4. 1930 -- Butterworth's maximally flat filter: `butterworth_lowpass`, `butterworth_highpass`, `frequency_response` [`filters/plot_01_butterworth.py`]
5. 1934 -- The Wiener-Khinchin theorem: `correlation`, `power_spectral_density` [`spectral/plot_02_wiener_khinchin.py`]
6. 1946 -- Gabor's time-frequency atoms: **new** `stft` [`spectral/plot_03_gabor_stft.py`]
7. 1947 -- Hurewicz's z-transform: `frequency_response`, `filter_signal` [`filters/plot_02_hurewicz_z_transform.py`]
8. 1948 -- Bennett's quantization noise: `quantize`, `quantization_step`, `snr_db` [`sampling/plot_01_bennett_quantization.py`]
9. 1949 -- Shannon's sampling theorem and aliasing: `sinc_reconstruct`, `alias_frequency`, `nyquist_frequency` [`sampling/plot_02_shannon_sampling.py`]
10. 1949 -- Wiener's optimal filter: **new** `wiener_filter` [`filters/plot_03_wiener_filter.py`]
11. 1958 -- Blackman and Tukey's spectral windows: `windowed_spectrum` [`spectral/plot_04_blackman_tukey_windows.py`]
12. 1958 -- Goertzel's algorithm: **new** `goertzel` [`spectral/plot_05_goertzel.py`]
13. 1960 -- Widrow and Hoff's LMS adaptive filter: **new** `lms_filter` [`filters/plot_04_widrow_hoff_lms.py`]
14. 1965 -- Cooley and Tukey's FFT: `fft_frequencies`, `dft` (count operations against a hand-written radix-2) [`fourier/plot_02_cooley_tukey_fft.py`]
15. 1966 -- Kaiser's window method for FIR design: `fir_lowpass` [`filters/plot_05_kaiser_fir.py`]
16. 1967 -- Welch's averaged periodogram: **new** `welch_psd` [`spectral/plot_06_welch_psd.py`]
17. 1970 -- Jackson: cascaded second-order sections: `butterworth_lowpass(..., output="sos")` (item 3) [`filters/plot_06_second_order_sections.py`]

`convolution`, `moving_average`, `rms`, `energy` and `power` are used in examples 7 and 14.

### structural_engineering

1. 1638 -- Galileo's cantilever: `rectangular_section_modulus`, `rectangular_second_moment` [`beams/plot_01_galileo_cantilever.py`]
2. 1678 -- Hooke's *Ut tensio, sic vis*: axial members: `axial_deformation`, `strain_energy_axial` [`members/plot_01_hooke_axial.py`]
3. 1744 -- Euler's column buckling: `euler_buckling_load`, `radius_of_gyration`, `slenderness_ratio` [`stability/plot_01_euler_buckling.py`]
4. 1784 -- Coulomb's torsion: `torsion_shear_stress`, `torsion_angle` [`members/plot_02_coulomb_torsion.py`]
5. 1826 -- Navier's beam theory: `bending_stress`, `SimplySupportedBeam` [`beams/plot_02_navier_beam.py`]
6. 1829 -- Poisson's ratio: `poisson_lateral_strain` [`members/plot_03_poisson_ratio.py`]
7. 1847 -- Whipple's method of joints: **new** `method_of_joints` [`trusses/plot_01_whipple_method_of_joints.py`]
8. 1857 -- Clapeyron's three-moment equation: **new** `three_moment_equation` [`beams/plot_03_clapeyron_three_moment.py`]
9. 1864 -- Maxwell's rigidity count and the reciprocal theorem: `maxwell_betti_reciprocity`, **new** `truss_determinacy` [`trusses/plot_02_maxwell_rigidity.py`]
10. 1867 -- Winkler's influence lines: **new** `influence_line` [`beams/plot_04_winkler_influence_lines.py`]
11. 1879 -- Castigliano's theorem: `castigliano_deflection` (integrating itself, item 3), `strain_energy_bending_constant_moment` [`energy/plot_01_castigliano.py`]
12. 1882 -- Mohr's circle: `mohr_circle` [`stress/plot_01_mohr_circle.py`]
13. 1889 -- Engesser's tangent-modulus load: **new** `tangent_modulus_load` [`stability/plot_02_engesser_tangent_modulus.py`]
14. 1914 -- Kazinczy's plastic hinge: **new** `plastic_section_modulus`, `plastic_collapse_load` [`plasticity/plot_01_kazinczy_plastic_hinge.py`]
15. 1930 -- Hardy Cross's moment distribution: **new** `moment_distribution` [`frames/plot_01_hardy_cross_moment_distribution.py`]
16. 1956 -- Turner, Clough, Martin and Topp: the direct stiffness method: `PlanarTruss` (with conditioning diagnostics, item 3) [`trusses/plot_03_direct_stiffness.py`]
17. 1956 -- Plane frames by direct stiffness: **new** `PlanarFrame` [`frames/plot_02_planar_frame.py`]

Entry 17 is a candidate for dropping if 16 carries the stiffness story on its own. `hooke_stress`, `thermal_strain` and `parallel_axis` duplicate `materials` and `mechanical_engineering` and are not given entries here.

### thermodynamics

1. 1824 -- Carnot's *Réflexions*: `carnot_efficiency`, `heat_engine_work` [`cycles/plot_01_carnot_efficiency.py`]
2. 1834 -- Clapeyron's ideal gas law: `ideal_gas_pressure`, `ideal_gas_temperature` [`properties/plot_01_clapeyron_ideal_gas.py`]
3. 1834-1850 -- The Clausius-Clapeyron relation: **new** `clausius_clapeyron_pressure` [`properties/plot_02_clausius_clapeyron.py`]
4. 1834 -- Perkins's vapor-compression refrigerator: `refrigeration_cop`, `carnot_cop_refrigerator` [`refrigeration/plot_01_perkins_vapor_compression.py`]
5. 1843-1850 -- Joule and the first law: `first_law_heat`, `specific_internal_energy_ideal_gas`, `specific_enthalpy_ideal_gas` [`laws/plot_01_joule_first_law.py`]
6. 1852 -- Kelvin's heat pump: `heat_pump_cop`, `carnot_cop_heat_pump` [`refrigeration/plot_02_kelvin_heat_pump.py`]
7. 1852 -- The Joule-Thomson effect: **new** `joule_thomson_coefficient_vdw` [`properties/plot_03_joule_thomson.py`]
8. 1854-1865 -- Clausius's inequality and entropy: `clausius_inequality`, `entropy_change_isothermal_ideal_gas`, `reversible_isothermal_work_ideal_gas`, `isentropic_pressure_volume`, `isentropic_temperature_ratio` [`laws/plot_02_clausius_entropy.py`]
9. 1859 -- Rankine's steam cycle: `rankine_ideal_efficiency` (with IF97 properties, item 3) [`cycles/plot_02_rankine_cycle.py`]
10. 1872 -- Brayton's constant-pressure cycle: `brayton_efficiency` [`cycles/plot_03_brayton_cycle.py`]
11. 1873 -- Van der Waals's equation of state: **new** `van_der_waals_pressure` [`properties/plot_04_van_der_waals.py`]
12. 1876 -- Otto's four-stroke cycle: `OttoCycle`, `otto_efficiency` [`cycles/plot_04_otto_cycle.py`]
13. 1888 -- Antoine's vapor-pressure equation: **new** `antoine_vapor_pressure` [`properties/plot_05_antoine_equation.py`]
14. 1889-1905 -- Gouy and Stodola: lost work: `exergy_destruction` [`exergy/plot_01_gouy_stodola.py`]
15. 1893 -- Diesel's compression-ignition cycle: `diesel_efficiency` [`cycles/plot_05_diesel_cycle.py`]
16. 1911 -- Carrier's psychrometric chart: **new** `humidity_ratio`, `relative_humidity`, `dew_point` [`psychrometrics/plot_01_carrier_psychrometrics.py`]
17. 1956 -- Rant names exergy: `exergy_heat` [`exergy/plot_02_rant_exergy.py`]
18. 1975 -- Curzon and Ahlborn: efficiency at maximum power: **new** `curzon_ahlborn_efficiency` [`cycles/plot_06_curzon_ahlborn.py`]
19. 1997 -- IAPWS-IF97 steam properties (regions 1, 2 and 4): **new** `steam_properties` [`properties/plot_06_iapws_if97.py`]

`thermal_efficiency` is used in examples 1 and 12-15.

## Done

Finished items move here with the release they shipped in, as in
CHANGELOG.md.

### Unreleased

**`aerodynamics` rebuilt.** The domain moved to the `core/systems/utils/
visualizers` layout. Its history page has 21 dated breakthroughs, from
Newton and Laplace's speed of sound to Oswald's drag polar, each with its
own gallery example in eight topic sections.

**`communications` rebuilt.** Same layout. 17 breakthroughs, from Hertz's
radio waves to OFDM, each with its own example in eight topic sections.

**`control_systems` partly rebuilt.** It has the new layout (`systems/`
split into transfer, frequency, stability, root locus, PID, state space,
optimal and sampled modules) and 13 examples. The history page is still
to be written; see the phase checklist.

### 0.1.0

Physical input contracts and nonfinite checks for the scalar models,
independent numerical benchmarks in all twelve domains, input/output units
in the docstrings and twelve assumptions/validity guides. 100% production
statement coverage is enforced in CI. Statement coverage does not show
that the numerics are correct: branch coverage, conditioning checks and
reference problems (items 2 and 3) are the next step.

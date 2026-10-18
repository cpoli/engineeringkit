# Changelog

The release workflow takes a release's notes from its `## [X.Y.Z] - date`
section, so keep that heading format.

## [0.1.0] - 2026-10-18

- Add release infrastructure: `CITATION.cff`, a tag-driven release workflow
  (PyPI trusted publishing, GitHub Release, docs rebuild), a docs workflow
  that deploys `gh-pages` on every push to `main`, and a CI matrix over
  Python 3.10-3.14 on Linux and macOS.
- Reset the public release version to 0.1.0 at the maintainer's request.
- Reach 100% production statement coverage: 725 tests, zero missed statements.
- Enforce 100% statement coverage in the configured coverage report and CI.
- Add analytical/reference tests for all previously uncovered numerical helpers,
  result summaries, plots, animations, and physical-input error paths.
- Fix the SDOF animation interval argument and integer-input digital filtering.
- Remove unreachable duplicate validation guards without changing input contracts.
- Retain all twelve domains, validated models and detailed documentation.
- `isentropic_pressure_volume`, `isentropic_temperature_ratio`,
  `otto_efficiency`, `diesel_efficiency` and `brayton_efficiency` accept
  arrays for `gamma` (and the Diesel and Brayton ratios), as the
  aerodynamics relations already did. Scalar calls are unchanged.
- `steam_properties` covers IAPWS-IF97 region 3 (near the critical point):
  states given by `T`, `s` or `h` above 350 °C and the B23 pressure, and
  wet states up to the critical pressure, instead of raising `ValueError`.
  `plot_saturation_dome` draws the dome up to the critical point by default.
- Complete the NumPy-style docstrings of the result-class methods, the
  shared `integrators`, `visualizers` and `units` packages, and the
  IAPWS-IF97 helpers.
- Rebuild `structural_engineering`: a dated history page of 17
  breakthroughs, each with its own gallery example; add
  `three_moment_equation`, `influence_line`, `truss_determinacy`,
  `method_of_joints`, `plastic_section_modulus`, `plastic_collapse_load`,
  `tangent_modulus_load`, `stress_transformation`, `moment_distribution`
  and `PlanarFrame`. `castigliano_deflection` also integrates
  `M dM/dP / EI` along a beam, and `PlanarTruss.solve` checks Maxwell's
  count, returns `cond(K_ff)` and warns near a mechanism. Existing calls
  are unchanged.
- Rebuild `reliability`: a dated history page of 17 breakthroughs, each
  with its own gallery example; add `life_table_hazard`,
  `mtbf_confidence_bounds`, `kaplan_meier`, `stress_strength_reliability`,
  `hasofer_lind_index`, `failure_probability_from_index`,
  `reliability_index_from_probability`, `fault_tree_probability`,
  `event_tree_probability`, `safety_integrity_level`,
  `markov_availability`, `markov_state_probabilities` and
  `age_replacement_cost_rate`. Existing functions keep their signatures;
  `safety_integrity_pfd` now accepts arrays.
- Rebuild `signal_processing`: a dated history page of 17 breakthroughs,
  each with its own gallery example; add `radix2_fft`,
  `fourier_partial_sum`, `haar_transform`, `inverse_haar_transform`,
  `welch_psd`, `stft`, `goertzel`, `kaiser_order`, `wiener_filter` and
  `lms_filter`. `butterworth_lowpass` and `butterworth_highpass` take
  `output="sos"`, and `filter_signal` and `frequency_response` accept
  second-order sections with `a=None`; existing calls are unchanged.
- Rebuild `control_systems`: a dated history page of 21 breakthroughs, each
  with its own gallery example; add `pade_delay`, `smith_predictor` and
  `jury_stable`. `stability_margins` now reads every crossover and returns
  the worst margins (same scalar keys) plus arrays of all of them, so
  conditionally stable and open-loop unstable loops come out right.
- Rebuild `mechanical_engineering`: a dated history page of 16
  breakthroughs, each with its own gallery example; add
  `belt_tension_ratio`, `petroff_friction_torque`, `sommerfeld_number`,
  `long_bearing_eccentricity`, `fourbar_output_angle`,
  `freudenstein_coefficients`, `freudenstein_link_lengths`, `dh_transform`,
  `forward_kinematics`, `tuned_mass_damper`, `holzer_residual`,
  `holzer_frequencies` and `jeffcott_response`. Existing functions keep
  their signatures; `damped_frequency`, `slider_crank_displacement` and
  `screw_torque_raise` now accept arrays.
- Rebuild `heat_transfer`: a dated history page of 17 breakthroughs, each
  with its own gallery example; add `conduction_1d`, `plane_wall_transient`,
  `dittus_boelter_nusselt`, `rayleigh_number`, `churchill_chu_nusselt`,
  `planck_spectral_emissive_power`, `view_factor_parallel_plates`,
  `view_factor_coaxial_disks`, `stefan_number`, `neumann_stefan_front`,
  `nusselt_film_condensation` and `zuber_critical_heat_flux`. Existing
  functions keep their signatures; the radial resistances, fin efficiency
  and radiation helpers now accept arrays.
- Rebuild `fluid_engineering`: a dated history page of 17 breakthroughs,
  each with its own gallery example; add `hydrostatic_force_plane`,
  `center_of_pressure_depth`, `poiseuille_velocity`, `pipe_resistance`,
  `smooth_pipe_friction_factor`, `nikuradse_friction_factor`,
  `buckingham_pi`, `hydraulic_radius`, `manning_velocity`, `manning_flow`,
  `pipe_wave_speed`, `critical_closure_time`, `npsh_available`,
  `thoma_sigma` and `hardy_cross`, the `PiGroups` and `NetworkResult`
  result classes, and the `plot_moody_diagram` and `plot_grade_lines`
  visualizers. Existing functions keep their signatures; `venturi_flow`
  gains an optional discharge coefficient, and the friction factors,
  `euler_number` and `cavitation_number` now accept arrays. Colebrook's
  equation is solved to machine precision.
- Rebuild `materials`: a dated history page of 18 breakthroughs, each
  with its own gallery example; add `coffin_manson_cycles`,
  `strain_life_amplitude`, `crack_growth_life`, `larson_miller_parameter`,
  `larson_miller_rupture_time`, `ideal_shear_strength`, `peierls_stress`,
  `ramberg_osgood_strain`, `StressStrainCurve.offset_yield_strength` and
  `.summary()`, the `CrackGrowthResult` result class, `utils.power_law_fit`
  and the `plot_stress_strain`, `plot_yield_surfaces`, `plot_sn_curve` and
  `plot_crack_growth` visualizers. Existing functions keep their
  signatures; all now accept arrays, `tresca_equivalent` reduces along the
  last axis, and `thermal_strain` and `von_mises_plane_stress` validate
  their inputs. `StressStrainCurve.plot` labels stress in MPa.
- Rebuild `thermodynamics`: a dated history page of 19 breakthroughs, each
  with its own gallery example; add `van_der_waals_pressure`,
  `van_der_waals_critical_point`, `joule_thomson_coefficient_vdw`,
  `clausius_clapeyron_pressure`, `antoine_vapor_pressure`,
  `saturation_pressure`, `saturation_temperature`, `steam_properties`
  (IAPWS-IF97 regions 1, 2 and 4) with the `SteamState` result class,
  `curzon_ahlborn_efficiency`, `humidity_ratio`, `relative_humidity`,
  `dew_point` and `MOLAR_MASS_RATIO`, and the `plot_pv_diagram`,
  `plot_ts_diagram` and `plot_saturation_dome` visualizers. Existing
  functions keep their signatures.
- Remove the superseded per-domain `plot_<domain>_models.py` and
  `plot_<domain>_benchmark.py` examples, `examples/all_domains.py`,
  `examples/phase1_gallery.py`, `gallery/integrators.png` and the
  uncollected top-level `tests/test_phase1.py`.
- Add course material: four cross-domain tutorials, an exercise set with
  worked solutions for each of the twelve domains, and a quickstart
  notebook (Colab badge, served in JupyterLite). The docs build runs every
  tutorial and solution, and the test suite checks the numbers they print.

## Unpublished development milestones

Earlier 0.2–0.5 labels were development iterations, not external publications.
The initial public release is now 0.1.0. Full historical development notes are
available in earlier source versions; all implementation improvements are retained.

### Preserved development notes

### 0.5.0

- Add physical-domain/nonfinite input checks across the twelve engineering domains.
- Correct SeriesRLC charge/current mismatch and FM integration phase offset.
- Fix unsigned BPSK/QPSK mapping and validate binary input and matching shapes.
- Preserve tiny failure probabilities, large Poisson counts, and near-singular cycle/exchanger limits.
- Add Riccati matrix/weight, truss DOF/connectivity and stress-strain data validation.
- Add analytical, conservation, quadrature, distribution and invalid-domain benchmarks.
- Document scalar public API input/output units and twelve detailed model-validity guides.
- Add warnings-as-errors test policy and production-only coverage configuration.

### 0.4.0

- Implement shared Euler/RK4/adaptive integration with inspectable histories.
- Implement reusable state plots and displacement animations.
- Validate passive oscillator inputs and solver success; expose history/summary.
- Fix PID zero-integral cancellation, Routh zero rows and marginal stability.
- Validate root-locus polynomial degrees; fix complex SNR and mid-rise ADC levels.
- Handle zero m*L fin efficiency through its analytic limit.
- Format source, repair lint checks, include py.typed and expand Sphinx navigation/API.
- Add analytical, conservation, convergence, invalid-input and rendered-media tests.



### 0.3.0
- Completed all Phase 2 and Phase 3 engineering domains.
- Added heat transfer, materials, signal processing, mechanical engineering, aerodynamics, communications, and reliability implementations.
- Added domain-local numerical tests, documentation progression pages, and an all-domain executable gallery.


### 0.2.0
- Expanded all five Phase-1 domains with additional educational methods.
- Added direct-stiffness planar trusses, pole placement and Kalman estimation.
- Added RC/RL network helpers, Rankine/exergy relations, fluid similarity numbers and water hammer.
- Moved domain tests alongside subpackages and added history pages plus a runnable gallery.

### 0.1.0
- Initial educational foundation and five Phase-1 domains.

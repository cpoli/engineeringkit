# Changelog

The release workflow takes a release's notes from its `## [X.Y.Z] - date`
section, so keep that heading format.

## [Unreleased]

## [0.2.0] - 2026-10-30

- A weekly `Linkcheck` workflow checks every external link in the docs,
  and the history pages' DOIs are links to doi.org (ROADMAP B3).
- The CI `typecheck` job runs `python -m mypy` on the whole package with
  the `pyproject.toml` settings, instead of the integrators only
  (ROADMAP A3).
- `thermodynamics` covers the rest of a first engineering
  thermodynamics course (ROADMAP D12): steady-flow devices
  (`steady_flow_work`, `nozzle_exit_velocity`, `mixing_chamber`,
  `throttle_steam`, the turbine, compressor and nozzle isentropic
  efficiencies and `expansion_outlet_enthalpy` and
  `compression_outlet_enthalpy`); `rankine_cycle` (a `RankineResult`) with
  superheat, reheat, an open feedwater heater and component efficiencies,
  and `brayton_cycle` (a `BraytonResult`) with regeneration, intercooling
  and reheat; `stirling_efficiency` and `dual_efficiency`; the
  Peng-Robinson equation (`peng_robinson_pressure`,
  `peng_robinson_compressibility`, `peng_robinson_saturation_pressure`),
  `CubicFluid`, `R134A`, `refrigerant_properties` (a `RefrigerantState`)
  and `vapor_compression_cycle` (a `VaporCompressionResult`); Dalton's
  mixtures (`mole_fractions`, `mass_fractions`, `mixture_molar_mass`,
  `partial_pressures`, `entropy_of_mixing`); `compressibility_factor`,
  `pitzer_acentric_factor` and `pitzer_virial_compressibility`; Maxwell's
  relations as `cp_minus_cv` and `joule_thomson_coefficient`;
  psychrometric processes (`moist_air_enthalpy`,
  `adiabatic_saturation_temperature`, `adiabatic_mixing`, `cooling_coil`,
  `cooling_tower`); and combustion (`fuel_molar_mass`,
  `stoichiometric_air_fuel_ratio`, `combustion_products`, `heating_value`,
  `adiabatic_flame_temperature`, with NIST Shomate enthalpies in
  `utils.shomate`). `steam_properties` adds IAPWS-IF97 region 5 (1073.15 K
  to 2273.15 K up to 50 MPa) and accepts arrays, returning a `SteamState`
  of arrays; states it used to reject with `ValueError` above 1073.15 K,
  or because an input was an array, now return values. Scalar calls are
  otherwise unchanged. The history page gains six breakthroughs (Dalton
  1801, Stirling 1816, Hess 1840, Maxwell 1871, Pitzer 1955, Peng and
  Robinson 1976), each with its own gallery example, and seven existing
  entries a second example each.
- `aerodynamics` covers the rest of a first aerodynamics and flight
  course (ROADMAP D1): oblique shocks (`oblique_shock_deflection`,
  `oblique_shock_max_deflection`, `oblique_shock_angle` and
  `oblique_shock`, returning an `ObliqueShockResult`); the Prandtl-Meyer
  expansion (`prandtl_meyer_angle`, `inverse_prandtl_meyer`,
  `prandtl_meyer_expansion`); Fanno and Rayleigh duct flow (`fanno_flow`,
  `rayleigh_flow`, both returning a `DuctFlowState`, with
  `fanno_choking_length` and `rayleigh_choking_heat`); a general
  lifting-line solver by Glauert's Fourier series (`lifting_line`, a
  `LiftingLineResult`); the Hess-Smith panel method (`panel_method`, a
  `PanelResult`); `naca4_moment_coefficient` and `aerodynamic_center`;
  static longitudinal stability (`tail_volume_ratio`,
  `downwash_gradient`, `neutral_point`, `pitching_moment_slope`,
  `static_margin`); and climb, turn and take-off performance
  (`rate_of_climb`, `climb_angle`, `load_factor`, `turn_radius`,
  `turn_rate`, `takeoff_ground_roll`) with `breguet_range_propeller`.
  New plots: `plot_theta_beta_mach`, `plot_duct_flow`,
  `plot_spanwise_loading` and `plot_pressure_coefficient`. The history
  page gains six breakthroughs (Fanno 1904, Prandtl-Meyer 1907-1908,
  Meyer's oblique shock 1908, Rayleigh 1910, Glauert 1926, Hess and Smith
  1967), each with its own gallery example, and the thin-airfoil and
  Breguet entries a second example each (pitching moment, propeller
  range). Existing calls are unchanged.
- `communications` covers the rest of a first communication systems
  course (ROADMAP D2): the matched filter (`matched_filter`,
  `matched_filter_snr`); Gray-labelled M-PSK and square M-QAM
  (`psk_constellation`, `qam_constellation`, `minimum_distance_detect`)
  with exact symbol-error rates (`ser_mpsk_awgn`, `ser_mqam_awgn`);
  Nyquist's ISI criterion with raised-cosine pulses
  (`raised_cosine_pulse`, `raised_cosine_spectrum`,
  `root_raised_cosine_pulse`, `raised_cosine_bandwidth`); PCM, mu-law
  companding and delta modulation (`pcm_encode`, `pcm_decode`,
  `pcm_sqnr_db`, `mu_law_compress`, `mu_law_expand`, `delta_modulate`,
  `delta_demodulate`, `delta_slope_limit`); `entropy` and `huffman_code`
  (a `HuffmanCode`); CRCs (`crc_remainder`, `crc_encode`, `crc_check`,
  `CRC16_CCITT`); convolutional codes with Viterbi decoding
  (`convolutional_encode`, `convolutional_free_distance`,
  `viterbi_decode`, a `ViterbiDecodeResult`); Rayleigh fading and
  diversity (`rayleigh_pdf`, `ber_bpsk_rayleigh`,
  `ber_bpsk_mrc_rayleigh`, `selection_outage_probability`,
  `maximal_ratio_combine`); antennas (`effective_aperture`,
  `aperture_gain`, `directivity`, `half_wave_dipole_pattern`,
  `half_wave_dipole_directivity`, `half_wave_dipole_resistance`,
  `FREE_SPACE_IMPEDANCE`); and the noise performance of AM and FM
  (`baseband_snr`, `am_output_snr`, `fm_output_snr`,
  `fm_threshold_snr`). New plots: `plot_eye_diagram` and
  `plot_antenna_pattern`. The history page gains thirteen breakthroughs
  (Rayleigh 1880, Carter 1932, Crosby 1937, Reeves 1937-1938, North 1943,
  delta modulation 1946-1952, Huffman 1952, Kahn and Brennan 1954-1959,
  Elias 1955, Cahn 1959, Peterson and Brown 1961, Campopiano and Glazer
  1962, Viterbi 1967), each with its own gallery example, and the Nyquist
  and Friis transmission entries a second example each (raised-cosine
  pulses, effective aperture). Existing calls are unchanged.
- `heat_transfer` covers the rest of a first heat and mass transfer course
  (ROADMAP D6): `flat_plate_nusselt` (laminar, turbulent and mixed),
  `churchill_bernstein_nusselt`, `zukauskas_nusselt` and
  `tube_bank_max_velocity`; `laminar_tube_nusselt`, `entry_length` and
  `hausen_nusselt`; `lewis_number` and `chilton_colburn_mass_transfer`;
  `fin_heat_rate_convective_tip`; `conduction_2d` (a
  `Conduction2DResult`) and the shape factors
  `shape_factor_buried_sphere`, `shape_factor_buried_cylinder`,
  `shape_factor_eccentric_cylinders` and `shape_factor_parallel_cylinders`;
  `cylinder_transient` and `sphere_transient` (with
  `utils.cylinder_eigenvalues` and `utils.sphere_eigenvalues`);
  `radiation_enclosure` (an `EnclosureResult`), `two_surface_radiation`,
  `radiation_shield_heat_flux` and `view_factor_reciprocity`;
  `exchanger_effectiveness_crossflow`, `exchanger_effectiveness_shell_tube`
  and `lmtd_correction_factor`; `rohsenow_boiling_heat_flux`. Eight new
  history entries and gallery examples (17 to 25 breakthroughs). Existing
  functions are unchanged.
- `mechanical_engineering` covers the rest of a first machine-design and
  vibrations course (ROADMAP D8): `normal_modes` (returning a
  `ModalResult`), `modal_response` (an `MDOFResult`) and
  `harmonic_response`; `magnification_factor` and `phase_angle`;
  `rayleigh_frequency` and `dunkerley_frequency`; `lewis_form_factor`,
  `lewis_bending_stress` and `barth_velocity_factor`;
  `shaft_von_mises_stress`, `shaft_max_shear_stress`, `shaft_diameter`
  and `shaft_diameter_goodman`; `spring_index`, `spring_rate`,
  `wahl_factor`, `bergstrasser_factor` and `spring_shear_stress`;
  `bolt_stiffness`, `member_stiffness`, `bolted_joint` (a `BoltedJoint`)
  and `preload_torque`; `unbalance` and `balance_rotating_masses` (a
  `BalanceResult`); `fourbar_kinematics` and `slider_crank_kinematics`
  (`FourBarMotion`, `SliderCrankMotion`); `cam_follower_motion` (a
  `CamMotion`); and `lagrange_accelerations` and `simulate_lagrangian`.
  The history page gains eight breakthroughs (Bernoulli 1753, Lagrange
  1788, Aronhold and Kennedy 1872-1886, Lewis 1892, Dunkerley 1894, Dalby
  1902, Rötscher 1927, Wahl 1929), each with its own gallery example, and
  Rayleigh's entry a forced-vibration example. Existing calls are
  unchanged.
- `electrical_engineering` covers the rest of a first circuits and power
  course (ROADMAP D4): `mesh_analysis` returning a `MeshResult`;
  `rl_current`, `rlc_damping` and the over-, critically and underdamped
  `rlc_step_voltage`; ideal op amps (`inverting_amplifier_gain`,
  `noninverting_amplifier_gain`, `op_amp_integrator`); the Shockley diode
  (`thermal_voltage`, `shockley_diode_current`) and its load-line
  `diode_operating_point`; DC machine (`dc_motor_torque`,
  `dc_motor_speed`) and induction motor (`induction_motor_torque`,
  `induction_motor_breakdown`) torque-speed curves; the per-unit system
  (`per_unit_base`, `change_base`); `fault_currents` for three-phase,
  line-to-ground, line-to-line and double-line-to-ground faults;
  `newton_raphson_load_flow` with PV buses; and the swing equation
  (`swing_response` returning a `SwingResult`) with the equal-area
  criterion (`equal_area_clearing_angle`, `critical_clearing_time`). The
  power-system models stay in this domain rather than seeding a
  `power_systems` domain. The history page gains eight breakthroughs
  (Helmholtz 1851, Thomson 1853, the Hopkinsons 1886, Steinmetz 1897,
  Evans and Wagner 1926, Ragazzini 1947, Shockley 1949, Tinney and Hart
  1967), each with its own gallery example. Existing calls are unchanged.
- `control_systems` covers the rest of a first course (ROADMAP D3):
  lead, lag and lead-lag compensators (`lead_compensator`,
  `lag_compensator`, `lead_lag_compensator`) and their frequency-response
  design (`design_lead`, `design_lag`, `design_lead_lag` returning a
  `CompensatorResult`); `system_type` and `error_constants`;
  `sensitivity`, `complementary_sensitivity` and `sensitivity_peak`;
  the Nichols chart (`nichols`, `nichols_closed_loop`, `plot_nichols`);
  `discretize_tustin` with prewarping; `PID.simulate_saturated` with
  clamping or back-calculation anti-windup; describing functions
  (`describing_function`, `relay_describing_function`,
  `saturation_describing_function`, `dead_zone_describing_function`),
  `limit_cycle`, and `phase_portrait` with `plot_phase_portrait`.
  `ResponseResult` gains an optional actuator signal `u`. The history page
  gains five breakthroughs (Andronov, Vitt and Khaikin 1937, Bode 1945,
  Nichols 1947, Tustin 1947, Kochenburger 1950), each with its own gallery
  example. Existing calls are unchanged.
- `fluid_engineering` covers the rest of a first course (ROADMAP D5):
  the momentum equation for jets, vanes and pipe bends (`jet_force`,
  `bend_force`, `pelton_wheel_power`); pump and system curves, the
  operating point of pumps in series and parallel (`pump_curve_head`,
  `system_curve_head`, `pump_operating_point` returning an
  `OperatingPoint`), the affinity laws and Euler's turbomachine head;
  specific energy, critical and alternate depths and the hydraulic jump
  (`specific_energy`, `critical_depth`, `alternate_depths`,
  `sequent_depth`, `hydraulic_jump_head_loss`); orifices, tank draining
  and sharp-crested weirs (`orifice_flow`, `tank_drain_time`,
  `rectangular_weir_flow`, `rehbock_discharge_coefficient`,
  `v_notch_weir_flow`); Stokes drag, the sphere drag curve and terminal
  velocity; and Couette flow with the rotating-cylinder viscometer. The
  history page gains eight breakthroughs (Torricelli 1644, Euler 1754,
  Bélanger 1828, Stokes 1851, Thomson 1858, Pelton 1880, Couette 1890,
  Bakhmeteff 1912), each with its own gallery example. Existing calls are
  unchanged.
- `materials` covers the rest of a first course (ROADMAP D7): crystal
  structures and Bragg diffraction (`cubic_lattice_parameter`,
  `atomic_packing_factor`, `theoretical_density`, `interplanar_spacing`,
  `bragg_angle`, `diffraction_peaks`); Fick's laws and carburizing
  (`diffusion_coefficient`, `fick_flux`, `semi_infinite_concentration`,
  `diffusion_time`); `lever_rule` and Avrami kinetics (`avrami_fraction`,
  `avrami_time`, `avrami_fit`); `brinell_hardness`, `vickers_hardness` and
  `brinell_tensile_strength`; `gerber_safety_factor`,
  `soderberg_safety_factor`, `endurance_limit_estimate` and the Marin
  factors (`marin_surface_factor`, `marin_size_factor`,
  `marin_load_factor`, `marin_reliability_factor`,
  `marin_endurance_limit`); notches (`elliptical_hole_concentration`,
  `peterson_notch_sensitivity`, `neuber_notch_sensitivity`,
  `fatigue_notch_factor`, `neuber_notch_stress`); three-dimensional
  `principal_stresses` and `von_mises_stress`; `halpin_tsai` and
  classical laminate theory (`lamina_stiffness`, `transformed_stiffness`,
  `classical_laminate` returning a `Laminate`); and Weibull brittle
  strength (`weibull_failure_probability`, `weibull_size_effect`,
  `weibull_mean_strength`). The history page gains eight breakthroughs
  (Fick 1855, Brinell 1900, Bragg 1913, Avrami 1939, Weibull 1939, Neuber
  1961, Reissner and Stavsky 1961, Halpin and Tsai 1969), each with its own
  gallery example. Existing calls are unchanged.
- `signal_processing` covers the rest of a first course (ROADMAP D10):
  Chebyshev type I and II and elliptic low-pass filters
  (`chebyshev1_lowpass`, `chebyshev2_lowpass`, `elliptic_lowpass`) with
  minimum-order formulas (`iir_order`); the bilinear transform written out,
  with prewarping (`bilinear_transform`, `prewarp_frequency`);
  Parks-McClellan equiripple FIR design (`equiripple_lowpass`,
  `equiripple_order`); circular convolution and overlap-add
  (`circular_convolution`, `overlap_add`); sample-rate conversion
  (`downsample`, `upsample`, `decimate`, `interpolate`); the Hilbert
  transform and analytic signal (`analytic_signal`, `hilbert_transform`,
  `envelope`, `instantaneous_frequency`, `hilbert_fir`); autoregressive
  models (`levinson_durbin`, `yule_walker` returning an `ARModel` with its
  spectrum); the orthonormal DCT (`dct`, `idct`); and the real cepstrum
  (`real_cepstrum`). The guide covers zero padding. The history page gains
  eleven breakthroughs (Chebyshev 1854, Hilbert 1905, Yule 1927, Cauer
  1931, Durbin 1960, Bogert, Healy and Tukey 1963, Golden and Kaiser 1964,
  Stockham 1966, Parks and McClellan 1972, Schafer and Rabiner 1973,
  Ahmed, Natarajan and Rao 1974), each with its own gallery example.
  Existing calls are unchanged.
- Array-input policy (ROADMAP C4), written in CONTRIBUTING.md: scalar
  functions broadcast their `float` and `complex` parameters over arrays.
  `engineeringkit/tests/test_array_inputs.py` calls every scalar function
  with 1-D arrays built from its own test inputs. Now accept arrays:
  `carnot_efficiency`, `carnot_cop_heat_pump`, `carnot_cop_refrigerator`,
  `heat_engine_work`, `thermal_efficiency`, `heat_pump_cop`,
  `rankine_ideal_efficiency`, `log_mean_temperature_difference`,
  `exchanger_effectiveness_counterflow`, `power_factor`,
  `thevenin_from_source`, `wheatstone_balanced`,
  `distortionless_condition`, `maxwell_betti_reciprocity`,
  `ziegler_nichols`, `series_reliability`, `parallel_reliability`,
  `mtbf_confidence_bounds`, `poisson_failures_probability`,
  `quantization_step` and `quantize`. Scalar calls give the same values;
  some now return NumPy scalars (`np.float64`, still a `float`) instead of
  Python floats. Scalar-only parameters are documented in
  `naca4_zero_lift_angle`, `plane_wall_transient` (`Bi`),
  `neumann_stefan_front` (`Ste`), `crack_growth_life`,
  `jeffcott_response` (all but `omega`), `tuned_mass_damper` and
  `kaiser_order`. `carnot_efficiency` no longer rejects `T_cold <= 0`
  separately; the positivity check already did.
- `structural_engineering` covers the rest of a first mechanics of
  materials and structural analysis course (ROADMAP D11): `Beam`, a beam
  on any supports solved by Macaulay's singularity functions, and
  `superpose_beam_results`; `composite_section`, `i_section`,
  `t_section`, `channel_section`, `angle_section` and `box_section`
  returning a `SectionProperties` (transformed sections by modular
  ratio, principal axes, `first_moment`, `width_at`);
  `unsymmetrical_bending_stress` and `neutral_axis_angle`;
  `transverse_shear_stress`, `shear_flow` and `channel_shear_center`;
  `thin_walled_torsion_shear_stress` and `thin_walled_torsion_angle`;
  `thin_wall_vessel_stresses` and `thick_wall_cylinder_stresses`;
  `principal_stresses_3d`, `strain_rosette` and
  `plane_stress_from_strains`; `effective_length_factor`,
  `secant_formula_stress`, `secant_formula_load` and
  `johnson_critical_stress`; `unit_load_truss_deflection` and
  `unit_load_beam_deflection`; `slope_deflection` (a
  `SlopeDeflectionResult`); `parabolic_cable`, `catenary_cable` and
  `three_hinged_arch` (an `ArchResult`). `BeamResult` gains optional
  `support_moments` and `slope` fields and its `reactions` may hold more
  than two values. The history page gains twelve breakthroughs (Hooke
  1675, the catenary 1691, Young 1807, Lamé 1833, Jourawski 1856, Mohr
  1874, Johnson 1893, Bredt 1896, Maney 1915, Macaulay 1919, Maillart
  1921, Simmons and Ruge 1938), each with its own gallery example, and
  Mohr's circle example adds the three-dimensional circles. Existing
  calls are unchanged.

- `reliability` covers the rest of a first reliability engineering course
  (ROADMAP D9): lognormal and gamma lifetimes (`lognormal_reliability`,
  `lognormal_cdf`, `lognormal_hazard`, `lognormal_mttf`, `gamma_reliability`,
  `gamma_cdf`, `gamma_hazard`, `gamma_mttf`); Weibull parameter estimation
  from complete and censored data by median ranks (`median_ranks`,
  `adjusted_ranks`), rank regression (`weibull_rank_regression`) and
  maximum likelihood (`weibull_mle`), both returning a `WeibullFitResult`
  with the B-life and a probability plot (`plot_weibull_probability`);
  accelerated life testing (`arrhenius_acceleration_factor`,
  `arrhenius_fit`, `eyring_acceleration_factor`,
  `inverse_power_law_acceleration_factor`); minimal path and cut sets of
  general block diagrams with exact reliability and Esary-Proschan bounds
  (`minimal_path_sets`, `minimal_cut_sets`, `network_reliability`,
  `esary_proschan_bounds`); cold standby (`standby_reliability`,
  `standby_mttf`); FMEA (`risk_priority_number`, `fmea_ranking` returning
  an `FMEAResult`, `failure_mode_criticality`); Monte Carlo failure
  probability with importance sampling at the Hasofer-Lind design point
  (`monte_carlo_failure_probability`, a `MonteCarloResult`); and Bayesian
  updating of a failure rate (`bayesian_failure_rate`, a
  `FailureRatePosterior`, and `gamma_prior_from_bounds`). The history page
  gains nine breakthroughs (McAlister 1879, Arrhenius 1889, Eyring 1935,
  Metropolis and Ulam 1949, MIL-P-1629 1949, Johnson 1951, Raiffa and
  Schlaifer 1961, Esary and Proschan 1963, Cohen 1965), each with its own
  gallery example. Existing calls are unchanged.

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

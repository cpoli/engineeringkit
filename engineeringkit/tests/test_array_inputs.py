"""Array-input contract across every domain (policy in CONTRIBUTING.md).

A public *scalar function* is one whose parameters are all documented with a
scalar type (``float``, ``float or array_like``, ``complex``, ``int``,
``str``, ``bool`` or a choice set) and that returns a number, a tuple of
numbers, or a result object of numbers. Each one must broadcast its
``float`` and ``complex`` parameters over a 1-D array, unless it is listed in
``SCALAR_ONLY`` and its docstring says it is scalar-only. Parameters
documented as ``int``, ``str`` or ``bool`` (counts, orders, modes, flags)
stay scalars everywhere.

``CASES`` holds, for every scalar function, scalar calls taken from the
domain's own tests (array calls in those tests were split into their
elements); functions whose tests only call them with computed or
multi-dimensional inputs use their docstring example. A scalar function that no domain test calls yet is not required to
have a case until it has tests. The test stacks the cases into 1-D arrays, calls the function once
and compares against the scalar calls element by element.
"""

from __future__ import annotations

import dataclasses
import importlib
import inspect
import numbers
import re
from pathlib import Path

import numpy as np
import pytest

DOMAINS = (
    "aerodynamics",
    "communications",
    "control_systems",
    "electrical_engineering",
    "fluid_engineering",
    "heat_transfer",
    "materials",
    "mechanical_engineering",
    "reliability",
    "signal_processing",
    "structural_engineering",
    "thermodynamics",
)


def C(*args, **kwargs):
    return args, kwargs


CASES = {
    # aerodynamics
    "aerodynamics.aerodynamic_center": [C(6.28318530718, 0.0), C(6.0, -0.06)],
    "aerodynamics.area_mach_ratio": [C(1), C(0.5), C(2)],
    "aerodynamics.aspect_ratio": [C(10, 20)],
    "aerodynamics.best_glide_cl": [C(0.02, 8, 0.8)],
    "aerodynamics.boundary_layer_laminar": [C(1, 1000000.0)],
    "aerodynamics.breguet_range": [
        C(230, 1.6e-05, 10, 700000.0, 500000.0),
        C(230, 1.6e-05, 20, 700000.0, 500000.0),
        C(230, 1.6e-05, 10, 500000.0, 500000.0),
    ],
    "aerodynamics.breguet_range_propeller": [C(0.8, 7.6e-08, 13.0, 13100.0, 11700.0), C(0.8, 7.6e-08, 13.0, 10000.0, 10000.0)],
    "aerodynamics.climb_angle": [C(3000.0, 1500.0, 15000.0)],
    "aerodynamics.downwash_gradient": [C(4.8, 8.0)],
    "aerodynamics.drag": [C(1500.0, 2, 0.02)],
    "aerodynamics.drag_coefficient": [C(60, 1.2, 50, 2)],
    "aerodynamics.drag_polar": [C(0.634132367617, 0.02, 8, 0.8), C(0.570719130855, 0.02, 8, 0.8), C(0.697545604379, 0.02, 8, 0.8)],
    "aerodynamics.dynamic_pressure": [C(1.2, 50), C(1.2, 25.8198889747), C(1.225, 40.0)],
    "aerodynamics.fanno_choking_length": [C(0.5, 0.05, 0.02)],
    "aerodynamics.fanno_flow": [C(0.3, 1.4), C(0.8, 1.4), C(1.5, 1.4)],
    "aerodynamics.finite_wing_cl": [C(0.1, 8), C(0.1, 10.1859163579)],
    "aerodynamics.finite_wing_lift_slope": [C(10000000000.0), C(8)],
    "aerodynamics.glide_ratio": [C(0.634132367617, 0.04), C(0.570719130855, 0.0362), C(0.697545604379, 0.0442)],
    "aerodynamics.induced_drag_coefficient": [C(0.525196449085, 10.1859163579)],
    "aerodynamics.inverse_prandtl_meyer": [C(0.0)],
    "aerodynamics.isa_atmosphere": [C(0.0), C(11000.0), C(20000.0)],
    "aerodynamics.isentropic_factor": [C(1.0, 1.4)],
    "aerodynamics.kutta_joukowski_lift": [C(1.225, 40.0, 9.42477796077)],
    "aerodynamics.lift": [C(1500.0, 2, 0.5), C(400.0, 20, 1.5)],
    "aerodynamics.lift_coefficient": [C(1500, 1.2, 50, 2)],
    "aerodynamics.lifting_line": [C(8.0, 1.0, 0.0), C(8.0, 1.0, 0.1)],
    "aerodynamics.load_factor": [C(1.0471975512)],
    "aerodynamics.mach_angle": [C(1.0)],
    "aerodynamics.mach_number": [C(170, 340)],
    "aerodynamics.naca4_moment_coefficient": [C(0.02, 0.4), C(0.0, 0.4)],
    "aerodynamics.naca4_zero_lift_angle": [C(0.0, 0.4), C(0.02, 0.4)],
    "aerodynamics.neutral_point": [C(0.25, 4.9, 4.2, 0.55, 0.35, 0.9)],
    "aerodynamics.normal_shock_mach2": [C(2), C(1)],
    "aerodynamics.normal_shock_pressure_ratio": [C(1), C(1.60611226392)],
    "aerodynamics.oblique_shock": [C(2.0, 0.349065850399), C(2.5, 0.209439510239)],
    "aerodynamics.oblique_shock_angle": [C(2.0, 0.40096377004)],
    "aerodynamics.oblique_shock_deflection": [C(1.5, 0.729727656227), C(2.0, 0.523598775598), C(5.0, 0.20135792079)],
    "aerodynamics.oblique_shock_max_deflection": [C(2.0)],
    "aerodynamics.oswald_efficiency_estimate": [C(8)],
    "aerodynamics.prandtl_glauert_cp": [C(-1, 0.6)],
    "aerodynamics.prandtl_meyer_angle": [C(1.4, 1.3), C(2.5, 1.3), C(6.0, 1.3)],
    "aerodynamics.prandtl_meyer_expansion": [C(1.5, 0.349065850399)],
    "aerodynamics.rate_of_climb": [C(3000.0, 1500.0, 60.0, 15000.0)],
    "aerodynamics.rayleigh_choking_heat": [C(1.0, 300.0), C(0.2, 300.0)],
    "aerodynamics.rayleigh_flow": [C(2.0), C(1.0), C(0.2)],
    "aerodynamics.reynolds_number": [C(1, 10, 1, 1e-05)],
    "aerodynamics.skin_friction_laminar": [C(1000000.0)],
    "aerodynamics.skin_friction_turbulent": [C(100000.0)],
    "aerodynamics.speed_of_sound": [C(1.4, 287.0, 300.0), C(1.4, 287.0, 506.25), C(1.4, 287.05, 250.0)],
    "aerodynamics.stagnation_pressure": [C(100, 2)],
    "aerodynamics.stagnation_temperature": [C(300, 2)],
    "aerodynamics.stall_speed": [C(10700.0, 1.225, 16.2, 1.6), C(12000, 1.2, 20, 1.5)],
    "aerodynamics.static_margin": [C(0.561889654468, 0.3)],
    "aerodynamics.tail_volume_ratio": [C(3.0, 5.0, 16.0, 1.5)],
    "aerodynamics.takeoff_ground_roll": [C(10700.0, 1.225, 16.2, 1.6, 2400.0, 400.0)],
    "aerodynamics.thin_airfoil_cl": [C(0.1), C(0.05)],
    "aerodynamics.turn_radius": [C(60.0, 2.0)],
    "aerodynamics.turn_rate": [C(60.0, 2.0)],
    # communications
    "communications.am_output_snr": [C(1000.0)],
    "communications.aperture_gain": [C(0.0317869620723, 3000000000.0)],
    "communications.baseband_snr": [C(2e-12, 4e-21, 5000.0)],
    "communications.ber_bpsk_awgn": [C(1), C(10)],
    "communications.ber_bpsk_mrc_rayleigh": [C(10.0, 2), C(3.0, 2)],
    "communications.ber_bpsk_rayleigh": [C(10.0), C(10000.0), C(0.0)],
    "communications.ber_qpsk_awgn": [C(1)],
    "communications.binary_symmetric_capacity": [C(0.5), C(0)],
    "communications.carson_bandwidth": [C(75000.0, 15000.0)],
    "communications.db": [C(2), C(3), C(10)],
    "communications.dbm": [C(0.001)],
    "communications.delta_slope_limit": [C(0.05, 8000.0)],
    "communications.effective_aperture": [C(1.0, 3000000000.0), C(40.0, 3000000000.0), C(30.0, 3000000000.0)],
    "communications.entropy_binary": [C(0.11)],
    "communications.fm_output_snr": [C(1000.0, 1.0), C(1000.0, 5.0), C(1000.0, 2.5)],
    "communications.fm_threshold_snr": [C(5.0, 10.0)],
    "communications.free_space_path_loss": [C(200, 2000000000.0), C(100, 2000000000.0)],
    "communications.free_space_path_loss_db": [C(100, 2000000000.0)],
    "communications.friis_received_power": [C(1, 2, 3, 100, 2000000000.0), C(1.0, 30.0, 20.0, 5000.0, 3000000000.0)],
    "communications.from_db": [C(10.0)],
    "communications.from_dbm": [C(-40.6868706313)],
    "communications.half_wave_dipole_pattern": [C(1.57079632679), C(3.14159265359)],
    "communications.half_wave_dipole_resistance": [C(376.730313668), C(377.0)],
    "communications.hartley_capacity": [C(100, 4)],
    "communications.link_budget_db": [C(30, 3.01029995664, 4.7712125472, 78.4683831352)],
    "communications.matched_filter_snr": [C(55.2971845476, 0.3)],
    "communications.mu_law_compress": [C(3.0, 100, 3), C(0.03, 100, 3)],
    "communications.noise_figure": [C(10)],
    "communications.nyquist_symbol_rate": [C(100)],
    "communications.q_function": [C(5.41196100146), C(0)],
    "communications.raised_cosine_bandwidth": [C(1000.0, 0.5)],
    "communications.raised_cosine_pulse": [C(0.002, 0.001, 0.25), C(0.002000000001, 0.001, 0.25), C(0.3, 1.0, 0.35)],
    "communications.raised_cosine_spectrum": [C(0.1, 1.0, 0.5), C(0.5, 1.0, 0.0)],
    "communications.root_raised_cosine_pulse": [C(0.5, 1.0, 0.5), C(0.5000001, 1.0, 0.5), C(0.0, 1.0, 0.5)],
    "communications.selection_outage_probability": [C(2.0, 10.0, 2)],
    "communications.ser_mpsk_awgn": [C(19.9526231497, 8), C(0.0, 8), C(100.0, 8)],
    "communications.ser_mqam_awgn": [C(19.9526231497, 16), C(0.0, 16)],
    "communications.shannon_capacity": [C(1000000.0, 0), C(1000000.0, 1), C(1, 1)],
    "communications.spectral_efficiency": [C(400, 100)],
    "communications.symbol_rate": [C(400, 4)],
    "communications.thermal_noise_power": [C(1, 1)],
    "communications.wavelength": [C(3000000000.0), C(100000000.0), C(299792458)],
    "communications.mu_law_expand": [C(0.3), C(-0.7)],
    "communications.rayleigh_pdf": [C(0.7071067811865476)],
    # control_systems
    "control_systems.dead_zone_describing_function": [C(3.08121358871, 1.0), C(3.13502206257, 1.0), C(3.12811466216, 1.0)],
    "control_systems.nichols_closed_loop": [C(0.0, -180.0)],
    "control_systems.relay_describing_function": [C(0.156851516437, 1.0), C(0.15959067764, 1.0), C(0.159161232635, 1.0)],
    "control_systems.saturation_describing_function": [C(1.99871963866, 1.0), C(2.03362408469, 1.0), C(2.03309868143, 1.0)],
    "control_systems.second_order_metrics": [C(0.5, 2)],
    "control_systems.ziegler_nichols": [C(8.0, 4.0)],
    # electrical_engineering
    "electrical_engineering.capacitors_parallel": [C(2, 3)],
    "electrical_engineering.capacitors_series": [C(2, 2)],
    "electrical_engineering.change_base": [C(0.08j, 50000000.0, 13200.0, 100000000.0, 13800.0)],
    "electrical_engineering.critical_clearing_time": [C(4.0, 0.8, 1.6, frequency=60.0, p_max_post=1.6)],
    "electrical_engineering.dc_motor_speed": [C(240.0, 0.5, 2.0, 960.0), C(240.0, 0.5, 2.0, 760.0), C(240.0, 0.5, 2.0, 560.0)],
    "electrical_engineering.dc_motor_torque": [C(240.0, 0.5, 2.0, 120.0), C(240.0, 0.5, 2.0, 0.0)],
    "electrical_engineering.diode_operating_point": [C(3.0, 220.0, 2e-09, ideality=1.8), C(5.0, 1000.0, 1e-14, ideality=1.5)],
    "electrical_engineering.distortionless_condition": [C(0.1, 2.5e-07, 4e-05, 1e-10), C(0.1, 2.5e-07, 4e-06, 1e-10)],
    "electrical_engineering.equal_area_clearing_angle": [C(1.0, 2.0, 0.0, 2.0), C(1.0, 2.0, 1.9, 2.0)],
    "electrical_engineering.fault_currents": [C(1.0, 0.25j, 0.35j, 0.1j, "3ph")],
    "electrical_engineering.induction_motor_breakdown": [
        C(230.0, 157.079632679, 0.3, 0.6, 0.25, 0.6, 20.0),
        C(230.0, 157.079632679, 0.3, 0.6, 0.5, 0.6, 20.0),
    ],
    "electrical_engineering.induction_motor_torque": [C(230.0, 0.05, 157.079632679, 0.0, 1e-09, 0.25, 0.6, 1000000000.0)],
    "electrical_engineering.input_impedance": [C((30 + 40j), 50.0, 1.57079632679), C((30 + 40j), 50.0, 3.14159265359), C(0, 50.0, 0.3)],
    "electrical_engineering.inverting_amplifier_gain": [C(47000.0, 4700.0)],
    "electrical_engineering.maximum_power_transfer": [C(10, 5)],
    "electrical_engineering.noninverting_amplifier_gain": [C(0.0, 1000.0)],
    "electrical_engineering.norton_from_thevenin": [C(10, 5)],
    "electrical_engineering.ohms_law": [C(current=2, resistance=5)],
    "electrical_engineering.parallel_resistance": [C(2, 2)],
    "electrical_engineering.per_unit_base": [C(100000000.0, 138000.0), C(50000000.0, 13200.0), C(100000000.0, 13800.0)],
    "electrical_engineering.phase_components": [C(1, 0, 0), C(0, 1, 0)],
    "electrical_engineering.phasor": [C(2.0, 0.5)],
    "electrical_engineering.polar": [C((1.75516512378 + 0.958851077208j))],
    "electrical_engineering.power_factor": [C(8, 10)],
    "electrical_engineering.propagation_constant": [
        C(0, 2.5e-07, 0, 1e-10, 1000.0),
        C(0, 2.5e-07, 0, 1e-10, 1000000.0),
        C(0, 2.5e-07, 0, 1e-10, 1000000000.0),
    ],
    "electrical_engineering.rc_charging_voltage": [C(5, 1000, 0.001, 100)],
    "electrical_engineering.reactive_power": [C(10, 2, 1.57079632679)],
    "electrical_engineering.real_power": [C(10, 2, 0.8)],
    "electrical_engineering.reflection_coefficient": [C(50.0, 50.0), C(0, 50.0), C(150.0, 50.0)],
    "electrical_engineering.resonance_frequency": [C(0.1, 0.0001), C(1, 1)],
    "electrical_engineering.rl_current": [C(10.0, 5.0, 0.5, 0.1)],
    "electrical_engineering.rl_time_constant": [C(2, 4), C(0.5, 5.0)],
    "electrical_engineering.rlc_step_voltage": [C(1.0, 0.0, 1.0, 1.0, 3.14159265359)],
    "electrical_engineering.rms_sine": [C(1.41421356237)],
    "electrical_engineering.series_resistance": [C(2, 3)],
    "electrical_engineering.series_rlc_impedance": [C(10, 0.1, 0.0001, 316.227766017)],
    "electrical_engineering.shockley_diode_current": [C(-1.0, 1e-14)],
    "electrical_engineering.standing_wave_ratio": [C((0.5 + 0j))],
    "electrical_engineering.swing_response": [
        C(4.0, 0.8, 1.6, 0.0, 1.6, 0.207794518199, 2.0, frequency=60.0),
        C(4.0, 0.8, 1.6, 0.0, 1.6, 0.220647787366, 2.0, frequency=60.0),
        C(4.0, 0.8, 1.6, 0.0, 1.6, 0.171376922226, 2.0, frequency=60.0),
    ],
    "electrical_engineering.symmetrical_components": [C(1, (-0.5 - 0.866025403784j), (-0.5 + 0.866025403784j))],
    "electrical_engineering.thermal_voltage": [C(300.0)],
    "electrical_engineering.thevenin_from_source": [C(10, 5)],
    "electrical_engineering.three_phase_power": [C(10, 2, 0.8)],
    "electrical_engineering.transformer_secondary": [C(240, 100, 50)],
    "electrical_engineering.transmission_line_velocity": [C(2e-07, 8e-11), C(2.5e-07, 1e-10)],
    "electrical_engineering.transmission_line_z0": [C(2e-07, 8e-11)],
    "electrical_engineering.wheatstone_balanced": [C(1, 2, 2, 4), C(1, 2, 3, 4)],
    "electrical_engineering.rc_time_constant": [C(1e3, 1e-6)],
    "electrical_engineering.capacitor_impedance": [C(1e-3, 1e3)],
    "electrical_engineering.inductor_impedance": [C(1e-3, 1e3)],
    "electrical_engineering.rlc_damping": [C(2.0, 1.0, 1.0)],
    # fluid_engineering
    "fluid_engineering.affinity_laws": [
        C(0.05, 32.0, 21000.0, 0.5, 1.0),
        C(0.05, 32.0, 21000.0, 1.2, 0.9),
        C(0.05, 32.0, 21000.0, 1.0, 1.1),
    ],
    "fluid_engineering.alternate_depths": [C(2.0, 2.5), C(1.29085572827, 2.5)],
    "fluid_engineering.bend_force": [
        C(998.0, 0.08, 0.0314, 0.0154, 200000.0, 189773.062131, 0.0),
        C(998.0, 0.08, 0.0314, 0.0154, 200000.0, 189773.062131, 0.785398163397),
        C(998.0, 0.08, 0.0314, 0.0154, 200000.0, 189773.062131, 1.57079632679),
    ],
    "fluid_engineering.bernoulli_pressure": [C(100000.0, 1000.0, 1.0, 2.0), C(200000.0, 998.0, 2.54777070064, 5.19480519481)],
    "fluid_engineering.buoyancy_force": [C(1000, 0.01)],
    "fluid_engineering.cavitation_number": [C(3000, 1000, 1000, 2), C(90000.0, 2300.0, 998.0, 3.0)],
    "fluid_engineering.center_of_pressure_depth": [C(3.5, 4.5, 3.375), C(1.5, 3.0, 2.25)],
    "fluid_engineering.colebrook_friction_factor": [C(4000, 0), C(100000.0, 0.0001), C(100000000.0, 0.05)],
    "fluid_engineering.continuity_velocity": [C(0.02, 1.0, 0.01)],
    "fluid_engineering.couette_flow_rate": [C(0.002, 1.5, 0.08, 0.0), C(0.002, 1.5, 0.08, -20000.0), C(0.002, 1.5, 0.08, 30000.0)],
    "fluid_engineering.couette_shear_stress": [C(0.1, 2.0, 0.001), C(0.9, 0.3, 0.0001)],
    "fluid_engineering.couette_velocity": [
        C(0.001, 0.002, 1.5, 0.08, 0.0),
        C(2.60934714828e-05, 0.002, 1.5, 0.08, 0.0),
        C(0.00197390652852, 0.002, 1.5, 0.08, 0.0),
    ],
    "fluid_engineering.couette_viscometer_torque": [C(0.9, 3.0, 0.02, 0.025, 0.1), C(0.9, 3.0, 0.1, 0.1001, 0.1)],
    "fluid_engineering.critical_closure_time": [C(600.0, 1200.0)],
    "fluid_engineering.critical_depth": [C(3.0), C(2.5)],
    "fluid_engineering.darcy_head_loss": [C(0.32, 1.0, 0.004, 0.05), C(0.021, 250.0, 0.15, 1.69765272631)],
    "fluid_engineering.dynamic_pressure": [C(1000, 2)],
    "fluid_engineering.euler_number": [C(4000, 1000, 2)],
    "fluid_engineering.euler_turbomachine_head": [C(22.5, 16.0, 7.5, 2.0)],
    "fluid_engineering.froude_number": [
        C(3.08707591961, 0.971793398711),
        C(1.30690514528, 1.91291618143),
        C(5.50654736291, 0.454004993554),
    ],
    "fluid_engineering.hydraulic_diameter": [C(3.14159265359, 6.28318530718), C(0.282743338823, 1.88495559215)],
    "fluid_engineering.hydraulic_jump_head_loss": [C(0.25, 1.68572275115)],
    "fluid_engineering.hydraulic_radius": [C(1.5, 4.0), C(3.0, 5.0), C(6.0, 7.0)],
    "fluid_engineering.hydrostatic_force_plane": [C(1000.0, 3.5, 4.5)],
    "fluid_engineering.hydrostatic_pressure": [C(1000, 1), C(1000.0, 3.5), C(1000.0, 2.03914020722)],
    "fluid_engineering.jet_force": [C(1000.0, 0.002, 30.0, 1.57079632679), C(1000.0, 0.002, 30.0, 3.14159265359)],
    "fluid_engineering.laminar_friction_factor": [C(200.0), C(1000)],
    "fluid_engineering.mach_number": [C(170, 340)],
    "fluid_engineering.manning_flow": [C(0.015, 1.5, 4.0, 0.0004), C(0.015, 3.0, 5.0, 0.0004), C(0.015, 6.0, 7.0, 0.0004)],
    "fluid_engineering.manning_velocity": [C(0.015, 0.375, 0.0004), C(0.015, 0.6, 0.0004), C(0.015, 0.857142857143, 0.0004)],
    "fluid_engineering.minor_head_loss": [C(2, 3)],
    "fluid_engineering.nikuradse_friction_factor": [C(0.001), C(0.01), C(0.033)],
    "fluid_engineering.npsh_available": [C(90000.0, 0.0, 2300.0, 998.0), C(90000.0, 3.0, 2300.0, 998.0)],
    "fluid_engineering.orifice_flow": [C(0.003, 3.0, 0.61), C(0.003, 2.97, 0.61), C(0.003, 2.99998121769, 0.61)],
    "fluid_engineering.pelton_wheel_power": [
        C(1000.0, 0.05, 60.0, 22.917960675, 2.87979326579),
        C(1000.0, 0.05, 60.0, 37.082039325, 2.87979326579),
        C(1000.0, 0.05, 60.0, 45.83592135, 2.87979326579),
    ],
    "fluid_engineering.pipe_resistance": [C(0.021, 250.0, 0.15)],
    "fluid_engineering.pipe_wave_speed": [C(2200000000.0, 998.0, 0.5, 0.01, 200000000000.0, c1=0)],
    "fluid_engineering.pitot_velocity": [C(2000, 1000)],
    "fluid_engineering.poiseuille_flow": [C(100.0, 0.002, 1.0, 0.001), C(2, 1, 1, 1), C(1, 1, 1, 1)],
    "fluid_engineering.poiseuille_velocity": [
        C(0.0015, 50.0, 0.003, 2.0, 0.001),
        C(3.91402072242e-05, 50.0, 0.003, 2.0, 0.001),
        C(0.00296085979278, 50.0, 0.003, 2.0, 0.001),
    ],
    "fluid_engineering.pressure_head": [C(9806.65, 1000)],
    "fluid_engineering.pump_curve_head": [C(0.0702834972619, 45.0, 2500.0, 30.0), C(0.0, 45.0, 2500.0, 30.0), C(1.0, 45.0, 2500.0, 30.0)],
    "fluid_engineering.pump_hydraulic_power": [C(1000, 0.01, 10), C(1000.0, 0.05, 32.0)],
    "fluid_engineering.pump_input_power": [C(1000, 0.01, 10, 0.8)],
    "fluid_engineering.pump_operating_point": [
        C(30.0, 3000.0, 5.0, 20000.0, parallel=2),
        C(30.0, 3000.0, 5.0, 200.0, parallel=2),
        C(40.0, 1000.0, 10.0, 2000.0, parallel=2),
    ],
    "fluid_engineering.pump_specific_speed": [C(10, 4, 1)],
    "fluid_engineering.rectangular_weir_flow": [C(0.7, 1.8), C(0.35, 1.8)],
    "fluid_engineering.rehbock_discharge_coefficient": [C(0.0, 1.0), C(0.4, 0.8)],
    "fluid_engineering.reynolds_number": [
        C(1000, 0.05, 0.004, 0.001),
        C(998.0, 0.00225008136111, 5e-05, 0.001),
        C(998.0, 0.0245989109215, 0.0002, 0.001),
    ],
    "fluid_engineering.sequent_depth": [C(0.25, 8.0), C(1.68572275115, 1.18643471985), C(1.0, 3.13155712067)],
    "fluid_engineering.smooth_pipe_friction_factor": [C(100000.0)],
    "fluid_engineering.specific_energy": [C(1.97163345513, 3.0), C(3.12836654487, 3.0), C(1.25673308975, 3.0)],
    "fluid_engineering.sphere_drag_coefficient": [C(4.90994261993), C(144.978082391), C(2244.89785555)],
    "fluid_engineering.stokes_drag": [C(0.001, 0.00225008136111, 5e-05), C(0.001, 0.00100200400802, 5e-05)],
    "fluid_engineering.stokes_settling_velocity": [
        C(5e-05, 2650.0, 998.0, 0.001),
        C(0.001, 50.0, 998.0, 0.001),
        C(5e-06, 2650.0, 998.0, 0.001),
    ],
    "fluid_engineering.swamee_jain_friction_factor": [C(100000.0, 0.0001)],
    "fluid_engineering.system_curve_head": [C(0.0702834972619, 15.0, 4000.0), C(0.0, 15.0, 4000.0), C(1.0, 15.0, 4000.0)],
    "fluid_engineering.tank_drain_time": [C(2.0, 0.003, 3.0, 0.75, 0.61), C(2.0, 0.003, 3.0, 1.5, 0.61), C(2.0, 0.003, 3.0, 0.0, 0.61)],
    "fluid_engineering.terminal_velocity": [
        C(0.0002, 2650.0, 998.0, 0.001),
        C(0.001, 2650.0, 998.0, 0.001),
        C(0.005, 2650.0, 998.0, 0.001),
    ],
    "fluid_engineering.thoma_sigma": [C(8.96083285352, 40.0)],
    "fluid_engineering.v_notch_weir_flow": [C(0.7), C(0.35)],
    "fluid_engineering.velocity_head": [C(4.42869055139), C(3.0)],
    "fluid_engineering.venturi_flow": [C(1500.0, 1000.0, 0.02, 0.01), C(20000.0, 998.0, 0.0177, 0.0044)],
    "fluid_engineering.water_hammer_pressure": [C(1000, 1200, 1)],
    "fluid_engineering.weber_number": [C(1000, 1, 0.01, 0.1)],
    # heat_transfer
    "heat_transfer.biot_number": [C(10, 0.01, 1)],
    "heat_transfer.blackbody_peak_wavelength": [C(1000)],
    "heat_transfer.chilton_colburn_mass_transfer": [C(9.61009355872, 1.16, 1018.80877743, 0.846153846154)],
    "heat_transfer.churchill_bernstein_nusselt": [C(0.0, 0.7), C(6071.0, 0.7)],
    "heat_transfer.churchill_chu_nusselt": [C(0.0, 0.7)],
    "heat_transfer.convection_heat_rate": [C(10, 2, 400, 300)],
    "heat_transfer.convection_resistance": [C(10, 2)],
    "heat_transfer.critical_insulation_radius_cylinder": [C(2, 10)],
    "heat_transfer.critical_insulation_radius_sphere": [C(2, 10)],
    "heat_transfer.cylinder_resistance": [C(1.0, 1.000001, 20.0, 2.0), C(0.05, 0.15, 1.0, 2.0)],
    "heat_transfer.cylinder_transient": [C(0.15, 2.0, 1.0), C(0.00391402072242, 2.0, 1.0), C(0.296085979278, 2.0, 1.0)],
    "heat_transfer.dittus_boelter_nusselt": [C(10000.0, 1.0)],
    "heat_transfer.entry_length": [C(1000.0, 0.01)],
    "heat_transfer.exchanger_effectiveness_counterflow": [C(0, 0), C(1e-12, 0), C(0, 0.5)],
    "heat_transfer.exchanger_effectiveness_crossflow": [C(0.25, 0.25), C(0.25, 0.5), C(0.25, 0.75)],
    "heat_transfer.exchanger_effectiveness_parallel": [C(2, 0.5)],
    "heat_transfer.exchanger_effectiveness_shell_tube": [C(1.5, 0.3, 1), C(1.5, 0.6, 1), C(1.5, 1.0, 1)],
    "heat_transfer.fin_efficiency_adiabatic": [C(0, 1), C(1, 0), C(1e-08, 1)],
    "heat_transfer.fin_heat_rate_adiabatic": [
        C(1, 1, 1, 1, 1, 400, 300),
        C(25.0, 0.04, 200.0, 5e-05, 0.1, 400.0, 300.0),
        C(25.0, 0.04, 200.0, 5e-05, 0.10125, 400.0, 300.0),
    ],
    "heat_transfer.fin_heat_rate_convective_tip": [
        C(25.0, 0.04, 200.0, 5e-05, 0.0, 400.0, 300.0),
        C(25.0, 0.04, 200.0, 5e-05, 0.1, 400.0, 300.0),
    ],
    "heat_transfer.fin_parameter": [C(1, 1, 1, 1)],
    "heat_transfer.flat_plate_nusselt": [C(1.0, 0.6, local=True), C(1.0, 0.7, local=True), C(1.0, 1.0, local=True)],
    "heat_transfer.fourier_heat_rate": [C(10, 2, 20, 0.5)],
    "heat_transfer.fourier_number": [C(0.01, 2, 0.1)],
    "heat_transfer.hausen_nusselt": [C(1000.0, 5.0, 0.01, 1000000000.0), C(1000.0, 5.0, 0.01, 0.1), C(1000.0, 5.0, 0.01, 1.0)],
    "heat_transfer.heat_exchanger_rate": [C(10, 2, 20)],
    "heat_transfer.lewis_number": [C(2.0, 4.0), C(2.2e-05, 2.6e-05)],
    "heat_transfer.lmtd_correction_factor": [C(0.717302495193, 0.3, 2), C(0.656708287928, 0.6, 2), C(0.578695223296, 1.0, 2)],
    "heat_transfer.log_mean_temperature_difference": [C(20, 20), C(20.00000001, 20), C(0.309227144821, 0.792768143446)],
    "heat_transfer.lumped_temperature": [
        C(400, 300, 10, 2, 1000, 1000, 0.01, 0),
        C(400, 300, 10, 2, 1000, 1000, 0.01, 1000000.0),
        C(400, 300, 10, 1, 1000, 1000, 0.01, 0),
    ],
    "heat_transfer.neumann_stefan_front": [C(1.0, 1.0, 0.5), C(0.0, 1.0, 0.5)],
    "heat_transfer.nusselt_film_condensation": [
        C(960.0, 0.6, 0.0003, 0.68, 2260000.0, 373.0, 363.0, 0.2),
        C(960.0, 0.6, 0.0003, 0.68, 2260000.0, 373.0, 363.0, 3.2),
    ],
    "heat_transfer.overall_u": [C(0.1, 0.1)],
    "heat_transfer.planck_spectral_emissive_power": [C(1e-09, 10.0), C(1.0, 300.0)],
    "heat_transfer.plane_wall_resistance": [C(1e-06, 20.0, 12.5663706144)],
    "heat_transfer.plane_wall_transient": [C(0.2, 2.0, 1.0), C(0.00521869429657, 2.0, 1.0), C(0.394781305703, 2.0, 1.0)],
    "heat_transfer.radiation_coefficient": [C(0.8, 400, 300)],
    "heat_transfer.radiation_heat_rate": [C(1, 1, 300, 300), C(0, 1, 400, 300), C(0.8, 2, 400, 300)],
    "heat_transfer.rayleigh_number": [C(1.0, 1.0, 1.0, 1.0, 1.0, g=2.0)],
    "heat_transfer.rohsenow_boiling_heat_flux": [C(0.000279, 2257000.0, 957.9, 0.5955, 0.0589, 4217.0, 1.76, 10.0, C_sf=0.0065)],
    "heat_transfer.semi_infinite_temperature": [C(400, 300, 0, 0.01, 1)],
    "heat_transfer.shape_factor_buried_cylinder": [C(0.3, 0.6, 2.0)],
    "heat_transfer.shape_factor_buried_sphere": [C(1.0, 1000000000.0)],
    "heat_transfer.shape_factor_eccentric_cylinders": [C(0.1, 0.3, 0.0, 2.0)],
    "heat_transfer.shape_factor_parallel_cylinders": [C(0.3, 0.3, 1.2, 2.0)],
    "heat_transfer.sphere_resistance": [C(1, 2, 1)],
    "heat_transfer.sphere_transient": [C(0.15, 2.0, 1.0), C(0.00391402072242, 2.0, 1.0), C(0.296085979278, 2.0, 1.0)],
    "heat_transfer.stefan_boltzmann_emissive_power": [C(300), C(800.0)],
    "heat_transfer.stefan_number": [C(1.0, 2.0, 4.0)],
    "heat_transfer.tube_bank_max_velocity": [C(1.0, 0.01, 0.03, 0.005, "staggered"), C(1.0, 0.01, 0.03, 0.03, "staggered")],
    "heat_transfer.two_surface_radiation": [C(500.0, 300.0, 0.8, 0.6, 1.0, 4.0), C(500.0, 300.0, 0.3, 0.3, 1.0, 1.0)],
    "heat_transfer.view_factor_coaxial_disks": [C(1, 1, 1e-06), C(0.1, 0.1, 100.0), C(0.5, 1.0, 0.4)],
    "heat_transfer.view_factor_parallel_plates": [C(1, 1, 0.0001), C(2.0, 1.0, 0.5), C(1.0, 2.0, 0.5)],
    "heat_transfer.view_factor_reciprocity": [C(0.831935614725, 0.785398163397, 3.14159265359)],
    "heat_transfer.zuber_critical_heat_flux": [C(957.9, 0.5955, 2257000.0, 0.0589)],
    "heat_transfer.zukauskas_nusselt": [C(13943.0, 0.71, 0.0313, 0.0343, "staggered", rows=7, Pr_s=0.7)],
    # materials
    "materials.arrhenius_creep_rate": [C(1, 1, 1, 1, 2)],
    "materials.avrami_fraction": [C(1200.0, 1.01580591857e-10, 3.1)],
    "materials.avrami_time": [C(0.5, 1.01580591857e-10, 3.1)],
    "materials.basquin_cycles": [C(200, 1000, -0.1), C(300, 1000, -0.1), C(429621159.08, 1000000000.0, -0.09)],
    "materials.bragg_angle": [C(1.01328401744e-10, 1.79e-10)],
    "materials.brinell_hardness": [C(4903.325, 0.01, 0.00162)],
    "materials.brinell_tensile_strength": [C(240.975133378)],
    "materials.bulk_modulus": [C(200000000000.0, 0.3)],
    "materials.coffin_manson_cycles": [C(0.0021481057954, 0.6, -0.6)],
    "materials.constrained_thermal_stress": [C(100, 1e-05, 100)],
    "materials.critical_crack_length": [C(50000000.0, 100000000.0, 1.12)],
    "materials.cubic_lattice_parameter": [C(1.28e-10, "fcc")],
    "materials.density_rule_of_mixtures": [C(100, 200, 0.5)],
    "materials.diffusion_coefficient": [C(0.00012, 130000.0, 823.0)],
    "materials.diffusion_time": [C(0.0005, 0.8, 1.6e-11, 0.25, 1.2), C(0.00011, 0.745476568116, 1e-11, 0.2, 1.0)],
    "materials.elliptical_hole_concentration": [C(1.0, 1.0), C(4.0, 1.0)],
    "materials.endurance_limit_estimate": [C(600000000.0), C(2000000000.0)],
    "materials.factor_of_safety": [C(100, 20)],
    "materials.fatigue_notch_factor": [C(2.5, 0.285714285714), C(2.5, 0.8), C(2.5, 0.975609756098)],
    "materials.fick_flux": [C(3e-11, -80.0)],
    "materials.gerber_safety_factor": [C(0.0, 100000000.0, 200000000.0, 600000000.0)],
    "materials.goodman_safety_factor": [C(10, 20, 100, 200)],
    "materials.griffith_fracture_stress": [C(3.14159265359, 2, 1), C(70000000000.0, 1.0, 1e-06)],
    "materials.hall_petch": [C(100, 2, 4)],
    "materials.halpin_tsai": [
        C(72000000000.0, 3500000000.0, 0.0, 0.0),
        C(72000000000.0, 3500000000.0, 0.2, 0.0),
        C(72000000000.0, 3500000000.0, 0.4, 0.0),
    ],
    "materials.hooke_stress": [C(100, 0.01), C(100, 0.02)],
    "materials.ideal_shear_strength": [C(26000000000.0, 2.86e-10, 2.86e-10)],
    "materials.interplanar_spacing": [C(2.866e-10, 2, 2, 0)],
    "materials.lame_parameters": [C(100, 0.25)],
    "materials.larson_miller_parameter": [C(900.0, 3600000.0), C(950.0, 221705.355984)],
    "materials.larson_miller_rupture_time": [C(20700.0, 900.0), C(20700.0, 950.0)],
    "materials.lever_rule": [C(35.0, 31.5, 42.5)],
    "materials.marin_endurance_limit": [C(300000000.0, ka=0.8, kb=0.9, kc=0.85, ke=0.814)],
    "materials.marin_reliability_factor": [C(0.99), C(0.5)],
    "materials.marin_size_factor": [C(0.032)],
    "materials.marin_surface_factor": [C(520000000.0, "machined")],
    "materials.neuber_notch_sensitivity": [C(0.0001, 0.0001), C(0.001, 0.0001), C(0.01, 0.0001)],
    "materials.neuber_notch_stress": [
        C(-300000000.0, 3.0, 200000000000.0, 400000000.0, 8.0),
        C(0.0, 3.0, 200000000000.0, 400000000.0, 8.0),
        C(50000000.0, 3.0, 200000000000.0, 400000000.0, 8.0),
    ],
    "materials.normal_strain": [C(0.01, 1)],
    "materials.normal_stress": [C(100, 2)],
    "materials.paris_crack_growth": [C(2, 3, 2)],
    "materials.peierls_stress": [
        C(26000000000.0, 0.35, 2.86e-10, 2.288e-10),
        C(26000000000.0, 0.35, 2.86e-10, 2.86e-10),
        C(26000000000.0, 0.35, 2.86e-10, 3.432e-10),
    ],
    "materials.peterson_notch_sensitivity": [C(0.0001, 0.00025), C(0.001, 0.00025), C(0.01, 0.00025)],
    "materials.ramberg_osgood_strain": [C(300000000.0, 70000000000.0, 300000000.0, 15), C(-300000000.0, 70000000000.0, 300000000.0, 15)],
    "materials.resilience_modulus": [C(200000000.0, 200000000000.0)],
    "materials.rule_of_mixtures_E_longitudinal": [
        C(200000000000.0, 3000000000.0, 0),
        C(200000000000.0, 3000000000.0, 0.25),
        C(200000000000.0, 3000000000.0, 0.5),
    ],
    "materials.rule_of_mixtures_E_transverse": [
        C(200000000000.0, 3000000000.0, 0),
        C(200000000000.0, 3000000000.0, 0.25),
        C(200000000000.0, 3000000000.0, 0.5),
    ],
    "materials.shear_modulus": [C(200000000000.0, 0.3), C(200000000000.0, -0.2), C(210000000000.0, 0.3)],
    "materials.soderberg_safety_factor": [
        C(80000000.0, 0.0, 200000000.0, 450000000.0),
        C(80000000.0, 100000000.0, 200000000.0, 450000000.0),
        C(80000000.0, 300000000.0, 200000000.0, 450000000.0),
    ],
    "materials.strain_life_amplitude": [C(5967.54876998, 200000000000.0, 1000000000.0, -0.09, 0.6, -0.6)],
    "materials.stress_intensity": [C(1.12, 100000000.0, 0.0634386731074), C(1, 100000000.0, 0.0795774715459), C(1.0, 211100412.282, 1e-06)],
    "materials.theoretical_density": [C(4, 0.0635, 4.74531328121e-29)],
    "materials.thermal_strain": [C(1e-05, 100)],
    "materials.vickers_hardness": [C(9.80665, 0.001)],
    "materials.von_mises_plane_stress": [C(100, 100, 0), C(0, 0, 100), C(80.0, 20.0, 40.0)],
    "materials.von_mises_stress": [C(50.0, -20.0, 10.0, 30.0, -15.0, 5.0)],
    "materials.weibull_failure_probability": [
        C(400000000.0, 300000000.0, 10.0, 1e-06, 1e-06),
        C(200474893.451, 300000000.0, 10.0, 0.001, 1e-06),
    ],
    "materials.weibull_mean_strength": [C(300000000.0, 10.0)],
    "materials.weibull_size_effect": [C(400000000.0, 1e-06, 0.001, 10.0)],
    "materials.crack_growth_life": [C(1e-29, 3.0, 100e6, 1e-3, 1e-2)],
    "materials.semi_infinite_concentration": [C(0.0, 1.0, 1e-11, 0.2, 1.0), C(1e-6, 1.0, 1e-11, 0.2, 1.0)],
    # mechanical_engineering
    "mechanical_engineering.angular_momentum": [C(2, 3)],
    "mechanical_engineering.barth_velocity_factor": [C(9.0, "hobbed"), C(0.0, "hobbed")],
    "mechanical_engineering.bearing_l10_life": [C(2, 1)],
    "mechanical_engineering.belt_speed": [C(1, 60)],
    "mechanical_engineering.belt_tension_ratio": [C(0.0, 5.0), C(0.25, 6.28318530718)],
    "mechanical_engineering.bergstrasser_factor": [C(8.0)],
    "mechanical_engineering.bolt_stiffness": [C(207000000000.0, 0.0001, 0.02, 8e-05, 0.01)],
    "mechanical_engineering.bolted_joint": [
        C(800000000.0, 3200000000.0, 25000.0, 10000.0),
        C(800000000.0, 3200000000.0, 25000.0, 15000.0),
        C(800000000.0, 3200000000.0, 25000.0, 62500.0),
    ],
    "mechanical_engineering.cam_follower_motion": [C(0.5, 1.0, 1.0, "uniform")],
    "mechanical_engineering.centripetal_force": [C(2, 3, 1)],
    "mechanical_engineering.damped_frequency": [C(4, 1, 4)],
    "mechanical_engineering.disk_inertia": [C(2, 1)],
    "mechanical_engineering.flywheel_energy": [C(2, 3, 1)],
    "mechanical_engineering.fourbar_kinematics": [C(1.0, 3.0, 2.0, 3.0, 0.4, 3.0)],
    "mechanical_engineering.fourbar_output_angle": [C(1.0, 3.0, 2.0, 3.0, 1.57079632679), C(1.0, 3.0, 2.0, 3.0, 0.3)],
    "mechanical_engineering.gear_output_speed": [C(100.0, 20, 60)],
    "mechanical_engineering.gear_output_torque": [C(10.0, 20, 60, 0.9)],
    "mechanical_engineering.jeffcott_response": [C(10.0, 1000000.0, 0.0001, 316.227766017)],
    "mechanical_engineering.kinetic_energy": [C(2, 3)],
    "mechanical_engineering.lewis_bending_stress": [C(2000.0, 0.025, 0.004, 0.244502314815)],
    "mechanical_engineering.lewis_form_factor": [C(12), C(400), C(23)],
    "mechanical_engineering.long_bearing_eccentricity": [C(0.0265686661514)],
    "mechanical_engineering.magnification_factor": [C(1.3, 0.1), C(0.0, 0.3), C(1.0, 0.0)],
    "mechanical_engineering.member_stiffness": [C(207000000000.0, 0.01, 0.03)],
    "mechanical_engineering.natural_frequency": [C(100, 4)],
    "mechanical_engineering.natural_frequency_hz": [C(4, 1)],
    "mechanical_engineering.parallel_axis": [C(1.0, 2, 1)],
    "mechanical_engineering.petroff_friction_torque": [C(0.05, 0.02, 0.04, 2e-05, 100.0), C(0.05, 0.02, 0.04, 2e-05, -100.0)],
    "mechanical_engineering.phase_angle": [C(1.3, 0.1), C(100.0, 0.1)],
    "mechanical_engineering.power_torque": [C(27.0, 33.3333333333), C(10.0, 100.0)],
    "mechanical_engineering.preload_torque": [C(25000.0, 0.012, 0.18)],
    "mechanical_engineering.rotating_unbalance_force": [C(2, 0.1, 10)],
    "mechanical_engineering.rotational_energy": [C(2, 3)],
    "mechanical_engineering.screw_torque_raise": [C(100, 0.02, 0.01, 0)],
    "mechanical_engineering.shaft_diameter": [C(420.0, 610.0, 400000000.0, 2.0)],
    "mechanical_engineering.shaft_diameter_goodman": [
        C(180000000.0, 620000000.0, safety_factor=2.0, Kf=1.6, Kfs=1.3, Ma=300.0, Mm=50.0, Ta=40.0, Tm=500.0)
    ],
    "mechanical_engineering.shaft_max_shear_stress": [C(0.035, 420.0, 610.0), C(0.033536634313, 420.0, 610.0)],
    "mechanical_engineering.shaft_von_mises_stress": [
        C(0.035, 420.0, 610.0),
        C(0.0325137798448, 420.0, 610.0),
        C(0.0418401895586, 480.0, 52.0),
    ],
    "mechanical_engineering.slider_crank_displacement": [C(1, 3, 0.0), C(1, 3, 3.14159265359)],
    "mechanical_engineering.solid_sphere_inertia": [C(5, 1)],
    "mechanical_engineering.sommerfeld_number": [C(0.04, 250.0, 2079973.93257, 0.025, 1.0, 3e-05)],
    "mechanical_engineering.spring_index": [C(0.032, 0.004)],
    "mechanical_engineering.spring_rate": [C(0.004, 0.032, 8.0, 79300000000.0)],
    "mechanical_engineering.spring_shear_stress": [C(150.0, 0.032, 0.004, "none")],
    "mechanical_engineering.transmissibility": [C(0, 0.1), C(1, 0)],
    "mechanical_engineering.tuned_mass_damper": [C(10.0, 10000.0, 0.2, damping_ratio=0.05), C(10.0, 10000.0, 0.2, damping_ratio=0.3)],
    "mechanical_engineering.wahl_factor": [C(8.0), C(1000000.0)],
    "mechanical_engineering.work_constant_force": [C(10, 2, 1.0471975512)],
    "mechanical_engineering.damping_ratio": [C(4.0, 4.0, 1.0)],
    "mechanical_engineering.slider_crank_kinematics": [C(0.05, 0.2, 0.0, 100.0), C(0.05, 0.2, 1.0, 100.0)],
    # reliability
    "reliability.age_replacement_cost_rate": [
        C(100.0, 1.0, 10.0, 1000.0, 2.5),
        C(500.0, 1.0, 10.0, 1000.0, 2.5),
        C(2000.0, 1.0, 10.0, 1000.0, 2.5),
    ],
    "reliability.arrhenius_acceleration_factor": [
        C(1.1215236438e-19, 330.0, 330.0),
        C(0.0, 300.0, 400.0),
        C(1.1215236438e-19, 300.0, 350.0),
    ],
    "reliability.availability": [C(100, 0), C(100, 25), C(1000.0, 20.0)],
    "reliability.bathtub_hazard": [C(0)],
    "reliability.bayesian_failure_rate": [C(2.0, 500000.0, 4, 2000000.0)],
    "reliability.expected_failures": [C(0.5, 2)],
    "reliability.exponential_failure_probability": [C(1e-12, 1)],
    "reliability.exponential_hazard": [C(0.5)],
    "reliability.exponential_mttf": [C(0.5)],
    "reliability.exponential_reliability": [C(2.0, 0.2)],
    "reliability.eyring_acceleration_factor": [C(8.01088317e-20, 300.0, 360.0)],
    "reliability.failure_mode_criticality": [C(0.5, 1.0, 1e-05, 100.0), C(0.3, 1.0, 1e-05, 100.0), C(0.2, 1.0, 1e-05, 100.0)],
    "reliability.failure_probability_from_index": [C(-1.0), C(0.0), C(3.0)],
    "reliability.failure_rate_from_mtbf": [C(10)],
    "reliability.gamma_hazard": [C(0.0, 0.5, 2.0), C(0.0, 1.0, 2.0), C(0.0, 3.0, 2.0)],
    "reliability.gamma_mttf": [C(0.5, 300.0), C(1.0, 300.0), C(2.5, 300.0)],
    "reliability.gamma_prior_from_bounds": [C(1e-06, 3e-06, 0.95), C(2e-05, 2.2e-05, 0.9), C(1e-06, 0.0001, 0.99)],
    "reliability.gamma_reliability": [C(1.0, 0.5, 300.0), C(233.06516869, 0.5, 300.0), C(0.004290645426, 0.5, 300.0)],
    "reliability.inverse_power_law_acceleration_factor": [C(100.0, 100.0, 2.0), C(100.0, 200.0, 2.0), C(100.0, 400.0, 2.0)],
    "reliability.k_of_n_reliability": [C(5, 3, 0.8)],
    "reliability.lognormal_hazard": [C(0.0, 7.60090245954, 0.6), C(485165195.41, 0.0, 0.5)],
    "reliability.lognormal_mttf": [C(7.60090245954, 0.6)],
    "reliability.lognormal_reliability": [C(0.0, 7.60090245954, 0.6), C(1.0, 7.60090245954, 0.6), C(233.06516869, 7.60090245954, 0.6)],
    "reliability.maintainability": [C(0, 0.5), C(2, 0.5)],
    "reliability.markov_availability": [C(1000000000.0, 0.001, 0.05)],
    "reliability.mtbf_confidence_bounds": [C(4947.33909795, 4, 0.9), C(4098.87281837, 4, 0.9), C(5954.18616487, 4, 0.9)],
    "reliability.parallel_reliability": [C(0.7, 0.8, 0.9), C(0.9, 0.8, 0.7)],
    "reliability.poisson_failures_probability": [C(0, 0, 1), C(0, 2, 1)],
    "reliability.reliability_growth_duane": [C(4, 2, 0.5)],
    "reliability.reliability_index_from_probability": [C(0.005)],
    "reliability.renewal_interval_cost_rate": [C(10, 20, 0.1, 2), C(0.0, 10.0, 0.001, 1.0)],
    "reliability.repair_rate": [C(2)],
    "reliability.safety_integrity_level": [C(0.00876), C(0.001)],
    "reliability.safety_integrity_pfd": [C(0.001, 10), C(2e-06, 8760.0)],
    "reliability.series_reliability": [C(0.7, 0.8, 0.9), C(0.7, 0.80001, 0.9), C(0.7, 0.79999, 0.9)],
    "reliability.standby_mttf": [C(0.0002, 1)],
    "reliability.standby_reliability": [C(5000.0, 0.0002, 4), C(1165325.84345, 0.0002, 4), C(21.45322713, 0.0002, 4)],
    "reliability.stress_strength_reliability": [C(500, 30, 400, 40)],
    "reliability.unavailability": [C(100, 25)],
    "reliability.weibull_cdf": [C(1e-12, 1, 1), C(22.9487231975, 73.5260741878, 1.93267798498), C(354.663665916, 1000.0, 2.5)],
    "reliability.weibull_hazard": [C(0, 2, 0.5), C(0, 2, 1), C(0, 2, 2)],
    "reliability.weibull_mttf": [C(3.0, 0.5), C(3.0, 1), C(3.0, 2)],
    "reliability.weibull_reliability": [C(1.0, 3.0, 0.5), C(233.06516869, 3.0, 0.5), C(0.004290645426, 3.0, 0.5)],
    "reliability.gamma_cdf": [C(1.0, 1.0, 1.0), C(2.0, 3.0, 0.5)],
    "reliability.lognormal_cdf": [C(2.718281828459045, 0.0, 1.0)],
    # signal_processing
    "signal_processing.alias_frequency": [C(70, 100)],
    "signal_processing.db20": [C(10)],
    "signal_processing.equiripple_order": [C(0.01, 0.001, 50.0, 1000.0)],
    "signal_processing.iir_order": [
        C(1000.0, 1500.0, 1.0, 60.0, 8000.0, "butterworth"),
        C(300.0, 350.0, 0.1, 40.0, 8000.0, "butterworth"),
        C(2000.0, 3500.0, 3.0, 20.0, 8000.0, "butterworth"),
    ],
    "signal_processing.kaiser_order": [C(60.0, 50.0, 1000.0), C(40.0, 50.0, 1000.0), C(15.0, 50.0, 1000.0)],
    "signal_processing.nyquist_frequency": [C(100)],
    "signal_processing.prewarp_frequency": [C(100.0, 1000.0)],
    "signal_processing.quantization_step": [C(-1, 1, 8)],
    "signal_processing.quantize": [C(-100, -1, 1, 3), C(100, -1, 1, 3)],
    "signal_processing.sampling_period": [C(10)],
    # structural_engineering
    "structural_engineering.angle_section": [C(0.15, 0.1, 0.01)],
    "structural_engineering.axial_deformation": [C(1000.0, 2.0, 0.01, 200000000000.0)],
    "structural_engineering.bending_stress": [C(10, 2, 4)],
    "structural_engineering.box_section": [C(0.2, 0.1, 0.01)],
    "structural_engineering.catenary_cable": [C(2000.0, 100.0, 12.0), C(2000.0, 100.0, 0.5)],
    "structural_engineering.channel_section": [C(0.2, 0.08, 0.01, 0.008)],
    "structural_engineering.channel_shear_center": [C(0.2, 0.08, 0.01, 0.006), C(0.2, 0.08, 0.01, 0.01), C(0.2, 0.0, 0.01, 0.006)],
    "structural_engineering.euler_buckling_load": [C(1, 1, 1), C(200000000000.0, 6.4e-06, 3.0), C(200000000000.0, 6.4e-06, 30.0)],
    "structural_engineering.hooke_stress": [C(200000000000.0, 0.001)],
    "structural_engineering.i_section": [C(0.3, 0.15, 0.012, 0.008)],
    "structural_engineering.influence_line": [
        C(10.0, 1.0, "moment", at=3.0),
        C(10.0, 3.0, "moment", at=3.0),
        C(10.0, 7.5, "moment", at=3.0),
    ],
    "structural_engineering.johnson_critical_stress": [C(250000000.0, 200000000000.0, 0.0), C(250000000.0, 200000000000.0, 200.0)],
    "structural_engineering.maxwell_betti_reciprocity": [C(0.001, 0.001)],
    "structural_engineering.mohr_circle": [C(100, -20, 35), C(80.0, 40.0, 0.0), C(50.0, -20.0, 25.0)],
    "structural_engineering.neutral_axis_angle": [C(2000.0, -5000.0, 2.02625e-06, 5.57625e-06, -1.96875e-06)],
    "structural_engineering.parabolic_cable": [C(2000.0, 100.0, 12.0), C(2000.0, 100.0, 0.5)],
    "structural_engineering.parallel_axis": [C(1, 2, 3)],
    "structural_engineering.plane_stress_from_strains": [C(0.00042, -0.00013, 0.00026, 70000000000.0, 0.33)],
    "structural_engineering.plastic_collapse_load": [C(100000.0, 6.0, "simple", "uniform")],
    "structural_engineering.plastic_section_modulus": [C(0.1, 0.3)],
    "structural_engineering.poisson_lateral_strain": [C(0.01, 0.3)],
    "structural_engineering.principal_stresses_3d": [C(50.0, -20.0, 30.0, 25.0, -10.0, 15.0)],
    "structural_engineering.radius_of_gyration": [C(4, 1)],
    "structural_engineering.rectangular_second_moment": [C(12, 1), C(0.05, 0.25), C(0.01, 0.25)],
    "structural_engineering.rectangular_section_modulus": [C(0.1, 0.3), C(0.2, 0.3)],
    "structural_engineering.secant_formula_load": [
        C(250000000.0, 0.004, 0.01, 0.05, 0.04, 3.0, 200000000000.0),
        C(250000000.0, 0.004, 0.0, 0.05, 0.04, 3.0, 200000000000.0),
        C(250000000.0, 0.004, 0.0, 0.05, 0.04, 30.0, 200000000000.0),
    ],
    "structural_engineering.secant_formula_stress": [
        C(100000.0, 0.004, 0.0, 0.05, 0.04, 3.0, 200000000000.0),
        C(1000.0, 0.004, 0.01, 0.05, 0.04, 0.001, 200000000000.0),
        C(1403675.6667, 0.004, 0.001, 0.05, 0.04, 3.0, 200000000000.0),
    ],
    "structural_engineering.shear_flow": [C(30000.0, 0.0002592, 8.8709184e-05)],
    "structural_engineering.slenderness_ratio": [C(10, 1, 2)],
    "structural_engineering.strain_energy_axial": [
        C(1001.0, 2.0, 0.01, 200000000000.0),
        C(999.0, 2.0, 0.01, 200000000000.0),
        C(1000, 2, 0.01, 200000000000.0),
    ],
    "structural_engineering.strain_energy_bending_constant_moment": [C(10, 2, 5, 4)],
    "structural_engineering.stress_transformation": [
        C(80000000.0, -20000000.0, 30000000.0, 0.270209750135),
        C(80000000.0, -20000000.0, 30000000.0, 1.05560791353),
    ],
    "structural_engineering.t_section": [C(0.2, 0.15, 0.02, 0.01)],
    "structural_engineering.tangent_modulus_load": [C(200000000000.0, 8e-06, 2.0)],
    "structural_engineering.thermal_strain": [C(1e-05, 100)],
    "structural_engineering.thin_wall_vessel_stresses": [C(1000000.0, 1.0, 0.005)],
    "structural_engineering.thin_walled_torsion_angle": [C(2000.0, 1.5, 80000000000.0, 0.00785398163397, 157.079632679)],
    "structural_engineering.thin_walled_torsion_shear_stress": [
        C(2000.0, 0.00785398163397, 0.002),
        C(2000.0, 0.0171, 0.005),
        C(2000.0, 0.0171, 0.01),
    ],
    "structural_engineering.torsion_angle": [C(10, 2, 5, 4), C(2000.0, 1.5, 80000000000.0, 1.57079632679e-06)],
    "structural_engineering.torsion_shear_stress": [C(10, 2, 4), C(2000.0, 0.05, 1.57079632679e-06)],
    "structural_engineering.transverse_shear_stress": [
        C(30000.0, 0.001125, 0.000225, 0.1),
        C(30000.0, 0.0001323, 8.8709184e-05, 0.15),
        C(30000.0, 3.5207803032e-06, 8.8709184e-05, 0.15),
    ],
    "structural_engineering.unsymmetrical_bending_stress": [
        C(2000.0, -5000.0, -0.830746878172, 0.556650361006, 2.02625e-06, 5.57625e-06, -1.96875e-06)
    ],
    "structural_engineering.thick_wall_cylinder_stresses": [C(1.5, 1.0, 2.0, 1e6), C(1.0, 1.0, 2.0, 1e6, 2e5)],
    # thermodynamics
    "thermodynamics.brayton_efficiency": [C(6, 1.3), C(6, 1.4), C(6, 1.66666666667)],
    "thermodynamics.carnot_cop_heat_pump": [C(600, 300)],
    "thermodynamics.carnot_cop_refrigerator": [C(600, 300)],
    "thermodynamics.carnot_efficiency": [C(600, 300), C(800.0, 300.0)],
    "thermodynamics.clausius_clapeyron_pressure": [C(368.15, 373.15, 101418.0, 2257000.0, 461.5)],
    "thermodynamics.curzon_ahlborn_efficiency": [C(800.0, 300.0)],
    "thermodynamics.dew_point": [C(0.00525740566156), C(0.0133123109583), C(0.0243796231657)],
    "thermodynamics.diesel_efficiency": [C(8, 1), C(8, 1.000000001), C(8, 2)],
    "thermodynamics.entropy_change_isothermal_ideal_gas": [C(2.0, 350.0, 0.01, 0.04)],
    "thermodynamics.exergy_destruction": [C(300, 2)],
    "thermodynamics.exergy_heat": [C(100, 600, 300)],
    "thermodynamics.first_law_heat": [C(10, 20)],
    "thermodynamics.heat_engine_work": [C(1000, 500.0)],
    "thermodynamics.heat_pump_cop": [C(150, 50)],
    "thermodynamics.humidity_ratio": [C(293.15, 1.0), C(293.15, 0.0)],
    "thermodynamics.ideal_gas_pressure": [C(2.0, 350.0, 0.025), C(2.0, 350.0, 0.0103914020722), C(2.0, 350.0, 0.0396085979278)],
    "thermodynamics.ideal_gas_temperature": [C(582012.383271, 0.01, 2.0)],
    "thermodynamics.isentropic_pressure_volume": [C(100, 1, 2, 1.4), C(100, 1, 2, 1.3), C(100, 1, 2, 1.66666666667)],
    "thermodynamics.isentropic_temperature_ratio": [C(3, 1.3), C(3, 1.4), C(3, 1.66666666667)],
    "thermodynamics.joule_thomson_coefficient_vdw": [
        C(300.0, 0.137, 3.87e-05, 29.1),
        C(851.540704954, 0.137, 3.87e-05, 29.1),
        C(1703.08140991, 0.137, 3.87e-05, 29.1),
    ],
    "thermodynamics.otto_efficiency": [C(8.0, 1.4), C(8, 1.3), C(8, 1.66666666667)],
    "thermodynamics.rankine_ideal_efficiency": [
        C(210000.0, 3200000.0, 2100000.0, 200000.0),
        C(201904.798315, 3375058.44185, 2089640.34743, 191812.295194),
        C(200, 3200, 2400, 190),
    ],
    "thermodynamics.refrigeration_cop": [C(100, 50)],
    "thermodynamics.relative_humidity": [C(303.15, 0.00525740566156), C(303.15, 0.0133123109583), C(303.15, 0.0243796231657)],
    "thermodynamics.reversible_isothermal_work_ideal_gas": [C(2.0, 350.0, 0.01, 0.04), C(1, 300, 1, 2)],
    "thermodynamics.saturation_pressure": [C(450.0), C(450.01), C(449.99)],
    "thermodynamics.saturation_temperature": [C(100000.0), C(1000000.0), C(10000000.0)],
    "thermodynamics.specific_enthalpy_ideal_gas": [C(300, 1000)],
    "thermodynamics.specific_internal_energy_ideal_gas": [C(300, 700)],
    "thermodynamics.steam_properties": [C(5000100.0, s=6743.13045961), C(4999900.0, s=6743.13045961), C(10000.0, s=8155.22062555)],
    "thermodynamics.thermal_efficiency": [C(100, 50)],
    "thermodynamics.van_der_waals_critical_point": [C(0.364, 4.267e-05)],
    "thermodynamics.van_der_waals_pressure": [
        C(1, 300, 0.024, 0, 0),
        C(1, 300, 10.0, 0.364, 4.267e-05),
        C(1.0, 303.998036511, 0.00012801, 0.364, 4.267e-05),
    ],
    "thermodynamics.antoine_vapor_pressure": [C(300.0, 10.196212, 1730.63, -39.724), C(350.0, 10.196212, 1730.63, -39.724)],
    "thermodynamics.adiabatic_flame_temperature": [
        C(1, 4, 0.0, h_f_fuel=-74870.0, excess_air=0.0),
        C(8, 18, 0.0, h_f_fuel=-249950.0, excess_air=3.0),
    ],
    "thermodynamics.adiabatic_mixing": [C(1.0, 305.15, 0.015, 3.0, 287.15, 0.006), C(2.0, 300.15, 0.01, 1.0, 290.15, 0.005)],
    "thermodynamics.adiabatic_saturation_temperature": [C(303.15, 0.01), C(298.15, 0.00653)],
    "thermodynamics.brayton_cycle": [
        C(8.0, 300.0, 1300.0, regenerator_effectiveness=0.0),
        C(4.0, 300.0, 1300.0, regenerator_effectiveness=1.0),
    ],
    "thermodynamics.combustion_products": [C(8, 18, 0.0, excess_air=0.5), C(1, 4, 0.0, excess_air=0.0)],
    "thermodynamics.compressibility_factor": [C(101325.0, 273.15, 0.022414), C(2.5e6, 510.0, 0.0015)],
    "thermodynamics.compression_outlet_enthalpy": [C(300000.0, 540000.0, 0.8), C(300000.0, 540000.0, 1.0)],
    "thermodynamics.compressor_isentropic_efficiency": [C(300000.0, 600000.0, 540000.0), C(300000.0, 540000.0, 540000.0)],
    "thermodynamics.cooling_coil": [C(303.15, 0.3, 298.15), C(303.15, 0.8, 287.15)],
    "thermodynamics.cooling_tower": [C(100.0, 308.15, 295.15, 293.15, 0.6, 303.15), C(50.0, 310.15, 298.15, 295.15, 0.5, 305.15)],
    "thermodynamics.cp_minus_cv": [C(350.0, 0.50225, 0.00285714285714, 5e-06), C(277.13, 0.001, 0.0, 4.6e-10)],
    "thermodynamics.dual_efficiency": [C(16.0, 1.0, 2.0), C(16.0, 1.8, 1.0), C(16.0, 1.5, 1.4, 1.4)],
    "thermodynamics.expansion_outlet_enthalpy": [C(3000000.0, 2500000.0, 0.6), C(3000000.0, 2500000.0, 0.85)],
    "thermodynamics.fuel_molar_mass": [C(1, 4), C(8, 18)],
    "thermodynamics.heating_value": [C(1, 4, 0.0, h_f_fuel=-74870.0), C(8, 18, 0.0, h_f_fuel=-249950.0)],
    "thermodynamics.joule_thomson_coefficient": [C(350.0, 0.50225, 0.00285714285714, 1005.0), C(300.0, 0.0249, 0.00666666666667, 1040.0)],
    "thermodynamics.moist_air_enthalpy": [C(303.15, 0.012), C(273.15, 0.0)],
    "thermodynamics.nozzle_exit_velocity": [C(3000000.0, 2900000.0, 50.0), C(3000000.0, 2950000.0, 0.0)],
    "thermodynamics.nozzle_isentropic_efficiency": [C(3000000.0, 2955000.0, 2950000.0, 0.0)],
    "thermodynamics.peng_robinson_compressibility": [C(374.2, 10.0, 374.2, 4.059e6, 0.327), C(510.0, 2.5e6, 425.1, 3.796e6, 0.2)],
    "thermodynamics.peng_robinson_pressure": [C(300.0, 0.024, 190.56, 4.599e6, 0.011), C(374.2, 0.000236, 374.2, 4.059e6, 0.327)],
    "thermodynamics.peng_robinson_saturation_pressure": [C(280.0, 374.2, 4.059e6, 0.327), C(233.15, 374.2, 4.059e6, 0.327)],
    "thermodynamics.pitzer_acentric_factor": [C(486000.0, 4.86e6), C(1.0e6, 22.064e6)],
    "thermodynamics.pitzer_virial_compressibility": [C(510.0, 2.5e6, 425.1, 3.796e6, 0.2), C(600.0, 1e6, 425.1, 3.796e6, 0.2)],
    "thermodynamics.rankine_cycle": [C(15e6, 10e3, 873.15), C(10e6, 10e3, 773.15)],
    "thermodynamics.steady_flow_work": [C(3000000.0, 2500000.0, q=0.0, V_in=0.0, z_in=0.0), C(0.0, 0.0, q=-10.0, V_in=10.0, z_in=1.0)],
    "thermodynamics.stirling_efficiency": [C(900.0, 300.0, 4.0, 1.4, 0.0), C(900.0, 300.0, 4.0, 1.4, 0.5), C(900.0, 300.0, 2.0, 1.4, 1.0)],
    "thermodynamics.stoichiometric_air_fuel_ratio": [C(1, 4, 0), C(8, 18, 0), C(1, 4, 1)],
    "thermodynamics.turbine_isentropic_efficiency": [C(3000000.0, 2600000.0, 2500000.0)],
}

# Scalar parameters, but the result is an array, a system or an airfoil, not a number.
NOT_SCALAR_RESULT = {
    "aerodynamics.naca4_airfoil",
    "control_systems.pade_delay",
    "mechanical_engineering.dh_transform",
    "signal_processing.fft_frequencies",
    "signal_processing.butterworth_lowpass",
    "signal_processing.butterworth_highpass",
    "signal_processing.fir_lowpass",
    "communications.psk_constellation",
    "control_systems.lag_compensator",
    "control_systems.lead_compensator",
    "control_systems.lead_lag_compensator",
    "heat_transfer.conduction_2d",
    "materials.diffraction_peaks",
    "materials.lamina_stiffness",
    "materials.principal_stresses",
    "signal_processing.chebyshev1_lowpass",
    "signal_processing.chebyshev2_lowpass",
    "signal_processing.elliptic_lowpass",
    "signal_processing.equiripple_lowpass",
}

# Added by ROADMAP D items still in progress on 2026-10-23. Each one fails the
# array test; the item that added it vectorizes it or moves it to SCALAR_ONLY.
PENDING = {
    "aerodynamics.lifting_line",
    "aerodynamics.naca4_moment_coefficient",
    "aerodynamics.oblique_shock",
    "electrical_engineering.diode_operating_point",
    "electrical_engineering.rlc_damping",
    "electrical_engineering.rlc_step_voltage",
    "electrical_engineering.swing_response",
    "fluid_engineering.alternate_depths",
    "fluid_engineering.jet_force",
    "fluid_engineering.pump_operating_point",
    "fluid_engineering.terminal_velocity",
    "heat_transfer.cylinder_transient",
    "heat_transfer.exchanger_effectiveness_crossflow",
    "heat_transfer.lmtd_correction_factor",
    "heat_transfer.rohsenow_boiling_heat_flux",
    "heat_transfer.sphere_transient",
    "heat_transfer.tube_bank_max_velocity",
    "heat_transfer.zukauskas_nusselt",
    "mechanical_engineering.bolt_stiffness",
    "mechanical_engineering.bolted_joint",
    "mechanical_engineering.member_stiffness",
    "mechanical_engineering.shaft_diameter_goodman",
    "mechanical_engineering.slider_crank_kinematics",
    "signal_processing.equiripple_order",
    "signal_processing.iir_order",
    "structural_engineering.angle_section",
    "structural_engineering.box_section",
    "structural_engineering.catenary_cable",
    "structural_engineering.channel_section",
    "structural_engineering.i_section",
    "structural_engineering.secant_formula_load",
    "structural_engineering.t_section",
    "structural_engineering.thick_wall_cylinder_stresses",
}

# Parameters that stay scalar-only, and why; each docstring says so. The
# function's other float parameters still broadcast.
SCALAR_ONLY = {
    # one quadrature of the camber line per airfoil
    "aerodynamics.naca4_zero_lift_angle": ("m", "p"),
    # the eigenvalues are roots found for one Biot number
    "heat_transfer.plane_wall_transient": ("Bi",),
    # the front coefficient is a root found for one Stefan number
    "heat_transfer.neumann_stefan_front": ("Ste",),
    # integrates one crack-growth history
    "materials.crack_growth_life": ("C", "m", "delta_sigma", "a0", "af", "Y"),
    # one posterior distribution per update
    "reliability.bayesian_failure_rate": ("prior_shape", "prior_rate", "total_time"),
    # the prior shape is a root found for one pair of bounds
    "reliability.gamma_prior_from_bounds": ("median", "upper", "confidence"),
    # one rotor, swept over omega
    "mechanical_engineering.jeffcott_response": ("m", "k", "e", "c"),
    # designs one absorber
    "mechanical_engineering.tuned_mass_damper": ("m", "k", "mass_ratio", "frequency_ratio", "damping_ratio"),
    # sizes one FIR filter
    "signal_processing.kaiser_order": ("attenuation_db", "transition_width", "fs"),
    # one cycle solved state by state, with Brent inversions of IF97 at each state
    "thermodynamics.rankine_cycle": (
        "p_boiler",
        "p_condenser",
        "T_turbine_inlet",
        "reheat_pressure",
        "T_reheat",
        "feedwater_pressure",
        "turbine_efficiency",
        "pump_efficiency",
    ),
    # one gas-turbine cycle with its list of states
    "thermodynamics.brayton_cycle": (
        "pressure_ratio",
        "T_min",
        "T_max",
        "compressor_efficiency",
        "turbine_efficiency",
        "regenerator_effectiveness",
        "cp",
        "gamma",
        "p_min",
    ),
    # a cubic root chosen for one state
    "thermodynamics.peng_robinson_compressibility": ("T", "p", "Tc", "pc", "omega"),
    # Brent's method between the spinodals of one isotherm
    "thermodynamics.peng_robinson_saturation_pressure": ("T", "Tc", "pc", "omega"),
    # Brent's method for one air state
    "thermodynamics.adiabatic_saturation_temperature": ("T", "w", "p"),
    # Brent's method on one product enthalpy balance
    "thermodynamics.adiabatic_flame_temperature": ("c", "h", "o", "h_f_fuel", "excess_air", "T_air"),
}

_SCALAR_TYPE = re.compile(r"(float|complex|int|str|bool|\{)")
_ARRAY_TYPE = re.compile(r"(float|complex)")


def _function(name):
    domain, attr = name.split(".")
    return getattr(importlib.import_module(f"engineeringkit.{domain}"), attr)


def _parameter_types(f):
    doc = inspect.getdoc(f) or ""
    section = re.search(r"Parameters\n-+\n(.*?)(?:\n\n\w[\w ]*\n-+\n|\Z)", doc, re.S)
    types = {}
    for names, kind in re.findall(r"^(\S[^:\n]*?) : (.+)$", section.group(1) if section else "", re.M):
        for p in names.split(","):
            types[p.strip().lstrip("*")] = kind
    return types


def _scalar_functions():
    names = set()
    for domain in DOMAINS:
        module = importlib.import_module(f"engineeringkit.{domain}")
        for attr in module.__all__:
            f = getattr(module, attr)
            if inspect.isfunction(f) and not attr.startswith("plot_"):
                kinds = [_parameter_types(f).get(p, "array_like") for p in inspect.signature(f).parameters]
                if all(_SCALAR_TYPE.match(k) for k in kinds) and any(_ARRAY_TYPE.match(k) for k in kinds):
                    names.add(f"{domain}.{attr}")
    return names - NOT_SCALAR_RESULT


def _leaves(result, path=()):
    """Flatten a result into ``{path: number or array}``."""
    if isinstance(result, tuple):
        return {k: v for i, item in enumerate(result) for k, v in _leaves(item, (*path, i)).items()}
    if isinstance(result, dict):
        return {k: v for key, item in result.items() for k, v in _leaves(item, (*path, key)).items()}
    if dataclasses.is_dataclass(result):
        return {k: v for f in dataclasses.fields(result) for k, v in _leaves(getattr(result, f.name), (*path, f.name)).items()}
    if isinstance(result, (numbers.Number, np.ndarray, np.bool_)):
        return {path: result}
    return {}


def _stack(values):
    return np.array(values, dtype=np.result_type(*values, float))


def _check_broadcasts(name, scalar=()):
    """Call ``name`` once on its cases stacked into 1-D arrays and compare with the scalar calls.

    Parameters in ``scalar`` keep the first case's value, and result fields
    that depend only on them may stay scalars.
    """
    f = _function(name)
    signature, types = inspect.signature(f), _parameter_types(f)
    calls = [signature.bind(*args, **kwargs) for args, kwargs in (CASES[name] * 3)[:3]]
    stacked = signature.bind(*CASES[name][0][0], **CASES[name][0][1])
    for p, first in stacked.arguments.items():
        assert all(p in call.arguments for call in calls), f"{name}: every case must pass {p}"
        if p in scalar or not _ARRAY_TYPE.match(types[p]):
            for call in calls:
                assert p in scalar or call.arguments[p] == first, f"{name}: {p} is documented {types[p]!r} and must not vary"
                call.arguments[p] = first
        elif signature.parameters[p].kind is inspect.Parameter.VAR_POSITIONAL:
            assert len({len(call.arguments[p]) for call in calls}) == 1, f"{name}: every case must pass as many *{p}"
            stacked.arguments[p] = tuple(_stack(column) for column in zip(*(call.arguments[p] for call in calls)))
        else:
            stacked.arguments[p] = _stack([call.arguments[p] for call in calls])
    got = _leaves(f(*stacked.args, **stacked.kwargs))
    expected = [_leaves(f(*call.args, **call.kwargs)) for call in calls]
    assert expected[0], f"{name} returned no numbers"
    assert set(got) == set(expected[0]), f"{name}: array result has fields {sorted(got)}, scalar result {sorted(expected[0])}"
    magnitudes = np.abs([v for e in expected for v in np.ravel(list(e.values())) if np.isfinite(v)])
    atol = 1e-12 * magnitudes.max(initial=0)
    for path, value in got.items():
        want = np.array([e[path] for e in expected])
        if scalar and np.ndim(value) == 0:
            value = np.broadcast_to(value, want.shape)
        assert np.shape(value) == want.shape, f"{name}{list(path)}: shape {np.shape(value)}, expected {want.shape}"
        np.testing.assert_allclose(value, want, rtol=1e-9, atol=atol, err_msg=name)


def _tested(name):
    """Whether any domain test calls ``name`` (functions still being written get their case with their tests)."""
    domain, attr = name.split(".")
    sources = [path.read_text() for path in (Path(__file__).parents[1] / domain / "tests").glob("test_*.py")]
    return any(re.search(rf"\b{attr}\(", text) for text in sources)


def test_every_scalar_function_has_cases():
    scalar = _scalar_functions()
    assert sorted(name for name in scalar - set(CASES) if _tested(name)) == []
    assert sorted(set(CASES) - scalar) == []
    assert set(SCALAR_ONLY) | PENDING <= set(CASES)
    assert not set(SCALAR_ONLY) & PENDING


@pytest.mark.parametrize(
    "name",
    [
        pytest.param(name, marks=pytest.mark.xfail(name in PENDING, reason="pending vectorization", strict=True))
        for name in sorted(set(CASES) - set(SCALAR_ONLY))
    ],
)
def test_scalar_function_broadcasts_over_1d_arrays(name):
    _check_broadcasts(name)


@pytest.mark.parametrize("name", sorted(SCALAR_ONLY))
def test_scalar_only_parameters_are_documented_and_needed(name):
    f = _function(name)
    assert "scalar-only" in (inspect.getdoc(f) or ""), f"{name} must say in its docstring which parameters are scalar-only"
    with pytest.raises(Exception):
        _check_broadcasts(name)
    types = _parameter_types(f)
    if any(_ARRAY_TYPE.match(types[p]) and p not in SCALAR_ONLY[name] for p in inspect.signature(f).parameters):
        _check_broadcasts(name, scalar=SCALAR_ONLY[name])

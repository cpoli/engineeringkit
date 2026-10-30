"""engineeringkit.electrical_engineering: circuits from Ohm's law to the digital load flow.

Ohm's law, series and parallel equivalents, the Wheatstone bridge, modified
nodal analysis and mesh analysis of resistive networks; Thévenin, Norton
and Jacobi's maximum power transfer; RC and RL transients and the over-,
critically and underdamped RLC step; complex impedance, series resonance
and the ideal transformer; ideal op-amp circuits, Shockley's diode equation
and the load line; DC machine and induction motor torque-speed curves;
lossless and lossy transmission lines, Heaviside's distortionless
condition, reflections and input impedance for the Smith chart; and real,
reactive and three-phase power, the per-unit system, Fortescue's
symmetrical components and fault currents, Gauss-Seidel and Newton-Raphson
load flow, and the swing equation with the equal-area criterion. Models are
hand-written; NumPy solves the linear systems and SciPy computes the RLC
step response, integrals, roots and swing curves.
"""

from engineeringkit import __version__
from engineeringkit.electrical_engineering.core.base import CircuitResult, LoadFlowResult, MeshResult, NodalResult, SwingResult
from engineeringkit.electrical_engineering.systems.ac_circuits import (
    SeriesRLC,
    capacitor_impedance,
    inductor_impedance,
    resonance_frequency,
    rms_sine,
    series_rlc_impedance,
    transformer_secondary,
)
from engineeringkit.electrical_engineering.systems.dc_circuits import (
    capacitors_parallel,
    capacitors_series,
    mesh_analysis,
    nodal_analysis,
    ohms_law,
    parallel_resistance,
    series_resistance,
    wheatstone_balanced,
)
from engineeringkit.electrical_engineering.systems.electronics import (
    diode_operating_point,
    inverting_amplifier_gain,
    noninverting_amplifier_gain,
    op_amp_integrator,
    shockley_diode_current,
    thermal_voltage,
)
from engineeringkit.electrical_engineering.systems.machines import (
    dc_motor_speed,
    dc_motor_torque,
    induction_motor_breakdown,
    induction_motor_torque,
)
from engineeringkit.electrical_engineering.systems.networks import maximum_power_transfer, norton_from_thevenin, thevenin_from_source
from engineeringkit.electrical_engineering.systems.power import (
    admittance_matrix,
    change_base,
    fault_currents,
    gauss_seidel_load_flow,
    newton_raphson_load_flow,
    per_unit_base,
    phase_components,
    power_factor,
    reactive_power,
    real_power,
    symmetrical_components,
    three_phase_power,
)
from engineeringkit.electrical_engineering.systems.stability import critical_clearing_time, equal_area_clearing_angle, swing_response
from engineeringkit.electrical_engineering.systems.transients import (
    rc_charging_voltage,
    rc_time_constant,
    rl_current,
    rl_time_constant,
    rlc_damping,
    rlc_step_voltage,
)
from engineeringkit.electrical_engineering.systems.transmission_lines import (
    distortionless_condition,
    input_impedance,
    propagation_constant,
    reflection_coefficient,
    standing_wave_ratio,
    transmission_line_velocity,
    transmission_line_z0,
)
from engineeringkit.electrical_engineering.utils.phasors import phasor, polar

__all__ = [
    "__version__",
    "CircuitResult",
    "NodalResult",
    "MeshResult",
    "LoadFlowResult",
    "SwingResult",
    "phasor",
    "polar",
    "ohms_law",
    "series_resistance",
    "parallel_resistance",
    "capacitors_parallel",
    "capacitors_series",
    "wheatstone_balanced",
    "nodal_analysis",
    "mesh_analysis",
    "thevenin_from_source",
    "norton_from_thevenin",
    "maximum_power_transfer",
    "rc_time_constant",
    "rl_time_constant",
    "rc_charging_voltage",
    "rl_current",
    "rlc_damping",
    "rlc_step_voltage",
    "capacitor_impedance",
    "inductor_impedance",
    "series_rlc_impedance",
    "resonance_frequency",
    "rms_sine",
    "transformer_secondary",
    "SeriesRLC",
    "inverting_amplifier_gain",
    "noninverting_amplifier_gain",
    "op_amp_integrator",
    "thermal_voltage",
    "shockley_diode_current",
    "diode_operating_point",
    "dc_motor_torque",
    "dc_motor_speed",
    "induction_motor_torque",
    "induction_motor_breakdown",
    "transmission_line_z0",
    "transmission_line_velocity",
    "propagation_constant",
    "distortionless_condition",
    "reflection_coefficient",
    "standing_wave_ratio",
    "input_impedance",
    "real_power",
    "reactive_power",
    "power_factor",
    "three_phase_power",
    "symmetrical_components",
    "phase_components",
    "per_unit_base",
    "change_base",
    "fault_currents",
    "admittance_matrix",
    "gauss_seidel_load_flow",
    "newton_raphson_load_flow",
    "equal_area_clearing_angle",
    "critical_clearing_time",
    "swing_response",
]

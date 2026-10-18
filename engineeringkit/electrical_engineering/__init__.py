"""engineeringkit.electrical_engineering: circuits from Ohm's law to the digital load flow.

Ohm's law, series and parallel equivalents, the Wheatstone bridge and
modified nodal analysis of resistive networks; Thévenin, Norton and Jacobi's
maximum power transfer; RC and RL transients; complex impedance, series
resonance and the ideal transformer; lossless and lossy transmission lines,
Heaviside's distortionless condition, reflections and input impedance for
the Smith chart; and real, reactive and three-phase power, Fortescue's
symmetrical components and Gauss-Seidel load flow. Models are hand-written;
NumPy solves the linear systems and SciPy computes the RLC step response.
"""

from engineeringkit import __version__
from engineeringkit.electrical_engineering.core.base import CircuitResult, LoadFlowResult, NodalResult
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
    nodal_analysis,
    ohms_law,
    parallel_resistance,
    series_resistance,
    wheatstone_balanced,
)
from engineeringkit.electrical_engineering.systems.networks import maximum_power_transfer, norton_from_thevenin, thevenin_from_source
from engineeringkit.electrical_engineering.systems.power import (
    admittance_matrix,
    gauss_seidel_load_flow,
    phase_components,
    power_factor,
    reactive_power,
    real_power,
    symmetrical_components,
    three_phase_power,
)
from engineeringkit.electrical_engineering.systems.transients import rc_charging_voltage, rc_time_constant, rl_time_constant
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
    "LoadFlowResult",
    "phasor",
    "polar",
    "ohms_law",
    "series_resistance",
    "parallel_resistance",
    "capacitors_parallel",
    "capacitors_series",
    "wheatstone_balanced",
    "nodal_analysis",
    "thevenin_from_source",
    "norton_from_thevenin",
    "maximum_power_transfer",
    "rc_time_constant",
    "rl_time_constant",
    "rc_charging_voltage",
    "capacitor_impedance",
    "inductor_impedance",
    "series_rlc_impedance",
    "resonance_frequency",
    "rms_sine",
    "transformer_secondary",
    "SeriesRLC",
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
    "admittance_matrix",
    "gauss_seidel_load_flow",
]

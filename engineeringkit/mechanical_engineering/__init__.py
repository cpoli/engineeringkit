"""engineeringkit.mechanical_engineering: dynamics, machines, mechanisms and vibration.

Moments of inertia, work, energy and power, and Lagrange's equations; gear
trains, Lewis's gear-tooth bending stress and the Euler-Eytelwein
belt-friction law; Coulomb's power screw, Petroff's and Sommerfeld's journal
bearings, the Lundberg-Palmgren rolling-bearing life, shafts in combined
bending and torsion, helical springs with Wahl's factor and preloaded bolted
joints; the slider-crank, Grashof's condition, four-bar velocity and
acceleration analysis, Freudenstein's synthesis, cam displacement laws and
Denavit-Hartenberg forward kinematics; the damped and forced oscillator,
normal modes and modal superposition, Rayleigh's and Dunkerley's estimates,
the tuned vibration absorber, Holzer's torsional method, the Jeffcott rotor
and the balancing of rotating masses. Models
are hand-written closed forms; SciPy's ODE integrator and root finder are used
where they are the standard tool.
"""

from engineeringkit import __version__
from engineeringkit.mechanical_engineering.core.base import (
    BalanceResult,
    BoltedJoint,
    CamMotion,
    FourBarMotion,
    HolzerResult,
    JeffcottResult,
    MDOFResult,
    ModalResult,
    SDOFResult,
    SliderCrankMotion,
    TunedMassDamper,
)
from engineeringkit.mechanical_engineering.systems.balancing import balance_rotating_masses, unbalance
from engineeringkit.mechanical_engineering.systems.bolted_joints import bolt_stiffness, bolted_joint, member_stiffness, preload_torque
from engineeringkit.mechanical_engineering.systems.cams import cam_follower_motion
from engineeringkit.mechanical_engineering.systems.dynamics import (
    angular_momentum,
    centripetal_force,
    disk_inertia,
    flywheel_energy,
    kinetic_energy,
    parallel_axis,
    power_torque,
    rotational_energy,
    solid_sphere_inertia,
    work_constant_force,
)
from engineeringkit.mechanical_engineering.systems.lagrangian import lagrange_accelerations, simulate_lagrangian
from engineeringkit.mechanical_engineering.systems.machine_elements import (
    bearing_l10_life,
    long_bearing_eccentricity,
    petroff_friction_torque,
    screw_torque_raise,
    sommerfeld_number,
)
from engineeringkit.mechanical_engineering.systems.mdof import (
    dunkerley_frequency,
    harmonic_response,
    modal_response,
    normal_modes,
    rayleigh_frequency,
)
from engineeringkit.mechanical_engineering.systems.mechanisms import (
    dh_transform,
    forward_kinematics,
    fourbar_grashof,
    fourbar_kinematics,
    fourbar_output_angle,
    freudenstein_coefficients,
    freudenstein_link_lengths,
    slider_crank_displacement,
    slider_crank_kinematics,
)
from engineeringkit.mechanical_engineering.systems.rotordynamics import jeffcott_response, rotating_unbalance_force
from engineeringkit.mechanical_engineering.systems.shafts import (
    shaft_diameter,
    shaft_diameter_goodman,
    shaft_max_shear_stress,
    shaft_von_mises_stress,
)
from engineeringkit.mechanical_engineering.systems.springs import (
    bergstrasser_factor,
    spring_index,
    spring_rate,
    spring_shear_stress,
    wahl_factor,
)
from engineeringkit.mechanical_engineering.systems.transmission import (
    barth_velocity_factor,
    belt_speed,
    belt_tension_ratio,
    gear_output_speed,
    gear_output_torque,
    gear_ratio,
    lewis_bending_stress,
    lewis_form_factor,
)
from engineeringkit.mechanical_engineering.systems.vibration import (
    damped_frequency,
    damping_ratio,
    holzer_frequencies,
    holzer_residual,
    magnification_factor,
    natural_frequency,
    natural_frequency_hz,
    phase_angle,
    simulate_sdof,
    transmissibility,
    tuned_mass_damper,
)

__all__ = [
    "__version__",
    "SDOFResult",
    "TunedMassDamper",
    "JeffcottResult",
    "HolzerResult",
    "ModalResult",
    "MDOFResult",
    "BoltedJoint",
    "BalanceResult",
    "FourBarMotion",
    "SliderCrankMotion",
    "CamMotion",
    "kinetic_energy",
    "rotational_energy",
    "work_constant_force",
    "power_torque",
    "angular_momentum",
    "centripetal_force",
    "flywheel_energy",
    "disk_inertia",
    "solid_sphere_inertia",
    "parallel_axis",
    "lagrange_accelerations",
    "simulate_lagrangian",
    "gear_ratio",
    "gear_output_speed",
    "gear_output_torque",
    "belt_speed",
    "belt_tension_ratio",
    "lewis_form_factor",
    "lewis_bending_stress",
    "barth_velocity_factor",
    "screw_torque_raise",
    "petroff_friction_torque",
    "sommerfeld_number",
    "long_bearing_eccentricity",
    "bearing_l10_life",
    "shaft_von_mises_stress",
    "shaft_max_shear_stress",
    "shaft_diameter",
    "shaft_diameter_goodman",
    "spring_index",
    "wahl_factor",
    "bergstrasser_factor",
    "spring_rate",
    "spring_shear_stress",
    "bolt_stiffness",
    "member_stiffness",
    "bolted_joint",
    "preload_torque",
    "slider_crank_displacement",
    "slider_crank_kinematics",
    "fourbar_grashof",
    "fourbar_output_angle",
    "fourbar_kinematics",
    "freudenstein_coefficients",
    "freudenstein_link_lengths",
    "dh_transform",
    "forward_kinematics",
    "natural_frequency",
    "natural_frequency_hz",
    "damping_ratio",
    "damped_frequency",
    "transmissibility",
    "magnification_factor",
    "phase_angle",
    "simulate_sdof",
    "tuned_mass_damper",
    "holzer_residual",
    "holzer_frequencies",
    "normal_modes",
    "modal_response",
    "harmonic_response",
    "rayleigh_frequency",
    "dunkerley_frequency",
    "rotating_unbalance_force",
    "jeffcott_response",
    "unbalance",
    "balance_rotating_masses",
    "cam_follower_motion",
]

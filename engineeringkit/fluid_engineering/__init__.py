"""engineeringkit.fluid_engineering: hydrostatics, pipe flow, similarity, open channels, transients, pumps and networks.

Stevin's hydrostatic pressure and the force on submerged surfaces,
Archimedes's buoyancy; Bernoulli's equation and continuity, with Pitot and
Venturi meters; Reynolds's regimes, Hagen-Poiseuille laminar flow,
Darcy-Weisbach losses and the turbulent friction factors of Nikuradse,
Colebrook and Swamee-Jain; Froude, Euler, Weber and Mach numbers and
Buckingham's Pi theorem by exact rational row reduction; Manning's
open-channel formula; Joukowsky's water hammer with Korteweg's wave speed;
pump power, specific speed, NPSH and Thoma's cavitation parameter; and
Hardy Cross's balancing of looped pipe networks. Models are hand-written
closed forms; SciPy supplies root bracketing.
"""

from engineeringkit import __version__
from engineeringkit.fluid_engineering.core.base import NetworkResult, PiGroups
from engineeringkit.fluid_engineering.systems.bernoulli import (
    bernoulli_pressure,
    continuity_velocity,
    dynamic_pressure,
    velocity_head,
)
from engineeringkit.fluid_engineering.systems.hydrostatics import (
    buoyancy_force,
    center_of_pressure_depth,
    hydrostatic_force_plane,
    hydrostatic_pressure,
    pressure_head,
)
from engineeringkit.fluid_engineering.systems.measurement import (
    pitot_velocity,
    venturi_flow,
)
from engineeringkit.fluid_engineering.systems.networks import (
    hardy_cross,
)
from engineeringkit.fluid_engineering.systems.open_channel import (
    hydraulic_radius,
    manning_flow,
    manning_velocity,
)
from engineeringkit.fluid_engineering.systems.pipe_flow import (
    colebrook_friction_factor,
    darcy_head_loss,
    hydraulic_diameter,
    laminar_friction_factor,
    minor_head_loss,
    nikuradse_friction_factor,
    pipe_resistance,
    poiseuille_flow,
    poiseuille_velocity,
    reynolds_number,
    smooth_pipe_friction_factor,
    swamee_jain_friction_factor,
)
from engineeringkit.fluid_engineering.systems.pumps import (
    cavitation_number,
    npsh_available,
    pump_hydraulic_power,
    pump_input_power,
    pump_specific_speed,
    thoma_sigma,
)
from engineeringkit.fluid_engineering.systems.similarity import (
    buckingham_pi,
    euler_number,
    froude_number,
    mach_number,
    weber_number,
)
from engineeringkit.fluid_engineering.systems.transients import (
    critical_closure_time,
    pipe_wave_speed,
    water_hammer_pressure,
)
from engineeringkit.fluid_engineering.utils.constants import STANDARD_GRAVITY, G

__all__ = [
    "__version__",
    "PiGroups",
    "NetworkResult",
    "G",
    "STANDARD_GRAVITY",
    "hydrostatic_pressure",
    "pressure_head",
    "buoyancy_force",
    "hydrostatic_force_plane",
    "center_of_pressure_depth",
    "bernoulli_pressure",
    "continuity_velocity",
    "dynamic_pressure",
    "velocity_head",
    "pitot_velocity",
    "venturi_flow",
    "reynolds_number",
    "hydraulic_diameter",
    "poiseuille_flow",
    "poiseuille_velocity",
    "laminar_friction_factor",
    "darcy_head_loss",
    "minor_head_loss",
    "pipe_resistance",
    "smooth_pipe_friction_factor",
    "nikuradse_friction_factor",
    "colebrook_friction_factor",
    "swamee_jain_friction_factor",
    "froude_number",
    "weber_number",
    "euler_number",
    "mach_number",
    "buckingham_pi",
    "hydraulic_radius",
    "manning_velocity",
    "manning_flow",
    "water_hammer_pressure",
    "pipe_wave_speed",
    "critical_closure_time",
    "pump_hydraulic_power",
    "pump_input_power",
    "pump_specific_speed",
    "cavitation_number",
    "npsh_available",
    "thoma_sigma",
    "hardy_cross",
]

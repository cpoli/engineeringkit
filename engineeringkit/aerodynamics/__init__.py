"""engineeringkit.aerodynamics: from dynamic pressure to finite wings and compressible flow.

Force coefficients and dynamic pressure; Kutta-Joukowski circulation lift;
thin-airfoil theory (lift, pitching moment, aerodynamic centre), NACA
four-digit sections and the Hess-Smith panel method; Prandtl's lifting line,
with Glauert's Fourier solution for any planform, induced drag and the
parabolic drag polar; static longitudinal stability; climb, turn, take-off
and jet and propeller Breguet range; Reynolds-number similarity and
Blasius/Prandtl flat-plate boundary layers; perfect-gas sound speed,
isentropic stagnation states, the Mach cone, Prandtl-Glauert,
Rankine-Hugoniot normal shocks, oblique shocks, Prandtl-Meyer expansions,
the de Laval area-Mach relation, and Fanno and Rayleigh duct flow; and the
International Standard Atmosphere. Models are hand-written closed forms or
small linear solves (the formula is the pedagogical subject); SciPy is used
only for quadrature and the inverse Prandtl-Meyer root.
"""

from engineeringkit import __version__
from engineeringkit.aerodynamics.core.base import (
    Airfoil,
    AtmosphereState,
    DuctFlowState,
    LiftingLineResult,
    ObliqueShockResult,
    PanelResult,
)
from engineeringkit.aerodynamics.systems.airfoils import naca4_airfoil
from engineeringkit.aerodynamics.systems.atmosphere import GAMMA_AIR, R_AIR, isa_atmosphere
from engineeringkit.aerodynamics.systems.boundary_layer import (
    boundary_layer_laminar,
    reynolds_number,
    skin_friction_laminar,
    skin_friction_turbulent,
)
from engineeringkit.aerodynamics.systems.compressible import (
    area_mach_ratio,
    mach_angle,
    mach_number,
    normal_shock_mach2,
    normal_shock_pressure_ratio,
    prandtl_glauert_cp,
    speed_of_sound,
    stagnation_pressure,
    stagnation_temperature,
)
from engineeringkit.aerodynamics.systems.duct_flow import fanno_choking_length, fanno_flow, rayleigh_choking_heat, rayleigh_flow
from engineeringkit.aerodynamics.systems.lift import (
    aerodynamic_center,
    drag,
    drag_coefficient,
    dynamic_pressure,
    kutta_joukowski_lift,
    lift,
    lift_coefficient,
    naca4_moment_coefficient,
    naca4_zero_lift_angle,
    stall_speed,
    thin_airfoil_cl,
)
from engineeringkit.aerodynamics.systems.lifting_line import lifting_line
from engineeringkit.aerodynamics.systems.panels import panel_method
from engineeringkit.aerodynamics.systems.performance import (
    breguet_range_propeller,
    climb_angle,
    load_factor,
    rate_of_climb,
    takeoff_ground_roll,
    turn_radius,
    turn_rate,
)
from engineeringkit.aerodynamics.systems.stability import (
    downwash_gradient,
    neutral_point,
    pitching_moment_slope,
    static_margin,
    tail_volume_ratio,
)
from engineeringkit.aerodynamics.systems.supersonic import (
    inverse_prandtl_meyer,
    oblique_shock,
    oblique_shock_angle,
    oblique_shock_deflection,
    oblique_shock_max_deflection,
    prandtl_meyer_angle,
    prandtl_meyer_expansion,
)
from engineeringkit.aerodynamics.systems.wings import (
    aspect_ratio,
    best_glide_cl,
    breguet_range,
    drag_polar,
    finite_wing_cl,
    finite_wing_lift_slope,
    glide_ratio,
    induced_drag_coefficient,
    oswald_efficiency_estimate,
)
from engineeringkit.aerodynamics.utils.gas import isentropic_factor

__all__ = [
    "__version__",
    "Airfoil",
    "AtmosphereState",
    "ObliqueShockResult",
    "DuctFlowState",
    "LiftingLineResult",
    "PanelResult",
    "dynamic_pressure",
    "lift",
    "drag",
    "lift_coefficient",
    "drag_coefficient",
    "stall_speed",
    "kutta_joukowski_lift",
    "thin_airfoil_cl",
    "naca4_zero_lift_angle",
    "naca4_moment_coefficient",
    "aerodynamic_center",
    "lifting_line",
    "panel_method",
    "oblique_shock_deflection",
    "oblique_shock_max_deflection",
    "oblique_shock_angle",
    "oblique_shock",
    "prandtl_meyer_angle",
    "inverse_prandtl_meyer",
    "prandtl_meyer_expansion",
    "fanno_flow",
    "fanno_choking_length",
    "rayleigh_flow",
    "rayleigh_choking_heat",
    "tail_volume_ratio",
    "downwash_gradient",
    "neutral_point",
    "pitching_moment_slope",
    "static_margin",
    "rate_of_climb",
    "climb_angle",
    "load_factor",
    "turn_radius",
    "turn_rate",
    "takeoff_ground_roll",
    "breguet_range_propeller",
    "naca4_airfoil",
    "aspect_ratio",
    "induced_drag_coefficient",
    "drag_polar",
    "finite_wing_lift_slope",
    "finite_wing_cl",
    "oswald_efficiency_estimate",
    "glide_ratio",
    "best_glide_cl",
    "breguet_range",
    "reynolds_number",
    "skin_friction_laminar",
    "skin_friction_turbulent",
    "boundary_layer_laminar",
    "speed_of_sound",
    "mach_number",
    "mach_angle",
    "stagnation_temperature",
    "stagnation_pressure",
    "prandtl_glauert_cp",
    "normal_shock_mach2",
    "normal_shock_pressure_ratio",
    "area_mach_ratio",
    "isentropic_factor",
    "isa_atmosphere",
    "R_AIR",
    "GAMMA_AIR",
]

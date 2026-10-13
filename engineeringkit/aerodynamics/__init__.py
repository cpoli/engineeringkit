"""engineeringkit.aerodynamics: from dynamic pressure to finite wings and compressible flow.

Force coefficients and dynamic pressure; Kutta-Joukowski circulation lift;
thin-airfoil theory and NACA four-digit sections; Prandtl's lifting line,
induced drag, the parabolic drag polar and Breguet range; Reynolds-number
similarity and Blasius/Prandtl flat-plate boundary layers; perfect-gas sound
speed, isentropic stagnation states, the Mach cone, Prandtl-Glauert,
Rankine-Hugoniot normal shocks and the de Laval area-Mach relation; and the
International Standard Atmosphere. Every model is a hand-written closed-form
relation (the formula is the pedagogical subject); SciPy is used only for the
thin-airfoil quadrature.
"""

from engineeringkit import __version__
from engineeringkit.aerodynamics.core.base import Airfoil, AtmosphereState
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
from engineeringkit.aerodynamics.systems.lift import (
    drag,
    drag_coefficient,
    dynamic_pressure,
    kutta_joukowski_lift,
    lift,
    lift_coefficient,
    naca4_zero_lift_angle,
    stall_speed,
    thin_airfoil_cl,
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
    "dynamic_pressure",
    "lift",
    "drag",
    "lift_coefficient",
    "drag_coefficient",
    "stall_speed",
    "kutta_joukowski_lift",
    "thin_airfoil_cl",
    "naca4_zero_lift_angle",
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

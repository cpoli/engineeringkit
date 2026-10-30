"""engineeringkit.heat_transfer: conduction, convection, radiation, phase change and heat exchangers.

Fourier's law and thermal-resistance networks, with composite walls, radial
shells, critical insulation and Harper-Brown fins (adiabatic or convecting
tip); finite-difference conduction (explicit, implicit, Crank-Nicolson, and
steady 2-D) and conduction shape factors; Newton cooling, lumped and
semi-infinite transients and Heisler's plane-wall, cylinder and sphere
series; forced convection over flat plates, cylinders and tube banks and in
laminar tubes, Dittus-Boelter and Churchill-Chu, and the Chilton-Colburn
analogy; blackbody laws (Planck, Stefan-Boltzmann, Wien), view factors and
Oppenheim's gray-surface network with shields; Neumann's solution of the
Stefan problem, Nusselt film condensation, Rohsenow's nucleate boiling and
Zuber's critical heat flux; and LMTD and effectiveness-NTU exchangers in
every common arrangement, with the LMTD correction factor. Models are
hand-written closed forms; SciPy supplies erf, Bessel functions, root
bracketing and the banded and sparse solvers."""

from engineeringkit import __version__
from engineeringkit.heat_transfer.core.base import Conduction2DResult, ConductionResult, EnclosureResult, WallResult
from engineeringkit.heat_transfer.systems.conduction import (
    composite_wall,
    conduction_1d,
    critical_insulation_radius_cylinder,
    critical_insulation_radius_sphere,
    cylinder_resistance,
    fin_efficiency_adiabatic,
    fin_heat_rate_adiabatic,
    fin_heat_rate_convective_tip,
    fin_parameter,
    fourier_heat_rate,
    overall_u,
    plane_wall_resistance,
    sphere_resistance,
)
from engineeringkit.heat_transfer.systems.convection import (
    churchill_chu_nusselt,
    convection_heat_rate,
    convection_resistance,
    dittus_boelter_nusselt,
    rayleigh_number,
)
from engineeringkit.heat_transfer.systems.exchangers import (
    exchanger_effectiveness_counterflow,
    exchanger_effectiveness_crossflow,
    exchanger_effectiveness_parallel,
    exchanger_effectiveness_shell_tube,
    heat_exchanger_rate,
    lmtd_correction_factor,
    log_mean_temperature_difference,
)
from engineeringkit.heat_transfer.systems.forced_convection import (
    chilton_colburn_mass_transfer,
    churchill_bernstein_nusselt,
    entry_length,
    flat_plate_nusselt,
    hausen_nusselt,
    laminar_tube_nusselt,
    lewis_number,
    tube_bank_max_velocity,
    zukauskas_nusselt,
)
from engineeringkit.heat_transfer.systems.multidimensional import (
    conduction_2d,
    shape_factor_buried_cylinder,
    shape_factor_buried_sphere,
    shape_factor_eccentric_cylinders,
    shape_factor_parallel_cylinders,
)
from engineeringkit.heat_transfer.systems.phase_change import (
    neumann_stefan_front,
    nusselt_film_condensation,
    rohsenow_boiling_heat_flux,
    stefan_number,
    zuber_critical_heat_flux,
)
from engineeringkit.heat_transfer.systems.radiation import (
    blackbody_peak_wavelength,
    planck_spectral_emissive_power,
    radiation_coefficient,
    radiation_enclosure,
    radiation_heat_rate,
    radiation_shield_heat_flux,
    stefan_boltzmann_emissive_power,
    two_surface_radiation,
    view_factor_coaxial_disks,
    view_factor_parallel_plates,
    view_factor_reciprocity,
)
from engineeringkit.heat_transfer.systems.transient import (
    biot_number,
    cylinder_transient,
    fourier_number,
    lumped_temperature,
    plane_wall_transient,
    semi_infinite_temperature,
    sphere_transient,
)
from engineeringkit.heat_transfer.utils.constants import STEFAN_BOLTZMANN, WIEN_DISPLACEMENT

__all__ = [
    "__version__",
    "WallResult",
    "ConductionResult",
    "Conduction2DResult",
    "EnclosureResult",
    "STEFAN_BOLTZMANN",
    "WIEN_DISPLACEMENT",
    "fourier_heat_rate",
    "plane_wall_resistance",
    "cylinder_resistance",
    "sphere_resistance",
    "overall_u",
    "critical_insulation_radius_cylinder",
    "critical_insulation_radius_sphere",
    "composite_wall",
    "fin_parameter",
    "fin_efficiency_adiabatic",
    "fin_heat_rate_adiabatic",
    "conduction_1d",
    "convection_heat_rate",
    "convection_resistance",
    "dittus_boelter_nusselt",
    "rayleigh_number",
    "churchill_chu_nusselt",
    "biot_number",
    "fourier_number",
    "lumped_temperature",
    "semi_infinite_temperature",
    "plane_wall_transient",
    "radiation_heat_rate",
    "radiation_coefficient",
    "stefan_boltzmann_emissive_power",
    "blackbody_peak_wavelength",
    "planck_spectral_emissive_power",
    "view_factor_parallel_plates",
    "view_factor_coaxial_disks",
    "stefan_number",
    "neumann_stefan_front",
    "nusselt_film_condensation",
    "zuber_critical_heat_flux",
    "log_mean_temperature_difference",
    "heat_exchanger_rate",
    "exchanger_effectiveness_parallel",
    "exchanger_effectiveness_counterflow",
    "fin_heat_rate_convective_tip",
    "conduction_2d",
    "shape_factor_buried_sphere",
    "shape_factor_buried_cylinder",
    "shape_factor_eccentric_cylinders",
    "shape_factor_parallel_cylinders",
    "flat_plate_nusselt",
    "churchill_bernstein_nusselt",
    "tube_bank_max_velocity",
    "zukauskas_nusselt",
    "laminar_tube_nusselt",
    "entry_length",
    "hausen_nusselt",
    "lewis_number",
    "chilton_colburn_mass_transfer",
    "exchanger_effectiveness_crossflow",
    "exchanger_effectiveness_shell_tube",
    "lmtd_correction_factor",
    "rohsenow_boiling_heat_flux",
    "view_factor_reciprocity",
    "two_surface_radiation",
    "radiation_shield_heat_flux",
    "radiation_enclosure",
    "cylinder_transient",
    "sphere_transient",
]

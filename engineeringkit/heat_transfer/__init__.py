"""engineeringkit.heat_transfer: conduction, convection, radiation, phase change and heat exchangers.

Fourier's law and thermal-resistance networks, with composite walls, radial
shells, critical insulation and Harper-Brown fins; finite-difference
conduction (explicit, implicit, Crank-Nicolson); Newton cooling, lumped and
semi-infinite transients and Heisler's plane-wall series; Dittus-Boelter
and Churchill-Chu convection correlations; blackbody laws (Planck,
Stefan-Boltzmann, Wien), gray-surface exchange and view factors; Neumann's
solution of the Stefan problem, Nusselt film condensation and Zuber's
critical heat flux; and LMTD and effectiveness-NTU exchangers. Models are
hand-written closed forms; SciPy supplies erf, root bracketing and the
banded solver.
"""

from engineeringkit import __version__
from engineeringkit.heat_transfer.core.base import ConductionResult, WallResult
from engineeringkit.heat_transfer.systems.conduction import (
    composite_wall,
    conduction_1d,
    critical_insulation_radius_cylinder,
    critical_insulation_radius_sphere,
    cylinder_resistance,
    fin_efficiency_adiabatic,
    fin_heat_rate_adiabatic,
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
    exchanger_effectiveness_parallel,
    heat_exchanger_rate,
    log_mean_temperature_difference,
)
from engineeringkit.heat_transfer.systems.phase_change import (
    neumann_stefan_front,
    nusselt_film_condensation,
    stefan_number,
    zuber_critical_heat_flux,
)
from engineeringkit.heat_transfer.systems.radiation import (
    blackbody_peak_wavelength,
    planck_spectral_emissive_power,
    radiation_coefficient,
    radiation_heat_rate,
    stefan_boltzmann_emissive_power,
    view_factor_coaxial_disks,
    view_factor_parallel_plates,
)
from engineeringkit.heat_transfer.systems.transient import (
    biot_number,
    fourier_number,
    lumped_temperature,
    plane_wall_transient,
    semi_infinite_temperature,
)
from engineeringkit.heat_transfer.utils.constants import STEFAN_BOLTZMANN, WIEN_DISPLACEMENT

__all__ = [
    "__version__",
    "WallResult",
    "ConductionResult",
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
]

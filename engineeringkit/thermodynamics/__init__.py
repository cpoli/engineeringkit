"""engineeringkit.thermodynamics: gases, the two laws, cycles, refrigeration, exergy, steam and moist air.

Ideal and van der Waals gases with the Joule-Thomson coefficient; Joule's
first law and Clausius's entropy; Carnot, Otto, Diesel, Brayton and Rankine
cycles and the Curzon-Ahlborn efficiency at maximum power; refrigerator and
heat-pump COPs; Gouy-Stodola lost work and the exergy of heat; vapour
pressure by Clausius-Clapeyron and Antoine; water and steam properties from
IAPWS-IF97 (regions 1 to 4); and Carrier's psychrometrics. Models are
hand-written closed forms; SciPy's ``brentq`` inverts IF97 for entropy and
enthalpy.
"""

from engineeringkit import __version__
from engineeringkit.thermodynamics.core.base import SteamState
from engineeringkit.thermodynamics.systems.cycles import (
    OttoCycle,
    brayton_efficiency,
    carnot_efficiency,
    curzon_ahlborn_efficiency,
    diesel_efficiency,
    heat_engine_work,
    otto_efficiency,
    rankine_ideal_efficiency,
    thermal_efficiency,
)
from engineeringkit.thermodynamics.systems.exergy import exergy_destruction, exergy_heat
from engineeringkit.thermodynamics.systems.gases import (
    R_UNIVERSAL,
    ideal_gas_pressure,
    ideal_gas_temperature,
    joule_thomson_coefficient_vdw,
    van_der_waals_critical_point,
    van_der_waals_pressure,
)
from engineeringkit.thermodynamics.systems.laws import (
    clausius_inequality,
    entropy_change_isothermal_ideal_gas,
    first_law_heat,
    isentropic_pressure_volume,
    isentropic_temperature_ratio,
    reversible_isothermal_work_ideal_gas,
    specific_enthalpy_ideal_gas,
    specific_internal_energy_ideal_gas,
)
from engineeringkit.thermodynamics.systems.psychrometrics import MOLAR_MASS_RATIO, dew_point, humidity_ratio, relative_humidity
from engineeringkit.thermodynamics.systems.refrigeration import (
    carnot_cop_heat_pump,
    carnot_cop_refrigerator,
    heat_pump_cop,
    refrigeration_cop,
)
from engineeringkit.thermodynamics.systems.saturation import antoine_vapor_pressure, clausius_clapeyron_pressure
from engineeringkit.thermodynamics.systems.steam import saturation_pressure, saturation_temperature, steam_properties

__all__ = [
    "__version__",
    "SteamState",
    "R_UNIVERSAL",
    "ideal_gas_pressure",
    "ideal_gas_temperature",
    "van_der_waals_pressure",
    "van_der_waals_critical_point",
    "joule_thomson_coefficient_vdw",
    "first_law_heat",
    "specific_internal_energy_ideal_gas",
    "specific_enthalpy_ideal_gas",
    "entropy_change_isothermal_ideal_gas",
    "reversible_isothermal_work_ideal_gas",
    "isentropic_pressure_volume",
    "isentropic_temperature_ratio",
    "clausius_inequality",
    "carnot_efficiency",
    "curzon_ahlborn_efficiency",
    "heat_engine_work",
    "thermal_efficiency",
    "otto_efficiency",
    "OttoCycle",
    "diesel_efficiency",
    "brayton_efficiency",
    "rankine_ideal_efficiency",
    "carnot_cop_refrigerator",
    "carnot_cop_heat_pump",
    "refrigeration_cop",
    "heat_pump_cop",
    "exergy_heat",
    "exergy_destruction",
    "clausius_clapeyron_pressure",
    "antoine_vapor_pressure",
    "saturation_pressure",
    "saturation_temperature",
    "steam_properties",
    "MOLAR_MASS_RATIO",
    "humidity_ratio",
    "relative_humidity",
    "dew_point",
]

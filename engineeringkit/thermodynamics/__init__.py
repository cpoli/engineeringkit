"""engineeringkit.thermodynamics: gases, the two laws, cycles, refrigeration, exergy, steam, moist air and combustion.

Ideal, van der Waals and Peng-Robinson gases, the compressibility factor and
Pitzer's acentric factor, the Joule-Thomson coefficient, Dalton's ideal-gas
mixtures and Maxwell's relations; Joule's first law for closed systems and
steady-flow devices, with isentropic efficiencies, and Clausius's entropy;
Carnot, Otto, Diesel, dual, Brayton, Stirling and Rankine cycles, Rankine
plants with reheat and feedwater heating, gas turbines with regeneration,
intercooling and reheat, and the Curzon-Ahlborn efficiency at maximum power;
refrigerator and heat-pump COPs and the vapour-compression cycle on R-134a;
Gouy-Stodola lost work and the exergy of heat; vapour pressure by
Clausius-Clapeyron and Antoine; water and steam properties from IAPWS-IF97
(all five regions); Carrier's psychrometrics and air-conditioning processes;
and combustion stoichiometry, heating values and the adiabatic flame
temperature. Models are hand-written closed forms; SciPy's ``brentq``
inverts the property equations.
"""

from engineeringkit import __version__
from engineeringkit.thermodynamics.core.base import (
    BraytonResult,
    CubicFluid,
    RankineResult,
    RefrigerantState,
    SteamState,
    VaporCompressionResult,
)
from engineeringkit.thermodynamics.systems.combustion import (
    MOLAR_MASS_AIR,
    adiabatic_flame_temperature,
    combustion_products,
    fuel_molar_mass,
    heating_value,
    stoichiometric_air_fuel_ratio,
)
from engineeringkit.thermodynamics.systems.cycles import (
    OttoCycle,
    brayton_efficiency,
    carnot_efficiency,
    curzon_ahlborn_efficiency,
    diesel_efficiency,
    dual_efficiency,
    heat_engine_work,
    otto_efficiency,
    rankine_ideal_efficiency,
    stirling_efficiency,
    thermal_efficiency,
)
from engineeringkit.thermodynamics.systems.exergy import exergy_destruction, exergy_heat
from engineeringkit.thermodynamics.systems.gases import (
    R_UNIVERSAL,
    compressibility_factor,
    ideal_gas_pressure,
    ideal_gas_temperature,
    joule_thomson_coefficient_vdw,
    pitzer_acentric_factor,
    pitzer_virial_compressibility,
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
from engineeringkit.thermodynamics.systems.mixtures import (
    entropy_of_mixing,
    mass_fractions,
    mixture_molar_mass,
    mole_fractions,
    partial_pressures,
)
from engineeringkit.thermodynamics.systems.peng_robinson import (
    R134A,
    peng_robinson_compressibility,
    peng_robinson_pressure,
    peng_robinson_saturation_pressure,
    refrigerant_properties,
)
from engineeringkit.thermodynamics.systems.power_cycles import (
    brayton_cycle,
    rankine_cycle,
)
from engineeringkit.thermodynamics.systems.psychrometrics import (
    MOLAR_MASS_RATIO,
    adiabatic_mixing,
    adiabatic_saturation_temperature,
    cooling_coil,
    cooling_tower,
    dew_point,
    humidity_ratio,
    moist_air_enthalpy,
    relative_humidity,
)
from engineeringkit.thermodynamics.systems.refrigeration import (
    carnot_cop_heat_pump,
    carnot_cop_refrigerator,
    heat_pump_cop,
    refrigeration_cop,
    vapor_compression_cycle,
)
from engineeringkit.thermodynamics.systems.relations import (
    cp_minus_cv,
    joule_thomson_coefficient,
)
from engineeringkit.thermodynamics.systems.saturation import antoine_vapor_pressure, clausius_clapeyron_pressure
from engineeringkit.thermodynamics.systems.steady_flow import (
    compression_outlet_enthalpy,
    compressor_isentropic_efficiency,
    expansion_outlet_enthalpy,
    mixing_chamber,
    nozzle_exit_velocity,
    nozzle_isentropic_efficiency,
    steady_flow_work,
    throttle_steam,
    turbine_isentropic_efficiency,
)
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
    "RankineResult",
    "BraytonResult",
    "CubicFluid",
    "RefrigerantState",
    "VaporCompressionResult",
    "dual_efficiency",
    "stirling_efficiency",
    "compressibility_factor",
    "pitzer_acentric_factor",
    "pitzer_virial_compressibility",
    "moist_air_enthalpy",
    "adiabatic_saturation_temperature",
    "adiabatic_mixing",
    "cooling_coil",
    "cooling_tower",
    "vapor_compression_cycle",
    "steady_flow_work",
    "nozzle_exit_velocity",
    "turbine_isentropic_efficiency",
    "compressor_isentropic_efficiency",
    "nozzle_isentropic_efficiency",
    "expansion_outlet_enthalpy",
    "compression_outlet_enthalpy",
    "mixing_chamber",
    "throttle_steam",
    "rankine_cycle",
    "brayton_cycle",
    "R134A",
    "peng_robinson_pressure",
    "peng_robinson_compressibility",
    "peng_robinson_saturation_pressure",
    "refrigerant_properties",
    "mole_fractions",
    "mass_fractions",
    "mixture_molar_mass",
    "partial_pressures",
    "entropy_of_mixing",
    "cp_minus_cv",
    "joule_thomson_coefficient",
    "MOLAR_MASS_AIR",
    "fuel_molar_mass",
    "stoichiometric_air_fuel_ratio",
    "combustion_products",
    "heating_value",
    "adiabatic_flame_temperature",
]

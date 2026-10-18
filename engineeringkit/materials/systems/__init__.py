"""Concrete materials models: elasticity, yielding, composites, fatigue, fracture, creep and plasticity."""

from engineeringkit.materials.systems.composites import (
    density_rule_of_mixtures,
    rule_of_mixtures_E_longitudinal,
    rule_of_mixtures_E_transverse,
)
from engineeringkit.materials.systems.creep import arrhenius_creep_rate, larson_miller_parameter, larson_miller_rupture_time
from engineeringkit.materials.systems.elasticity import (
    bulk_modulus,
    constrained_thermal_stress,
    hooke_stress,
    lame_parameters,
    normal_strain,
    normal_stress,
    resilience_modulus,
    shear_modulus,
    thermal_strain,
    toughness,
)
from engineeringkit.materials.systems.fatigue import (
    basquin_cycles,
    coffin_manson_cycles,
    goodman_safety_factor,
    miner_damage,
    strain_life_amplitude,
)
from engineeringkit.materials.systems.fracture import (
    crack_growth_life,
    critical_crack_length,
    griffith_fracture_stress,
    paris_crack_growth,
    stress_intensity,
)
from engineeringkit.materials.systems.plasticity import hall_petch, ideal_shear_strength, peierls_stress, ramberg_osgood_strain
from engineeringkit.materials.systems.yielding import factor_of_safety, tresca_equivalent, von_mises_plane_stress

__all__ = [
    "normal_stress",
    "normal_strain",
    "hooke_stress",
    "shear_modulus",
    "bulk_modulus",
    "lame_parameters",
    "thermal_strain",
    "constrained_thermal_stress",
    "resilience_modulus",
    "toughness",
    "von_mises_plane_stress",
    "tresca_equivalent",
    "factor_of_safety",
    "rule_of_mixtures_E_longitudinal",
    "rule_of_mixtures_E_transverse",
    "density_rule_of_mixtures",
    "basquin_cycles",
    "goodman_safety_factor",
    "coffin_manson_cycles",
    "strain_life_amplitude",
    "miner_damage",
    "griffith_fracture_stress",
    "stress_intensity",
    "critical_crack_length",
    "paris_crack_growth",
    "crack_growth_life",
    "arrhenius_creep_rate",
    "larson_miller_parameter",
    "larson_miller_rupture_time",
    "ideal_shear_strength",
    "peierls_stress",
    "hall_petch",
    "ramberg_osgood_strain",
]

"""engineeringkit.materials: elasticity, yielding, fracture, fatigue, creep, plasticity and composites.

Young's modulus, Lamé's isotropic constants and Duhamel's thermal stress;
Tresca's and von Mises's yield criteria; the Voigt and Reuss bounds on a
composite; Basquin's S-N law, Goodman's mean-stress line, the Coffin-Manson
strain-life law and Palmgren-Miner damage; Griffith's and Irwin's fracture
mechanics and Paris's crack growth; Norton creep and the Larson-Miller
parameter; Frenkel's ideal strength, the Peierls stress, Hall-Petch
strengthening and the Ramberg-Osgood curve. Models are hand-written closed
forms; SciPy's trapezoid rules integrate measured curves and crack growth.
"""

from engineeringkit import __version__
from engineeringkit.materials.core.base import CrackGrowthResult, StressStrainCurve
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
    "__version__",
    "StressStrainCurve",
    "CrackGrowthResult",
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

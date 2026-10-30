"""engineeringkit.structural_engineering: members, sections, beams, columns, stress, vessels, trusses, frames, cables and arches.

Hooke's law and Poisson's ratio for axial members, Coulomb's torsion and
Bredt's thin-walled tubes; section properties of rectangles and of I, T,
channel, angle, box and transformed sections; Navier's bending theory,
unsymmetrical bending, Jourawski's shear stress and the shear centre;
Macaulay's beams on any supports and superposition, Clapeyron's
three-moment equation and Winkler's influence lines; Euler's and
Engesser's column loads, the secant formula and Johnson's parabola; Mohr's
circle, three-dimensional principal stresses and strain rosettes; thin and
thick-walled pressure vessels; Castigliano's theorem, the unit-load method
and Maxwell-Betti reciprocity; Maxwell's count, Whipple's method of joints
and the direct stiffness method for trusses; plastic collapse loads;
Maney's slope-deflection, Hardy Cross's moment distribution and plane
frames; cables and three-hinged arches. Models are hand-written; NumPy's
linear algebra and SciPy's Simpson rule and Brent root finder are used
where they are the standard tool.
"""

from engineeringkit import __version__
from engineeringkit.structural_engineering.core.base import (
    ArchResult,
    BeamResult,
    FrameResult,
    JointMethodResult,
    MomentDistributionResult,
    SectionProperties,
    SlopeDeflectionResult,
    TrussResult,
)
from engineeringkit.structural_engineering.systems.beams import (
    Beam,
    SimplySupportedBeam,
    bending_stress,
    influence_line,
    neutral_axis_angle,
    superpose_beam_results,
    three_moment_equation,
    unsymmetrical_bending_stress,
)
from engineeringkit.structural_engineering.systems.cables import catenary_cable, parabolic_cable, three_hinged_arch
from engineeringkit.structural_engineering.systems.energy import (
    castigliano_deflection,
    maxwell_betti_reciprocity,
    strain_energy_bending_constant_moment,
    unit_load_beam_deflection,
    unit_load_truss_deflection,
)
from engineeringkit.structural_engineering.systems.frames import PlanarFrame, moment_distribution, slope_deflection
from engineeringkit.structural_engineering.systems.members import (
    axial_deformation,
    hooke_stress,
    poisson_lateral_strain,
    strain_energy_axial,
    thermal_strain,
    thin_walled_torsion_angle,
    thin_walled_torsion_shear_stress,
    torsion_angle,
    torsion_shear_stress,
)
from engineeringkit.structural_engineering.systems.plasticity import plastic_collapse_load
from engineeringkit.structural_engineering.systems.sections import (
    angle_section,
    box_section,
    channel_section,
    composite_section,
    i_section,
    parallel_axis,
    plastic_section_modulus,
    rectangular_second_moment,
    rectangular_section_modulus,
    t_section,
)
from engineeringkit.structural_engineering.systems.shear import channel_shear_center, shear_flow, transverse_shear_stress
from engineeringkit.structural_engineering.systems.stability import (
    effective_length_factor,
    euler_buckling_load,
    johnson_critical_stress,
    radius_of_gyration,
    secant_formula_load,
    secant_formula_stress,
    slenderness_ratio,
    tangent_modulus_load,
)
from engineeringkit.structural_engineering.systems.stress import (
    mohr_circle,
    plane_stress_from_strains,
    principal_stresses_3d,
    strain_rosette,
    stress_transformation,
)
from engineeringkit.structural_engineering.systems.trusses import PlanarTruss, method_of_joints, truss_determinacy
from engineeringkit.structural_engineering.systems.vessels import thick_wall_cylinder_stresses, thin_wall_vessel_stresses

__all__ = [
    "__version__",
    "BeamResult",
    "TrussResult",
    "JointMethodResult",
    "MomentDistributionResult",
    "FrameResult",
    "SectionProperties",
    "SlopeDeflectionResult",
    "ArchResult",
    "hooke_stress",
    "axial_deformation",
    "strain_energy_axial",
    "thermal_strain",
    "poisson_lateral_strain",
    "torsion_shear_stress",
    "torsion_angle",
    "thin_walled_torsion_shear_stress",
    "thin_walled_torsion_angle",
    "rectangular_second_moment",
    "rectangular_section_modulus",
    "plastic_section_modulus",
    "parallel_axis",
    "composite_section",
    "i_section",
    "t_section",
    "channel_section",
    "angle_section",
    "box_section",
    "bending_stress",
    "unsymmetrical_bending_stress",
    "neutral_axis_angle",
    "transverse_shear_stress",
    "shear_flow",
    "channel_shear_center",
    "SimplySupportedBeam",
    "Beam",
    "superpose_beam_results",
    "three_moment_equation",
    "influence_line",
    "euler_buckling_load",
    "radius_of_gyration",
    "slenderness_ratio",
    "tangent_modulus_load",
    "effective_length_factor",
    "secant_formula_stress",
    "secant_formula_load",
    "johnson_critical_stress",
    "mohr_circle",
    "stress_transformation",
    "principal_stresses_3d",
    "strain_rosette",
    "plane_stress_from_strains",
    "thin_wall_vessel_stresses",
    "thick_wall_cylinder_stresses",
    "strain_energy_bending_constant_moment",
    "castigliano_deflection",
    "unit_load_truss_deflection",
    "unit_load_beam_deflection",
    "maxwell_betti_reciprocity",
    "truss_determinacy",
    "method_of_joints",
    "PlanarTruss",
    "plastic_collapse_load",
    "slope_deflection",
    "moment_distribution",
    "PlanarFrame",
    "parabolic_cable",
    "catenary_cable",
    "three_hinged_arch",
]

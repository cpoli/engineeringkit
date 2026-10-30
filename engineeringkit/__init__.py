"""engineeringkit: inspectable computational engineering for learning and teaching.

Import as ``ek`` by convention::

    import engineeringkit as ek
    ek.structural_engineering.SimplySupportedBeam(5.0, 200e9, 8e-6).solve_point_load(1e4)

Twelve engineering domains share SI-unit conventions, explicit input
validation (``engineeringkit._validation``), ODE integrators, unit
conversion, constants and solver-independent visualizers.
"""

# Defined before the subpackage imports: each subpackage re-exports it via
# ``from engineeringkit import __version__`` while this module is still loading.
__version__ = "0.2.0"

from engineeringkit import (  # noqa: E402
    aerodynamics,
    communications,
    constants,
    control_systems,
    electrical_engineering,
    fluid_engineering,
    heat_transfer,
    integrators,
    materials,
    mechanical_engineering,
    reliability,
    signal_processing,
    structural_engineering,
    thermodynamics,
    units,
)

__all__ = [
    "aerodynamics",
    "communications",
    "constants",
    "control_systems",
    "electrical_engineering",
    "fluid_engineering",
    "heat_transfer",
    "integrators",
    "materials",
    "mechanical_engineering",
    "reliability",
    "signal_processing",
    "structural_engineering",
    "thermodynamics",
    "units",
]

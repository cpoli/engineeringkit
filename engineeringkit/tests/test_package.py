"""Package-level smoke tests: engineeringkit imports cleanly and every declared
subpackage is importable and internally consistent."""

import importlib

import engineeringkit as ek

DOMAINS = [
    "aerodynamics",
    "communications",
    "control_systems",
    "electrical_engineering",
    "fluid_engineering",
    "heat_transfer",
    "materials",
    "mechanical_engineering",
    "reliability",
    "signal_processing",
    "structural_engineering",
    "thermodynamics",
]


def test_top_level_import_exposes_declared_subpackages():
    for name in ek.__all__:
        assert hasattr(ek, name), f"engineeringkit.{name} listed in __all__ but not importable"


def test_version_is_a_string():
    assert isinstance(ek.__version__, str)
    assert ek.__version__.count(".") == 2


def test_domain_versions_match_top_level_version():
    for name in DOMAINS:
        assert getattr(ek, name).__version__ == ek.__version__


def test_every_domain_follows_the_core_systems_utils_visualizers_layout():
    for name in DOMAINS:
        for part in ("core", "systems", "utils", "visualizers"):
            importlib.import_module(f"engineeringkit.{name}.{part}")


def test_domain_all_names_resolve_and_are_documented():
    for name in DOMAINS:
        module = getattr(ek, name)
        for symbol in module.__all__:
            obj = getattr(module, symbol)
            if callable(obj):
                assert obj.__doc__, f"engineeringkit.{name}.{symbol} has no docstring"

"""Lightweight explicit SI conversion helpers. Public APIs document SI base units."""

import math

_LENGTH = {"m": 1.0, "mm": 1e-3, "cm": 1e-2, "km": 1e3, "in": 0.0254, "ft": 0.3048}
_PRESSURE = {"Pa": 1.0, "kPa": 1e3, "MPa": 1e6, "bar": 1e5, "psi": 6894.757293168}


def convert(value: float, from_unit: str, to_unit: str) -> float:
    if not math.isfinite(float(value)):
        raise ValueError("value must be finite")
    tables = (_LENGTH, _PRESSURE)
    for t in tables:
        if from_unit in t and to_unit in t:
            return float(value) * t[from_unit] / t[to_unit]
    raise ValueError(f"incompatible or unknown units: {from_unit!r}, {to_unit!r}")

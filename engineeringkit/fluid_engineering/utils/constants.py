"""Constants used by the fluid-engineering models (SI)."""

from __future__ import annotations

__all__ = ["STANDARD_GRAVITY", "G", "LAMINAR_LIMIT", "TURBULENT_LIMIT"]

STANDARD_GRAVITY = 9.80665
"""Standard acceleration of gravity :math:`g` [m/s^2] (exact, CGPM 1901)."""

G = STANDARD_GRAVITY
"""Short alias of :data:`STANDARD_GRAVITY`, kept for backward compatibility."""

LAMINAR_LIMIT = 2300.0
"""Reynolds number below which pipe flow is treated as laminar [-]."""

TURBULENT_LIMIT = 4000.0
"""Reynolds number from which pipe flow is treated as fully turbulent [-]."""

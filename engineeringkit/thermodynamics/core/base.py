"""Result containers for engineeringkit.thermodynamics.

Most thermodynamic models return plain floats or arrays. A steam state
carries several properties at once, plus which IF97 region produced them,
so it gets a frozen dataclass with the usual ``summary()``.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["SteamState"]


@dataclass(frozen=True)
class SteamState:
    """Output of :func:`engineeringkit.thermodynamics.steam_properties`."""

    p: float
    """float: Pressure [Pa]."""

    T: float
    """float: Temperature [K]."""

    v: float
    """float: Specific volume [m^3/kg]."""

    h: float
    """float: Specific enthalpy [J/kg], zero for liquid at the triple point."""

    s: float
    """float: Specific entropy [J/(kg K)], zero for liquid at the triple point."""

    u: float
    """float: Specific internal energy [J/kg]."""

    cp: float
    """float: Isobaric specific heat [J/(kg K)]; NaN inside the two-phase dome."""

    x: float
    """float: Vapour quality [-] inside the dome; NaN for a single-phase state."""

    region: int
    """int: IF97 region, 1 (compressed liquid), 2 (vapour), 3 (near-critical) or 4 (wet mixture)."""

    def summary(self) -> dict:
        """All properties as a dict of floats, plus the region."""
        return {
            "p": self.p,
            "T": self.T,
            "v": self.v,
            "h": self.h,
            "s": self.s,
            "u": self.u,
            "cp": self.cp,
            "x": self.x,
            "region": self.region,
        }

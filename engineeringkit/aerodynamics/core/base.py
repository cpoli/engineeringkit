"""Result containers for engineeringkit.aerodynamics.

Most aerodynamic models here are closed-form algebraic relations that return
plain floats or arrays. The two models that produce structured output, the
NACA four-digit section generator and the standard atmosphere, return the
frozen dataclasses below so plotting helpers have a stable interface.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["Airfoil", "AtmosphereState"]


@dataclass(frozen=True)
class Airfoil:
    """Coordinates of an airfoil section on a unit chord.

    All coordinates are fractions of the chord ``c``; multiply by ``c`` [m]
    for physical dimensions.
    """

    name: str
    """str: Section designation, e.g. ``"NACA 2412"``."""

    x: np.ndarray
    """ndarray: Chordwise stations of the mean camber line, from 0 to 1."""

    camber: np.ndarray
    """ndarray: Mean-camber-line ordinate at each station."""

    thickness: np.ndarray
    """ndarray: Half-thickness at each station."""

    upper: np.ndarray
    """ndarray: Upper-surface points, shape ``(n, 2)``."""

    lower: np.ndarray
    """ndarray: Lower-surface points, shape ``(n, 2)``."""

    def summary(self) -> dict:
        """Maximum thickness and camber (fractions of chord) and their locations."""
        i_t, i_c = int(np.argmax(self.thickness)), int(np.argmax(self.camber))
        return {
            "name": self.name,
            "max_thickness": float(2 * self.thickness[i_t]),
            "max_thickness_at": float(self.x[i_t]),
            "max_camber": float(self.camber[i_c]),
            "max_camber_at": float(self.x[i_c]),
        }

    def plot(self, ax=None):
        """Draw the section outline and camber line; return ``(fig, ax)``."""
        from engineeringkit.aerodynamics.visualizers.plots import plot_airfoil

        ax = plot_airfoil(self, ax=ax)
        return ax.figure, ax


@dataclass(frozen=True)
class AtmosphereState:
    """Thermodynamic state of the standard atmosphere at a geometric altitude."""

    altitude: np.ndarray
    """ndarray: Altitude above mean sea level [m]."""

    temperature: np.ndarray
    """ndarray: Static temperature [K]."""

    pressure: np.ndarray
    """ndarray: Static pressure [Pa]."""

    density: np.ndarray
    """ndarray: Density [kg/m^3]."""

    speed_of_sound: np.ndarray
    """ndarray: Speed of sound [m/s]."""

    def summary(self) -> dict:
        """Sea-level-relative density and pressure at the highest altitude."""
        return {
            "max_altitude_m": float(np.max(self.altitude)),
            "min_density_ratio": float(np.min(self.density) / 1.225),
            "min_pressure_ratio": float(np.min(self.pressure) / 101325.0),
        }

    def plot(self, ax=None):
        """Plot temperature, pressure and density against altitude; return ``(fig, axes)``."""
        from engineeringkit.aerodynamics.visualizers.plots import plot_atmosphere

        axes = plot_atmosphere(self, axes=ax)
        return axes[0].figure, axes

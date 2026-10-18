"""Physical constants used by the radiation and phase-change models (CODATA 2018, SI)."""

from __future__ import annotations

__all__ = ["STEFAN_BOLTZMANN", "WIEN_DISPLACEMENT", "PLANCK", "SPEED_OF_LIGHT", "BOLTZMANN", "STANDARD_GRAVITY"]

STEFAN_BOLTZMANN = 5.670374419e-8
"""Stefan-Boltzmann constant :math:`\\sigma` [W/(m^2 K^4)]."""

WIEN_DISPLACEMENT = 2.897771955e-3
"""Wien wavelength-displacement constant :math:`b` [m K]."""

PLANCK = 6.62607015e-34
"""Planck constant :math:`h` [J s] (exact)."""

SPEED_OF_LIGHT = 299792458.0
"""Speed of light in vacuum :math:`c` [m/s] (exact)."""

BOLTZMANN = 1.380649e-23
"""Boltzmann constant :math:`k` [J/K] (exact)."""

STANDARD_GRAVITY = 9.80665
"""Standard acceleration of gravity :math:`g` [m/s^2]."""

"""Result containers for engineeringkit.aerodynamics.

Most aerodynamic models here are closed-form algebraic relations that return
plain floats or arrays. The models that produce structured output (the NACA
four-digit section generator, the standard atmosphere, oblique shocks, Fanno
and Rayleigh duct flow, the lifting-line solver and the panel method) return
the frozen dataclasses below so plotting helpers have a stable interface.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["Airfoil", "AtmosphereState", "ObliqueShockResult", "DuctFlowState", "LiftingLineResult", "PanelResult"]


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


@dataclass(frozen=True)
class ObliqueShockResult:
    """State change across an attached oblique shock in a perfect gas."""

    mach1: float
    """float: Upstream Mach number [-]."""

    theta: float
    """float: Flow deflection (wedge half-angle) [rad]."""

    beta: float
    """float: Shock-wave angle measured from the upstream flow [rad]."""

    mach2: float
    """float: Downstream Mach number [-]."""

    pressure_ratio: float
    """float: Static pressure ratio :math:`p_2/p_1` [-]."""

    density_ratio: float
    """float: Density ratio :math:`\\rho_2/\\rho_1` [-]."""

    temperature_ratio: float
    """float: Static temperature ratio :math:`T_2/T_1` [-]."""

    stagnation_pressure_ratio: float
    """float: Total pressure ratio :math:`p_{02}/p_{01}` [-], at most one."""

    gamma: float
    """float: Ratio of specific heats [-]."""

    strong: bool
    """bool: ``True`` for the strong-shock branch, ``False`` for the weak one."""

    def summary(self) -> dict:
        """Wave angle in degrees and the downstream state."""
        return {
            "beta_deg": float(np.degrees(self.beta)),
            "mach2": float(self.mach2),
            "pressure_ratio": float(self.pressure_ratio),
            "temperature_ratio": float(self.temperature_ratio),
            "stagnation_pressure_ratio": float(self.stagnation_pressure_ratio),
        }

    def plot(self, ax=None):
        """Mark this shock on the theta-beta-M chart for its Mach number; return ``(fig, ax)``."""
        from engineeringkit.aerodynamics.visualizers.plots import plot_theta_beta_mach

        ax = plot_theta_beta_mach([self.mach1], gamma=self.gamma, ax=ax)
        ax.plot(np.degrees(self.theta), np.degrees(self.beta), "ko")
        return ax.figure, ax


@dataclass(frozen=True)
class DuctFlowState:
    """Fanno or Rayleigh flow state relative to the sonic (starred) reference state.

    Every ratio is the local value divided by the value at :math:`M = 1` on the
    same Fanno or Rayleigh line.
    """

    kind: str
    """str: ``"fanno"`` (adiabatic flow with friction) or ``"rayleigh"`` (frictionless flow with heat addition)."""

    mach: np.ndarray
    """ndarray: Mach number [-]."""

    temperature_ratio: np.ndarray
    """ndarray: :math:`T/T^*` [-]."""

    pressure_ratio: np.ndarray
    """ndarray: :math:`p/p^*` [-]."""

    density_ratio: np.ndarray
    """ndarray: :math:`\\rho/\\rho^*` [-]; also :math:`V^*/V` by continuity."""

    stagnation_pressure_ratio: np.ndarray
    """ndarray: :math:`p_0/p_0^*` [-]."""

    stagnation_temperature_ratio: np.ndarray
    """ndarray: :math:`T_0/T_0^*` [-]; identically one for Fanno flow."""

    friction_parameter: np.ndarray
    """ndarray: :math:`f L^*/D` [-] (Darcy friction factor) to reach :math:`M = 1`; zero for Rayleigh flow."""

    def summary(self) -> dict:
        """The ratios at the first Mach number in the state."""
        first = {name: float(np.ravel(getattr(self, name))[0]) for name in self._RATIOS}
        return {"kind": self.kind, "mach": float(np.ravel(self.mach)[0]), **first}

    _RATIOS = (
        "temperature_ratio",
        "pressure_ratio",
        "density_ratio",
        "stagnation_pressure_ratio",
        "stagnation_temperature_ratio",
        "friction_parameter",
    )

    def plot(self, ax=None):
        """Plot the property ratios against Mach number; return ``(fig, ax)``."""
        from engineeringkit.aerodynamics.visualizers.plots import plot_duct_flow

        ax = plot_duct_flow(self, ax=ax)
        return ax.figure, ax


@dataclass(frozen=True)
class LiftingLineResult:
    """Glauert's Fourier solution of Prandtl's lifting-line equation for one wing."""

    span: float
    """float: Wing span [m]."""

    area: float
    """float: Planform area [m^2]."""

    aspect_ratio: float
    """float: Aspect ratio :math:`b^2/S` [-]."""

    coefficients: np.ndarray
    """ndarray: Fourier coefficients :math:`A_1 \\ldots A_N` of the circulation [-]."""

    y: np.ndarray
    """ndarray: Spanwise collocation stations, from the left tip to the right [m]."""

    chord: np.ndarray
    """ndarray: Chord at each station [m]."""

    cl: np.ndarray
    """ndarray: Section lift coefficient at each station [-]."""

    induced_angle: np.ndarray
    """ndarray: Downwash angle :math:`\\alpha_i` at each station [rad]."""

    circulation_per_speed: np.ndarray
    """ndarray: :math:`\\Gamma/V_\\infty = \\tfrac12 c\\,c_l` at each station [m]."""

    CL: float
    """float: Wing lift coefficient :math:`\\pi AR A_1` [-]."""

    CDi: float
    """float: Induced drag coefficient :math:`\\pi AR \\sum n A_n^2` [-]."""

    span_efficiency: float
    """float: Span efficiency :math:`e = C_L^2/(\\pi AR C_{D,i})` [-]; NaN at zero lift and drag."""

    def summary(self) -> dict:
        """Wing coefficients, span efficiency and the largest section lift coefficient."""
        return {
            "aspect_ratio": float(self.aspect_ratio),
            "CL": float(self.CL),
            "CDi": float(self.CDi),
            "span_efficiency": float(self.span_efficiency),
            "max_section_cl": float(np.max(self.cl)),
        }

    def plot(self, ax=None):
        """Plot the spanwise lift distribution against the elliptic one; return ``(fig, ax)``."""
        from engineeringkit.aerodynamics.visualizers.plots import plot_spanwise_loading

        ax = plot_spanwise_loading(self, ax=ax)
        return ax.figure, ax


@dataclass(frozen=True)
class PanelResult:
    """Hess-Smith panel solution for a two-dimensional body in a uniform stream.

    Lengths are in units of the input coordinates; the coefficients are
    referred to the chord, the ``x`` extent of the body.
    """

    nodes: np.ndarray
    """ndarray: Panel end points, shape ``(N + 1, 2)``, clockwise from the trailing edge."""

    x: np.ndarray
    """ndarray: Panel control points (midpoints), ``x`` coordinates."""

    y: np.ndarray
    """ndarray: Panel control points (midpoints), ``y`` coordinates."""

    tangential_velocity: np.ndarray
    """ndarray: Surface speed :math:`V_t/V_\\infty` along each panel's clockwise tangent [-]."""

    cp: np.ndarray
    """ndarray: Pressure coefficient :math:`1 - (V_t/V_\\infty)^2` at each control point [-]."""

    source_strengths: np.ndarray
    """ndarray: Source strength per unit length on each panel, divided by :math:`V_\\infty` [-]."""

    vortex_strength: float
    """float: Uniform vortex-sheet strength (clockwise positive), divided by :math:`V_\\infty` [-]."""

    alpha: float
    """float: Angle of attack [rad]."""

    chord: float
    """float: Chord, the ``x`` extent of the nodes."""

    cl: float
    """float: Lift coefficient from the circulation, :math:`2\\Gamma/(V_\\infty c)` [-]."""

    cl_pressure: float
    """float: Lift coefficient from integrating the surface pressure [-]."""

    cm_quarter_chord: float
    """float: Pitching moment about the quarter chord from the surface pressure, nose up positive [-]."""

    def summary(self) -> dict:
        """Lift and moment coefficients and the suction peak."""
        return {
            "panels": len(self.cp),
            "cl": float(self.cl),
            "cl_pressure": float(self.cl_pressure),
            "cm_quarter_chord": float(self.cm_quarter_chord),
            "min_cp": float(np.min(self.cp)),
        }

    def plot(self, ax=None):
        """Plot :math:`-C_p` along the chord; return ``(fig, ax)``."""
        from engineeringkit.aerodynamics.visualizers.plots import plot_pressure_coefficient

        ax = plot_pressure_coefficient(self, ax=ax)
        return ax.figure, ax

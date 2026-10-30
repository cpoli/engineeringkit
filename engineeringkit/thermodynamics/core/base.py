"""Result containers for engineeringkit.thermodynamics.

Most thermodynamic models return plain floats or arrays. A steam or
refrigerant state carries several properties at once, and a cycle solved
state by state carries its states and energy balance, so they get frozen
dataclasses with the usual ``summary()`` (and ``plot()`` for the cycles).
A Peng-Robinson fluid is described by a frozen ``CubicFluid``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import interval as _interval
from engineeringkit._validation import positive as _positive

__all__ = ["SteamState", "RankineResult", "BraytonResult", "CubicFluid", "RefrigerantState", "VaporCompressionResult"]


@dataclass(frozen=True)
class SteamState:
    """Output of :func:`engineeringkit.thermodynamics.steam_properties`.

    Each field is a float for a scalar call and an array of the broadcast
    input shape for an array call.
    """

    p: float
    """float or ndarray: Pressure [Pa]."""

    T: float
    """float or ndarray: Temperature [K]."""

    v: float
    """float or ndarray: Specific volume [m^3/kg]."""

    h: float
    """float or ndarray: Specific enthalpy [J/kg], zero for liquid at the triple point."""

    s: float
    """float or ndarray: Specific entropy [J/(kg K)], zero for liquid at the triple point."""

    u: float
    """float or ndarray: Specific internal energy [J/kg]."""

    cp: float
    """float or ndarray: Isobaric specific heat [J/(kg K)]; NaN inside the two-phase dome."""

    x: float
    """float or ndarray: Vapour quality [-] inside the dome; NaN for a single-phase state."""

    region: int
    """int or ndarray: IF97 region, 1 (compressed liquid), 2 (vapour), 3 (near-critical), 4 (wet mixture) or 5 (above 1073.15 K)."""

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


def _isobar_or_line(ax, a, b, isobar):
    """Draw from state ``a`` to ``b``: along ``isobar(a, b)`` -> (s, T) if they share a pressure, else straight."""
    if np.isclose(a[2], b[2], rtol=1e-12):
        s, T = isobar(a, b)
    else:
        s, T = np.array([a[0], b[0]]), np.array([a[1], b[1]])
    ax.plot(s / 1e3, T - 273.15, "C0-")


@dataclass(frozen=True)
class RankineResult:
    """Output of :func:`engineeringkit.thermodynamics.rankine_cycle`.

    Energies are per kilogram of steam through the boiler.
    """

    labels: tuple
    """tuple of str: Name of each state, in the order the water meets them."""

    states: tuple
    """tuple of SteamState: The state points, matching ``labels``."""

    efficiency: float
    """float: Thermal efficiency, net work over heat added [-]."""

    turbine_work: float
    """float: Turbine work [J/kg], counting the bled steam only up to its bleed point."""

    pump_work: float
    """float: Work of all pumps [J/kg]."""

    heat_in: float
    """float: Heat added in the boiler and reheater [J/kg]."""

    heat_out: float
    """float: Heat rejected in the condenser [J/kg]."""

    bleed_fraction: float
    """float: Fraction of the boiler flow bled to the open feedwater heater [-]; zero without one."""

    @property
    def net_work(self):
        """float: Turbine work less pump work [J/kg]."""
        return self.turbine_work - self.pump_work

    @property
    def back_work_ratio(self):
        """float: Pump work over turbine work [-]."""
        return self.pump_work / self.turbine_work

    @property
    def turbine_exit_quality(self):
        """float: Vapour quality at the condenser inlet [-]; NaN if the steam is still superheated."""
        return self.states[self.labels.index("turbine exit")].x

    def summary(self):
        """Efficiency, energies per kilogram, back-work ratio, bleed fraction and exit quality."""
        return {
            "efficiency": self.efficiency,
            "net_work": self.net_work,
            "turbine_work": self.turbine_work,
            "pump_work": self.pump_work,
            "heat_in": self.heat_in,
            "heat_out": self.heat_out,
            "back_work_ratio": self.back_work_ratio,
            "bleed_fraction": self.bleed_fraction,
            "turbine_exit_quality": self.turbine_exit_quality,
        }

    def plot(self, ax=None):
        """Draw the cycle on the T-s diagram over the water dome; return ``(fig, ax)``.

        Heat exchanges follow their isobar; pumps and turbine stages are
        drawn as straight lines between their end states.
        """
        import matplotlib.pyplot as plt

        from engineeringkit.thermodynamics.systems.steam import steam_properties
        from engineeringkit.thermodynamics.visualizers.plots import plot_saturation_dome

        if ax is None:
            _, ax = plt.subplots()
        plot_saturation_dome(ax)

        def isobar(a, b):
            path = steam_properties(a[2], h=np.linspace(a[3], b[3], 40))
            return path.s, path.T

        # the main loop: feed train, boiler, every expansion state except the bleed, condenser
        order = [i for i, label in enumerate(self.labels) if label != "bleed"] + [0]
        points = [(st.s, st.T, st.p, st.h) for st in self.states]
        for i, j in zip(order[:-1], order[1:]):
            _isobar_or_line(ax, points[i], points[j], isobar)
        ax.plot([p[0] / 1e3 for p in points], [p[1] - 273.15 for p in points], "C0o", ms=4)
        ax.set(xlabel="entropy [kJ/(kg K)]", ylabel="temperature [°C]", title=f"Rankine cycle, η = {self.efficiency:.1%}")
        return ax.figure, ax


@dataclass(frozen=True)
class BraytonResult:
    """Output of :func:`engineeringkit.thermodynamics.brayton_cycle`.

    Energies are per kilogram of air (air-standard, constant specific heats).
    """

    labels: tuple
    """tuple of str: Name of each state, in the order the air meets them."""

    T: np.ndarray
    """ndarray: Temperature of each state [K]."""

    p: np.ndarray
    """ndarray: Pressure of each state [Pa]."""

    s: np.ndarray
    """ndarray: Entropy of each state relative to the compressor inlet [J/(kg K)]."""

    efficiency: float
    """float: Thermal efficiency, net work over heat added [-]."""

    compressor_work: float
    """float: Work of all compressor stages [J/kg]."""

    turbine_work: float
    """float: Work of all turbine stages [J/kg]."""

    heat_in: float
    """float: Heat added in the combustor and reheaters [J/kg]."""

    heat_regenerated: float
    """float: Heat passed from the exhaust to the compressed air [J/kg]."""

    cp: float
    """float: Isobaric specific heat used [J/(kg K)]."""

    @property
    def net_work(self):
        """float: Turbine work less compressor work [J/kg]."""
        return self.turbine_work - self.compressor_work

    @property
    def heat_out(self):
        """float: Heat rejected in the exhaust and intercoolers [J/kg]."""
        return self.heat_in - self.net_work

    @property
    def back_work_ratio(self):
        """float: Compressor work over turbine work [-]."""
        return self.compressor_work / self.turbine_work

    def summary(self):
        """Efficiency, energies per kilogram and the back-work ratio."""
        return {
            "efficiency": self.efficiency,
            "net_work": self.net_work,
            "compressor_work": self.compressor_work,
            "turbine_work": self.turbine_work,
            "heat_in": self.heat_in,
            "heat_out": self.heat_out,
            "heat_regenerated": self.heat_regenerated,
            "back_work_ratio": self.back_work_ratio,
        }

    def plot(self, ax=None):
        """Draw the cycle on T-s axes; return ``(fig, ax)``.

        Constant-pressure processes follow :math:`s = c_p\\ln T + \\text{const}`;
        compressor and turbine stages are straight lines.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()

        def isobar(a, b):
            T = np.linspace(a[1], b[1], 40)
            return a[0] + self.cp * np.log(T / a[1]), T

        points = list(zip(self.s, self.T, self.p))
        for a, b in zip(points, points[1:] + points[:1]):
            _isobar_or_line(ax, a, b, isobar)
        ax.plot(self.s / 1e3, self.T - 273.15, "C0o", ms=4)
        ax.set(xlabel="entropy [kJ/(kg K)]", ylabel="temperature [°C]", title=f"Brayton cycle, η = {self.efficiency:.1%}")
        return ax.figure, ax


@dataclass(frozen=True)
class CubicFluid:
    """A pure fluid described by the Peng-Robinson equation of state.

    Parameters
    ----------
    name : str
        Display name.
    Tc : float
        Critical temperature [K].
    pc : float
        Critical pressure [Pa].
    omega : float
        Pitzer's acentric factor [-].
    molar_mass : float
        Molar mass [kg/mol].
    cp0 : float
        Ideal-gas isobaric specific heat [J/(kg K)], taken as constant.
    """

    name: str
    Tc: float
    pc: float
    omega: float
    molar_mass: float
    cp0: float

    def __post_init__(self):
        _positive(Tc=self.Tc, pc=self.pc, molar_mass=self.molar_mass, cp0=self.cp0)
        _interval("omega", self.omega, -0.5, 1.5)

    def summary(self):
        """The fluid constants as a dict."""
        return {"name": self.name, "Tc": self.Tc, "pc": self.pc, "omega": self.omega, "molar_mass": self.molar_mass, "cp0": self.cp0}


@dataclass(frozen=True)
class RefrigerantState:
    """Output of :func:`engineeringkit.thermodynamics.refrigerant_properties`."""

    p: float
    """float: Pressure [Pa]."""

    T: float
    """float: Temperature [K]."""

    v: float
    """float: Specific volume [m^3/kg]."""

    h: float
    """float: Specific enthalpy [J/kg], 200 kJ/kg for saturated liquid at 0 °C (IIR reference)."""

    s: float
    """float: Specific entropy [J/(kg K)], 1 kJ/(kg K) for saturated liquid at 0 °C (IIR reference)."""

    x: float
    """float: Vapour quality [-] inside the dome; NaN for a single-phase state."""

    phase: str
    """str: ``"liquid"``, ``"vapour"``, ``"mixture"`` or ``"supercritical"``."""

    def summary(self):
        """All properties as a dict."""
        return {"p": self.p, "T": self.T, "v": self.v, "h": self.h, "s": self.s, "x": self.x, "phase": self.phase}


@dataclass(frozen=True)
class VaporCompressionResult:
    """Output of :func:`engineeringkit.thermodynamics.vapor_compression_cycle`.

    Energies are per kilogram of refrigerant.
    """

    fluid: CubicFluid
    """CubicFluid: The refrigerant."""

    labels: tuple
    """tuple of str: Compressor inlet, compressor outlet, condenser outlet and evaporator inlet."""

    states: tuple
    """tuple of RefrigerantState: The four state points, matching ``labels``."""

    q_evaporator: float
    """float: Heat absorbed in the evaporator (refrigerating effect) [J/kg]."""

    q_condenser: float
    """float: Heat rejected in the condenser [J/kg]."""

    compressor_work: float
    """float: Compressor work [J/kg]."""

    @property
    def cop_refrigerator(self):
        """float: Refrigerating effect over compressor work [-]."""
        return self.q_evaporator / self.compressor_work

    @property
    def cop_heat_pump(self):
        """float: Condenser heat over compressor work [-], one more than the refrigerator COP."""
        return self.q_condenser / self.compressor_work

    def summary(self):
        """COPs, energies per kilogram and the compressor outlet temperature."""
        return {
            "cop_refrigerator": self.cop_refrigerator,
            "cop_heat_pump": self.cop_heat_pump,
            "q_evaporator": self.q_evaporator,
            "q_condenser": self.q_condenser,
            "compressor_work": self.compressor_work,
            "compressor_outlet_T": self.states[1].T,
        }

    def plot(self, ax=None):
        """Draw the cycle on the pressure-enthalpy diagram over the fluid's dome; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        from engineeringkit.thermodynamics.systems.peng_robinson import refrigerant_properties

        if ax is None:
            _, ax = plt.subplots()
        p_dome = np.geomspace(0.02 * self.fluid.pc, 0.98 * self.fluid.pc, 40)
        for x in (0.0, 1.0):
            ax.plot([refrigerant_properties(self.fluid, p, x=x).h / 1e3 for p in p_dome], p_dome / 1e3, "k-", lw=1)
        h = [st.h / 1e3 for st in self.states] + [self.states[0].h / 1e3]
        p = [st.p / 1e3 for st in self.states] + [self.states[0].p / 1e3]
        ax.plot(h, p, "C0o-", ms=4)
        ax.set(
            yscale="log",
            xlabel="enthalpy [kJ/kg]",
            ylabel="pressure [kPa]",
            title=f"{self.fluid.name} vapour-compression cycle, COP = {self.cop_refrigerator:.2f}",
        )
        return ax.figure, ax

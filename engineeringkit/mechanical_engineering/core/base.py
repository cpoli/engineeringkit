"""Result containers for engineeringkit.mechanical_engineering.

Most mechanical models are closed-form relations that return plain floats or
arrays. Models whose output is a response curve or a set of modes, the
simulated oscillator, the vibration absorber, the Jeffcott rotor and Holzer's
torsional modes, return the frozen dataclasses below so they can be
summarised and plotted consistently.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import nonnegative as _nonnegative

__all__ = ["SDOFResult", "TunedMassDamper", "JeffcottResult", "HolzerResult"]


@dataclass(frozen=True)
class SDOFResult:
    """Time history of a single-degree-of-freedom oscillator."""

    t: np.ndarray
    """ndarray: Time samples [s]."""

    displacement: np.ndarray
    """ndarray: Displacement at each sample [m]."""

    velocity: np.ndarray
    """ndarray: Velocity at each sample [m/s]."""

    def summary(self) -> dict:
        """Peak absolute displacement [m] and final velocity [m/s]."""
        return {"peak_displacement_m": float(np.max(np.abs(self.displacement))), "final_velocity_m_s": float(self.velocity[-1])}

    @property
    def history(self) -> np.ndarray:
        """ndarray: Displacement and velocity as columns, shape ``(len(t), 2)``."""
        return np.column_stack((self.displacement, self.velocity))

    def plot(self):
        """Plot displacement against time; return ``(fig, ax)``."""
        from engineeringkit.visualizers import plot_history

        return plot_history(self.t, self.displacement, labels=["displacement [m]"])

    def animate(self, interval=30):
        """Animate the mass; see :func:`engineeringkit.visualizers.animate_displacement`.

        Parameters
        ----------
        interval : float, optional
            Delay between frames [ms], positive.

        Returns
        -------
        tuple
            ``(fig, animation)``.
        """
        from engineeringkit.visualizers import animate_displacement

        return animate_displacement(self.t, self.displacement, interval=interval)


@dataclass(frozen=True)
class TunedMassDamper:
    """A primary mass-spring system carrying a damped vibration absorber.

    Built by :func:`engineeringkit.mechanical_engineering.tuned_mass_damper`.
    The primary mass :math:`m_1` on stiffness :math:`k_1` is driven by a
    harmonic force; the absorber :math:`m_2` is attached to it through
    :math:`k_2` and a dashpot :math:`c_2`.
    """

    primary_mass: float
    """float: :math:`m_1` [kg]."""

    primary_stiffness: float
    """float: :math:`k_1` [N/m]."""

    absorber_mass: float
    """float: :math:`m_2` [kg]."""

    absorber_stiffness: float
    """float: :math:`k_2` [N/m]."""

    absorber_damping: float
    """float: :math:`c_2` [N s/m]."""

    mass_ratio: float
    """float: :math:`\\mu = m_2/m_1` [-]."""

    frequency_ratio: float
    """float: Absorber natural frequency over primary natural frequency [-]."""

    damping_ratio: float
    """float: :math:`c_2 / (2 m_2 \\omega_1)` [-], Den Hartog's convention."""

    @property
    def primary_frequency(self) -> float:
        """float: :math:`\\omega_1 = \\sqrt{k_1/m_1}` [rad/s]."""
        return float(np.sqrt(self.primary_stiffness / self.primary_mass))

    def response(self, r):
        """Complex steady-state amplitudes of both masses.

        Parameters
        ----------
        r : float or array_like
            Forcing frequency over :math:`\\omega_1` [-], nonnegative.

        Returns
        -------
        tuple of ndarray
            ``(X1, X2)``, each divided by the static deflection
            :math:`F_0/k_1` [-].
        """
        _nonnegative(r=r)
        w = np.asarray(r, dtype=float) * self.primary_frequency
        m1, k1, m2, k2, c2 = self.primary_mass, self.primary_stiffness, self.absorber_mass, self.absorber_stiffness, self.absorber_damping
        coupling = k2 + 1j * w * c2
        a11 = k1 + coupling - w * w * m1
        a22 = coupling - w * w * m2
        det = a11 * a22 - coupling * coupling
        return k1 * a22 / det, k1 * coupling / det

    def amplification(self, r):
        """Primary amplitude over static deflection, :math:`|X_1| k_1 / F_0`.

        Parameters
        ----------
        r : float or array_like
            Forcing frequency over :math:`\\omega_1` [-], nonnegative.

        Returns
        -------
        float or ndarray
            Amplification [-].
        """
        return np.abs(self.response(r)[0])

    def summary(self) -> dict:
        """Absorber design values and Den Hartog's fixed-point peak :math:`\\sqrt{1 + 2/\\mu}`."""
        return {
            "mass_ratio": self.mass_ratio,
            "frequency_ratio": self.frequency_ratio,
            "damping_ratio": self.damping_ratio,
            "absorber_stiffness_N_m": self.absorber_stiffness,
            "absorber_damping_N_s_m": self.absorber_damping,
            "fixed_point_amplification": float(np.sqrt(1 + 2 / self.mass_ratio)),
        }

    def plot(self, r=None, ax=None):
        """Plot amplification with and without the absorber; return ``(fig, ax)``.

        Parameters
        ----------
        r : array_like, optional
            Forcing frequency ratios [-]; 0.5 to 1.5 if omitted.
        ax : matplotlib.axes.Axes, optional
            Axes to draw on; a new figure is created if omitted.

        Returns
        -------
        tuple
            ``(fig, ax)``.
        """
        from engineeringkit.mechanical_engineering.visualizers.plots import plot_frequency_response

        r = np.linspace(0.5, 1.5, 1001) if r is None else np.asarray(r, dtype=float)
        with np.errstate(divide="ignore"):
            bare = 1 / np.abs(1 - r * r)
        ax = plot_frequency_response(r, np.where(np.isfinite(bare), bare, np.nan), ax=ax, label="no absorber")
        plot_frequency_response(r, self.amplification(r), ax=ax, label=f"absorber, mu = {self.mass_ratio:g}")
        return ax.figure, ax


@dataclass(frozen=True)
class JeffcottResult:
    """Steady synchronous whirl of an unbalanced disk on a flexible shaft."""

    omega: np.ndarray
    """ndarray: Spin speeds [rad/s]."""

    amplitude: np.ndarray
    """ndarray: Whirl radius of the shaft centre [m]."""

    phase: np.ndarray
    """ndarray: Angle by which the whirl lags the unbalance [rad], in [0, pi]."""

    critical_speed: float
    """float: :math:`\\sqrt{k/m}` [rad/s]."""

    eccentricity: float
    """float: Distance from shaft centre to disk mass centre [m]."""

    def summary(self) -> dict:
        """Critical speed [rad/s] and the largest whirl radius over ``eccentricity`` [-]."""
        return {"critical_speed_rad_s": self.critical_speed, "peak_amplitude_ratio": float(np.max(self.amplitude) / self.eccentricity)}

    def plot(self, ax=None):
        """Plot whirl radius over eccentricity against speed ratio; return ``(fig, ax)``."""
        from engineeringkit.mechanical_engineering.visualizers.plots import plot_frequency_response

        ax = plot_frequency_response(self.omega / self.critical_speed, self.amplitude / self.eccentricity, ax=ax)
        ax.set(xlabel="speed / critical speed", ylabel="whirl radius / eccentricity")
        return ax.figure, ax


@dataclass(frozen=True)
class HolzerResult:
    """Natural frequencies and mode shapes of a free-free torsional chain."""

    frequencies: np.ndarray
    """ndarray: Nonzero natural frequencies in increasing order [rad/s]."""

    mode_shapes: np.ndarray
    """ndarray: Disk twist amplitudes, one row per mode, largest entry 1 [-]."""

    def summary(self) -> dict:
        """Number of modes found and their frequencies [Hz]."""
        return {"modes": int(self.frequencies.size), "frequencies_hz": (self.frequencies / (2 * np.pi)).tolist()}

    def plot(self, ax=None):
        """Plot every mode shape against disk number; return ``(fig, ax)``."""
        from engineeringkit.mechanical_engineering.visualizers.plots import plot_mode_shapes

        ax = plot_mode_shapes(self.mode_shapes, self.frequencies, ax=ax)
        return ax.figure, ax

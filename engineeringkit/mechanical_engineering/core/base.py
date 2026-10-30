"""Result containers for engineeringkit.mechanical_engineering.

Most mechanical models are closed-form relations that return plain floats or
arrays. Models whose output is a response curve, a time history, a set of
modes or a linked set of quantities, the simulated oscillators, the vibration
absorber, the Jeffcott rotor, normal modes, the bolted joint, the balanced
rotor, linkage and cam motions, return the frozen dataclasses below so they
can be summarised and plotted consistently.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive

__all__ = [
    "SDOFResult",
    "TunedMassDamper",
    "JeffcottResult",
    "HolzerResult",
    "ModalResult",
    "MDOFResult",
    "BoltedJoint",
    "BalanceResult",
    "FourBarMotion",
    "SliderCrankMotion",
    "CamMotion",
]


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


@dataclass(frozen=True)
class ModalResult:
    """Natural frequencies and mass-normalised mode shapes of an undamped M-K system."""

    frequencies: np.ndarray
    """ndarray: Natural frequencies in increasing order [rad/s]; zero for rigid-body modes."""

    mode_shapes: np.ndarray
    """ndarray: One row per mode, scaled so :math:`\\phi^T M \\phi = 1` and the largest entry is positive."""

    def summary(self) -> dict:
        """Number of modes and their frequencies [Hz]."""
        return {"modes": int(self.frequencies.size), "frequencies_hz": (self.frequencies / (2 * np.pi)).tolist()}

    def plot(self, ax=None):
        """Plot every mode shape against degree of freedom; return ``(fig, ax)``."""
        from engineeringkit.mechanical_engineering.visualizers.plots import plot_mode_shapes

        ax = plot_mode_shapes(self.mode_shapes, self.frequencies, ax=ax)
        ax.set(xlabel="degree of freedom", ylabel="mass-normalised amplitude")
        return ax.figure, ax


@dataclass(frozen=True)
class MDOFResult:
    """Time history of a system with several (generalised) coordinates."""

    t: np.ndarray
    """ndarray: Time samples [s]."""

    displacement: np.ndarray
    """ndarray: Coordinates at each sample, shape ``(len(t), n)`` [m or rad]."""

    velocity: np.ndarray
    """ndarray: Coordinate rates at each sample, shape ``(len(t), n)`` [m/s or rad/s]."""

    def summary(self) -> dict:
        """Peak absolute value of each coordinate."""
        return {"coordinates": int(self.displacement.shape[1]), "peak_displacement": np.max(np.abs(self.displacement), axis=0).tolist()}

    def plot(self, ax=None):
        """Plot every coordinate against time; return ``(fig, ax)``."""
        from engineeringkit.visualizers import plot_history

        return plot_history(self.t, self.displacement, labels=[f"q{i + 1}" for i in range(self.displacement.shape[1])], ax=ax)


@dataclass(frozen=True)
class BoltedJoint:
    """A preloaded bolt and the members it clamps, carrying an external separating load.

    Built by :func:`engineeringkit.mechanical_engineering.bolted_joint`. Bolt
    and members act as springs in parallel under the external load
    :math:`P`: the bolt takes the share :math:`CP` of it, with joint constant
    :math:`C = k_b/(k_b + k_m)`, and the members lose the rest of their
    clamping force until they separate.
    """

    bolt_stiffness: float
    """float: :math:`k_b` [N/m]."""

    member_stiffness: float
    """float: :math:`k_m` [N/m]."""

    preload: float
    """float: Initial clamping force :math:`F_i` [N]."""

    load: float
    """float: External separating load :math:`P` [N]."""

    @property
    def joint_constant(self) -> float:
        """float: :math:`C = k_b/(k_b + k_m)` [-]."""
        return self.bolt_stiffness / (self.bolt_stiffness + self.member_stiffness)

    @property
    def separation_load(self) -> float:
        """float: External load :math:`P_0 = F_i/(1 - C)` at which the members part [N]."""
        return self.preload / (1 - self.joint_constant)

    @property
    def bolt_load(self) -> float:
        """float: Bolt tension :math:`F_i + CP` [N]; once the members separate it carries the whole load :math:`P`."""
        if self.load >= self.separation_load:
            return self.load
        return self.preload + self.joint_constant * self.load

    @property
    def member_load(self) -> float:
        """float: Member compression :math:`F_i - (1 - C)P` [N], zero once separated."""
        return max(0.0, self.preload - (1 - self.joint_constant) * self.load)

    def load_factor(self, proof_load) -> float:
        """Factor on the external load before the bolt reaches its proof load.

        Parameters
        ----------
        proof_load : float
            :math:`F_p = S_p A_t` [N], positive and above the preload.

        Returns
        -------
        float
            :math:`n_L = (F_p - F_i)/(CP)` [-].
        """
        _positive(proof_load=proof_load)
        if proof_load <= self.preload:
            raise ValueError("proof_load must exceed the preload")
        return float((proof_load - self.preload) / (self.joint_constant * self.load))

    def summary(self) -> dict:
        """Joint constant, bolt and member loads [N] and the separation factor :math:`P_0/P` [-]."""
        return {
            "joint_constant": self.joint_constant,
            "bolt_load_N": self.bolt_load,
            "member_load_N": self.member_load,
            "separation_load_N": self.separation_load,
            "separation_factor": self.separation_load / self.load,
        }

    def plot(self, ax=None):
        """Draw Rötscher's joint diagram, force against elongation; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        kb, km, fi = self.bolt_stiffness, self.member_stiffness, self.preload
        db, dm = fi / kb, fi / km
        ax.plot([0, db + dm], [0, (db + dm) * kb], label="bolt")
        ax.plot([db, db + dm], [fi, 0], label="members")
        delta = self.load / (kb + km) if self.load < self.separation_load else self.load / kb - db
        ax.plot([db + delta] * 2, [self.member_load, self.bolt_load], "k-", lw=2, label=f"external load P = {self.load:.3g} N")
        ax.axhline(fi, color="gray", lw=0.5)
        ax.set(xlabel="elongation [m]", ylabel="force [N]", xlim=(0, None), ylim=(0, None))
        ax.legend()
        return ax.figure, ax


@dataclass(frozen=True)
class BalanceResult:
    """Correction masses that balance a set of rotating masses on a shaft."""

    mass_radius: np.ndarray
    """ndarray: Correction mass times its radius in each correction plane [kg m]."""

    angle: np.ndarray
    """ndarray: Angular position of each correction [rad], in (-pi, pi]."""

    plane_positions: np.ndarray
    """ndarray: Axial positions of the correction planes [m]; empty for static (one-plane) balance."""

    def summary(self) -> dict:
        """Correction mass-radius products [kg m] and angles [deg]."""
        return {"mass_radius_kg_m": self.mass_radius.tolist(), "angle_deg": np.degrees(self.angle).tolist()}

    def plot(self, ax=None):
        """Draw each correction as a vector in the plane of rotation; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(subplot_kw={"projection": "polar"})
        for i, (mr, angle) in enumerate(zip(self.mass_radius, self.angle)):
            ax.annotate("", xy=(angle, mr), xytext=(0, 0), arrowprops={"arrowstyle": "->", "color": f"C{i}", "lw": 2})
            ax.plot([angle], [mr], "o", color=f"C{i}", label=f"plane {i + 1}: {mr:.3g} kg m")
        ax.set_rmax(1.15 * max(float(np.max(self.mass_radius)), 1e-300))
        ax.legend(loc="lower left", fontsize=8)
        return ax.figure, ax


@dataclass(frozen=True)
class FourBarMotion:
    """Positions, angular velocities and angular accelerations of a four-bar linkage."""

    theta2: np.ndarray
    """ndarray: Input crank angle [rad]."""

    theta3: np.ndarray
    """ndarray: Coupler angle [rad]."""

    theta4: np.ndarray
    """ndarray: Output link angle [rad]."""

    omega3: np.ndarray
    """ndarray: Coupler angular velocity [rad/s]."""

    omega4: np.ndarray
    """ndarray: Output angular velocity [rad/s]."""

    alpha3: np.ndarray
    """ndarray: Coupler angular acceleration [rad/s^2]."""

    alpha4: np.ndarray
    """ndarray: Output angular acceleration [rad/s^2]."""

    def summary(self) -> dict:
        """Peak absolute output angular velocity [rad/s] and acceleration [rad/s^2]."""
        return {"peak_omega4_rad_s": float(np.max(np.abs(self.omega4))), "peak_alpha4_rad_s2": float(np.max(np.abs(self.alpha4)))}

    def plot(self, ax=None):
        """Plot output angular velocity and acceleration against crank angle; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        degrees = np.degrees(np.atleast_1d(self.theta2))
        ax.plot(degrees, np.atleast_1d(self.omega4), label="omega4 [rad/s]")
        ax.plot(degrees, np.atleast_1d(self.alpha4), label="alpha4 [rad/s^2]")
        ax.set(xlabel="crank angle [deg]")
        ax.legend()
        return ax.figure, ax


@dataclass(frozen=True)
class SliderCrankMotion:
    """Slider and connecting-rod motion of an in-line slider-crank."""

    theta: np.ndarray
    """ndarray: Crank angle from top dead centre [rad]."""

    position: np.ndarray
    """ndarray: Slider position from the crank axis [m]."""

    velocity: np.ndarray
    """ndarray: Slider velocity [m/s]."""

    acceleration: np.ndarray
    """ndarray: Slider acceleration [m/s^2]."""

    rod_angle: np.ndarray
    """ndarray: Connecting-rod angle to the line of stroke [rad]."""

    rod_angular_velocity: np.ndarray
    """ndarray: Rate of change of ``rod_angle`` [rad/s]."""

    rod_angular_acceleration: np.ndarray
    """ndarray: Second rate of change of ``rod_angle`` [rad/s^2]."""

    def summary(self) -> dict:
        """Peak absolute slider velocity [m/s] and acceleration [m/s^2]."""
        return {
            "peak_velocity_m_s": float(np.max(np.abs(self.velocity))),
            "peak_acceleration_m_s2": float(np.max(np.abs(self.acceleration))),
        }

    def plot(self, ax=None):
        """Plot slider velocity and acceleration against crank angle on twin axes; return ``(fig, ax)``."""
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        degrees = np.degrees(np.atleast_1d(self.theta))
        ax.plot(degrees, np.atleast_1d(self.velocity), label="velocity")
        ax.set(xlabel="crank angle [deg]", ylabel="slider velocity [m/s]")
        twin = ax.twinx()
        twin.plot(degrees, np.atleast_1d(self.acceleration), "C1", label="acceleration")
        twin.set(ylabel="slider acceleration [m/s^2]")
        return ax.figure, ax


@dataclass(frozen=True)
class CamMotion:
    """Follower displacement and its derivatives with respect to cam angle over one rise."""

    theta: np.ndarray
    """ndarray: Cam angle from the start of the rise [rad]."""

    displacement: np.ndarray
    """ndarray: Follower lift [m]."""

    velocity: np.ndarray
    """ndarray: :math:`ds/d\\theta` [m/rad]; multiply by :math:`\\omega` for m/s."""

    acceleration: np.ndarray
    """ndarray: :math:`d^2s/d\\theta^2` [m/rad^2]; multiply by :math:`\\omega^2` for m/s^2."""

    jerk: np.ndarray
    """ndarray: :math:`d^3s/d\\theta^3` [m/rad^3]; multiply by :math:`\\omega^3` for m/s^3."""

    law: str
    """str: Name of the displacement law."""

    def summary(self) -> dict:
        """Law, lift [m] and peak absolute velocity, acceleration and jerk per radian of cam."""
        return {
            "law": self.law,
            "lift_m": float(np.max(self.displacement)),
            "peak_velocity_m_rad": float(np.max(np.abs(self.velocity))),
            "peak_acceleration_m_rad2": float(np.max(np.abs(self.acceleration))),
            "peak_jerk_m_rad3": float(np.max(np.abs(self.jerk))),
        }

    def plot(self, axs=None):
        """Draw the s-v-a-j diagram, one panel per derivative; return ``(fig, axs)``."""
        import matplotlib.pyplot as plt

        if axs is None:
            _, axs = plt.subplots(4, 1, sharex=True, figsize=(6, 7))
        degrees = np.degrees(self.theta)
        for ax, values, label in zip(
            axs, (self.displacement, self.velocity, self.acceleration, self.jerk), ("s [m]", "v [m/rad]", "a [m/rad^2]", "j [m/rad^3]")
        ):
            ax.plot(degrees, values)
            ax.set(ylabel=label)
            ax.grid(True, alpha=0.3)
        axs[0].set(title=f"{self.law} rise")
        axs[-1].set(xlabel="cam angle [deg]")
        return axs[0].figure, axs

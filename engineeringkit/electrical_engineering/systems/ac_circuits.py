r"""Sinusoidal steady state: complex impedance, resonance and the transformer.

Steinmetz replaced the differential equations of AC circuits by algebra on
complex numbers: a resistor, an inductor and a capacitor become the
impedances :math:`R`, :math:`j\omega L` and :math:`1/j\omega C`. Their
reactances cancel at the series resonance Lodge used to tune his receivers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

from engineeringkit._validation import finite as _finite
from engineeringkit._validation import nonnegative as _nonnegative
from engineeringkit._validation import positive as _positive
from engineeringkit.electrical_engineering.core.base import CircuitResult

__all__ = [
    "capacitor_impedance",
    "inductor_impedance",
    "series_rlc_impedance",
    "resonance_frequency",
    "rms_sine",
    "transformer_secondary",
    "SeriesRLC",
]


def capacitor_impedance(C, omega):
    r"""Capacitor impedance :math:`Z_C = 1/(j\omega C)`.

    Parameters
    ----------
    C : float or array_like
        Capacitance [F], positive.
    omega : float or array_like
        Angular frequency [rad/s], positive (an ideal capacitor is an open
        circuit at DC).

    Returns
    -------
    complex or ndarray of complex
        Impedance [ohm], purely negative imaginary.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> complex(capacitor_impedance(1e-3, 1e3))
    -1j
    """
    _finite(C=C, omega=omega)
    _positive(C=C, omega=omega)
    return 1 / (1j * np.asarray(omega) * C)


def inductor_impedance(L, omega):
    r"""Inductor impedance :math:`Z_L = j\omega L`.

    Parameters
    ----------
    L : float or array_like
        Inductance [H], positive.
    omega : float or array_like
        Angular frequency [rad/s], nonnegative.

    Returns
    -------
    complex or ndarray of complex
        Impedance [ohm], purely positive imaginary.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> complex(inductor_impedance(1e-3, 1e3))
    1j
    """
    _finite(L=L, omega=omega)
    _positive(L=L)
    _nonnegative(omega=omega)
    return 1j * np.asarray(omega) * L


def series_rlc_impedance(R, L, C, omega):
    r"""Series RLC impedance :math:`Z = R + j(\omega L - 1/\omega C)`.

    Parameters
    ----------
    R : float
        Resistance [ohm], nonnegative.
    L : float
        Inductance [H], positive.
    C : float
        Capacitance [F], positive.
    omega : float or array_like
        Angular frequency [rad/s], positive.

    Returns
    -------
    complex or ndarray of complex
        Impedance [ohm]; purely resistive at resonance.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    C. P. Steinmetz, "Complex quantities and their use in electrical
    engineering," *Proc. International Electrical Congress, Chicago 1893*
    (New York: AIEE, 1894) 33-74.

    Examples
    --------
    >>> complex(series_rlc_impedance(10.0, 1e-3, 1e-3, 1e3))
    (10+0j)
    """
    _finite(R=R, L=L, C=C, omega=omega)
    _nonnegative(R=R)
    return R + inductor_impedance(L, omega) + capacitor_impedance(C, omega)


def resonance_frequency(L, C):
    r"""Resonant frequency :math:`f_0 = 1/(2\pi\sqrt{LC})` of an LC circuit.

    Parameters
    ----------
    L : float or array_like
        Inductance [H], positive.
    C : float or array_like
        Capacitance [F], positive.

    Returns
    -------
    float or ndarray
        Resonant frequency [Hz].

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    O. J. Lodge, "Improvements in syntonized telegraphy without line
    wires," British Patent 11,575 (1897).

    Examples
    --------
    >>> round(float(resonance_frequency(1e-3, 1e-9)), 1)
    159154.9
    """
    _finite(L=L, C=C)
    _positive(L=L, C=C)
    return 1 / (2 * np.pi * np.sqrt(L * C))


def rms_sine(amplitude):
    r"""RMS value :math:`A/\sqrt2` of a sinusoid of peak amplitude :math:`A`.

    Parameters
    ----------
    amplitude : float or array_like
        Peak amplitude [V or A], nonnegative.

    Returns
    -------
    float or ndarray
        RMS value, same units as ``amplitude``.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    Examples
    --------
    >>> round(float(rms_sine(325.0)), 1)
    229.8
    """
    _finite(amplitude=amplitude)
    _nonnegative(amplitude=amplitude)
    return amplitude / np.sqrt(2)


def transformer_secondary(v_primary, n_primary, n_secondary):
    r"""Ideal-transformer secondary voltage :math:`V_2 = V_1 N_2/N_1`.

    Parameters
    ----------
    v_primary : float or array_like
        Primary voltage [V].
    n_primary : float
        Primary turns [-], positive.
    n_secondary : float
        Secondary turns [-], positive.

    Returns
    -------
    float or ndarray
        Secondary voltage [V]. Currents transform by the inverse ratio and
        impedances by its square.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.

    References
    ----------
    W. Stanley, "Induction-coil," U.S. Patent 349,611 (1886).

    Examples
    --------
    >>> transformer_secondary(240, 100, 50)
    120.0
    """
    _finite(v_primary=v_primary, n_primary=n_primary, n_secondary=n_secondary)
    _positive(n_primary=n_primary, n_secondary=n_secondary)
    return v_primary * n_secondary / n_primary


@dataclass(frozen=True)
class SeriesRLC:
    """Series RLC circuit driven by a voltage source, with zero initial energy.

    Parameters
    ----------
    resistance : float
        Series resistance [ohm], nonnegative (zero gives an undamped ring).
    inductance : float
        Series inductance [H], positive.
    capacitance : float
        Series capacitance [F], positive.

    Notes
    -----
    Inputs are raw SI numbers with no dimensional metadata; see
    :doc:`/guides/electrical_engineering`.
    """

    resistance: float
    """float: Series resistance [ohm]."""

    inductance: float
    """float: Series inductance [H]."""

    capacitance: float
    """float: Series capacitance [F]."""

    def __post_init__(self):
        _nonnegative(resistance=self.resistance)
        _positive(inductance=self.inductance, capacitance=self.capacitance)

    @property
    def resonance_frequency(self) -> float:
        """Undamped resonant frequency :math:`1/(2\\pi\\sqrt{LC})` [Hz]."""
        return float(resonance_frequency(self.inductance, self.capacitance))

    @property
    def quality_factor(self) -> float:
        """Quality factor :math:`Q = \\sqrt{L/C}/R` [-]; infinite for :math:`R = 0`."""
        if self.resistance == 0:
            return float("inf")
        return float(np.sqrt(self.inductance / self.capacitance) / self.resistance)

    def impedance(self, omega):
        """Complex impedance of the series circuit; see :func:`series_rlc_impedance`.

        Parameters
        ----------
        omega : float or array_like
            Angular frequency [rad/s], positive.

        Returns
        -------
        complex or ndarray
            Impedance [ohm].
        """
        return series_rlc_impedance(self.resistance, self.inductance, self.capacitance, omega)

    def transfer_function(self):
        """Admittance :math:`I(s)/V(s) = s/(Ls^2 + Rs + 1/C)` [A/V] as a SciPy ``TransferFunction``."""
        return signal.TransferFunction([1, 0], [self.inductance, self.resistance, 1 / self.capacitance])

    def step_response(self, T=None):
        """Loop current [A] after a 1 V step, from rest.

        Parameters
        ----------
        T : array_like, optional
            Strictly increasing, evenly spaced times [s]; SciPy chooses them
            if omitted.

        Returns
        -------
        CircuitResult
        """
        if T is not None:
            _nonnegative(T=T)
            if np.asarray(T).ndim != 1 or len(T) < 2 or np.any(np.diff(T) <= 0):
                raise ValueError("T must be a strictly increasing time vector")
        t, i = signal.step(self.transfer_function(), T=T)
        return CircuitResult(t, i)

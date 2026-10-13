"""engineeringkit.communications: modulation, noise, capacity, detection, coding and link budgets.

Hertz wavelengths and Friis free-space propagation with decibel link
budgets; Johnson-Nyquist thermal noise and Friis's cascaded noise factor;
Nyquist signalling, Hartley's line rate, Shannon entropy and capacity; AM,
suppressed-carrier and FM modulation with Carson's bandwidth rule; BPSK/QPSK
and their AWGN error rates; Hamming(7,4) error correction; and DFT-based
OFDM. Models are hand-written closed forms; NumPy's FFT and SciPy's erfc and
trapezoidal integration are used where they are the standard tool.
"""

from engineeringkit import __version__
from engineeringkit.communications.core.base import HammingDecodeResult
from engineeringkit.communications.systems.capacity import (
    binary_symmetric_capacity,
    entropy_binary,
    hartley_capacity,
    nyquist_symbol_rate,
    shannon_capacity,
    spectral_efficiency,
    symbol_rate,
)
from engineeringkit.communications.systems.coding import hamming74_decode, hamming74_encode, hamming_distance
from engineeringkit.communications.systems.detection import ber_bpsk_awgn, ber_qpsk_awgn, q_function
from engineeringkit.communications.systems.modulation import am_signal, bpsk, carson_bandwidth, dsb_sc, fm_signal, qpsk
from engineeringkit.communications.systems.noise import BOLTZMANN, cascaded_noise_factor, noise_figure, thermal_noise_power
from engineeringkit.communications.systems.ofdm import ofdm_demodulate, ofdm_modulate
from engineeringkit.communications.systems.propagation import (
    SPEED_OF_LIGHT,
    free_space_path_loss,
    free_space_path_loss_db,
    friis_received_power,
    link_budget_db,
    wavelength,
)
from engineeringkit.communications.utils.decibels import db, dbm, from_db, from_dbm

__all__ = [
    "__version__",
    "HammingDecodeResult",
    "db",
    "from_db",
    "dbm",
    "from_dbm",
    "SPEED_OF_LIGHT",
    "wavelength",
    "free_space_path_loss",
    "free_space_path_loss_db",
    "friis_received_power",
    "link_budget_db",
    "BOLTZMANN",
    "thermal_noise_power",
    "noise_figure",
    "cascaded_noise_factor",
    "nyquist_symbol_rate",
    "symbol_rate",
    "hartley_capacity",
    "spectral_efficiency",
    "shannon_capacity",
    "entropy_binary",
    "binary_symmetric_capacity",
    "am_signal",
    "dsb_sc",
    "fm_signal",
    "carson_bandwidth",
    "bpsk",
    "qpsk",
    "q_function",
    "ber_bpsk_awgn",
    "ber_qpsk_awgn",
    "hamming_distance",
    "hamming74_encode",
    "hamming74_decode",
    "ofdm_modulate",
    "ofdm_demodulate",
]

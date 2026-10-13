"""Concrete communication models: propagation, noise, capacity, modulation, detection, coding, OFDM."""

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

__all__ = [
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

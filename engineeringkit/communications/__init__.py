"""engineeringkit.communications: modulation, noise, capacity, detection, coding, fading and link budgets.

Hertz wavelengths and Friis free-space propagation with decibel link
budgets; antenna gain, effective aperture and the half-wave dipole;
Johnson-Nyquist thermal noise, Friis's cascaded noise factor and the output
SNR of AM and FM; Nyquist signalling, Hartley's line rate, Shannon entropy
and capacity, and Huffman coding; AM, suppressed-carrier and FM modulation
with Carson's bandwidth rule; raised-cosine pulse shaping; the matched
filter, BPSK/QPSK, M-PSK and M-QAM and their AWGN error rates; PCM,
mu-law companding and delta modulation; Hamming(7,4), CRC and
convolutional codes with Viterbi decoding; Rayleigh fading and diversity
combining; and DFT-based OFDM. Models are hand-written closed forms; NumPy's
FFT and SciPy's erfc, sine/cosine integrals, quadrature and trapezoidal
integration are used where they are the standard tool.
"""

from engineeringkit import __version__
from engineeringkit.communications.core.base import HammingDecodeResult, HuffmanCode, ViterbiDecodeResult
from engineeringkit.communications.systems.antennas import (
    FREE_SPACE_IMPEDANCE,
    aperture_gain,
    directivity,
    effective_aperture,
    half_wave_dipole_directivity,
    half_wave_dipole_pattern,
    half_wave_dipole_resistance,
)
from engineeringkit.communications.systems.capacity import (
    binary_symmetric_capacity,
    entropy_binary,
    hartley_capacity,
    nyquist_symbol_rate,
    shannon_capacity,
    spectral_efficiency,
    symbol_rate,
)
from engineeringkit.communications.systems.coding import (
    CRC16_CCITT,
    convolutional_encode,
    convolutional_free_distance,
    crc_check,
    crc_encode,
    crc_remainder,
    hamming74_decode,
    hamming74_encode,
    hamming_distance,
    viterbi_decode,
)
from engineeringkit.communications.systems.detection import (
    ber_bpsk_awgn,
    ber_qpsk_awgn,
    matched_filter,
    matched_filter_snr,
    minimum_distance_detect,
    q_function,
    ser_mpsk_awgn,
    ser_mqam_awgn,
)
from engineeringkit.communications.systems.fading import (
    ber_bpsk_mrc_rayleigh,
    ber_bpsk_rayleigh,
    maximal_ratio_combine,
    rayleigh_pdf,
    selection_outage_probability,
)
from engineeringkit.communications.systems.modulation import (
    am_signal,
    bpsk,
    carson_bandwidth,
    dsb_sc,
    fm_signal,
    psk_constellation,
    qam_constellation,
    qpsk,
)
from engineeringkit.communications.systems.noise import (
    BOLTZMANN,
    am_output_snr,
    baseband_snr,
    cascaded_noise_factor,
    fm_output_snr,
    fm_threshold_snr,
    noise_figure,
    thermal_noise_power,
)
from engineeringkit.communications.systems.ofdm import ofdm_demodulate, ofdm_modulate
from engineeringkit.communications.systems.pcm import (
    delta_demodulate,
    delta_modulate,
    delta_slope_limit,
    mu_law_compress,
    mu_law_expand,
    pcm_decode,
    pcm_encode,
    pcm_sqnr_db,
)
from engineeringkit.communications.systems.propagation import (
    SPEED_OF_LIGHT,
    free_space_path_loss,
    free_space_path_loss_db,
    friis_received_power,
    link_budget_db,
    wavelength,
)
from engineeringkit.communications.systems.pulse_shaping import (
    raised_cosine_bandwidth,
    raised_cosine_pulse,
    raised_cosine_spectrum,
    root_raised_cosine_pulse,
)
from engineeringkit.communications.systems.source_coding import entropy, huffman_code
from engineeringkit.communications.utils.decibels import db, dbm, from_db, from_dbm

__all__ = [
    "__version__",
    "HammingDecodeResult",
    "HuffmanCode",
    "ViterbiDecodeResult",
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
    "FREE_SPACE_IMPEDANCE",
    "effective_aperture",
    "aperture_gain",
    "directivity",
    "half_wave_dipole_pattern",
    "half_wave_dipole_directivity",
    "half_wave_dipole_resistance",
    "BOLTZMANN",
    "thermal_noise_power",
    "noise_figure",
    "cascaded_noise_factor",
    "baseband_snr",
    "am_output_snr",
    "fm_output_snr",
    "fm_threshold_snr",
    "nyquist_symbol_rate",
    "symbol_rate",
    "hartley_capacity",
    "spectral_efficiency",
    "shannon_capacity",
    "entropy_binary",
    "binary_symmetric_capacity",
    "entropy",
    "huffman_code",
    "am_signal",
    "dsb_sc",
    "fm_signal",
    "carson_bandwidth",
    "bpsk",
    "qpsk",
    "psk_constellation",
    "qam_constellation",
    "raised_cosine_pulse",
    "raised_cosine_spectrum",
    "root_raised_cosine_pulse",
    "raised_cosine_bandwidth",
    "q_function",
    "ber_bpsk_awgn",
    "ber_qpsk_awgn",
    "matched_filter",
    "matched_filter_snr",
    "minimum_distance_detect",
    "ser_mpsk_awgn",
    "ser_mqam_awgn",
    "pcm_encode",
    "pcm_decode",
    "pcm_sqnr_db",
    "mu_law_compress",
    "mu_law_expand",
    "delta_modulate",
    "delta_demodulate",
    "delta_slope_limit",
    "hamming_distance",
    "hamming74_encode",
    "hamming74_decode",
    "CRC16_CCITT",
    "crc_remainder",
    "crc_encode",
    "crc_check",
    "convolutional_encode",
    "convolutional_free_distance",
    "viterbi_decode",
    "rayleigh_pdf",
    "ber_bpsk_rayleigh",
    "ber_bpsk_mrc_rayleigh",
    "selection_outage_probability",
    "maximal_ratio_combine",
    "ofdm_modulate",
    "ofdm_demodulate",
]

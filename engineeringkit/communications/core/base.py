"""Result containers for engineeringkit.communications.

Most communication models return plain arrays or floats. Channel decoding
returns more than the data alone (which words were corrected, and where), so
it gets a frozen dataclass with the usual ``summary()``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["HammingDecodeResult"]


@dataclass(frozen=True)
class HammingDecodeResult:
    """Output of :func:`engineeringkit.communications.hamming74_decode`."""

    data: np.ndarray
    """ndarray of int: Decoded data bits, four per codeword."""

    syndromes: np.ndarray
    """ndarray of int: Per-word syndrome, the 1-based position of the corrected bit (0 if none)."""

    corrected: np.ndarray
    """ndarray of int: Codewords after correction."""

    def summary(self) -> dict:
        """Number of words decoded and number in which a bit was corrected."""
        return {"words": int(self.syndromes.size), "corrected_words": int(np.count_nonzero(self.syndromes))}

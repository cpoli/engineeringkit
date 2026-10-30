Examples
========

.. include:: /_generated/nav/communications.rst

This gallery walks through every public feature of
``engineeringkit.communications``, one breakthrough per script: decibels,
radio propagation and link budgets, antennas, thermal noise, analog
modulation and its noise performance, information, pulse shaping and
capacity, optimum detection, digitization, source and channel codes,
fading and diversity, and OFDM.

:doc:`Units, assumptions and numerical validation </guides/communications>`
gives the physical validity limits. Each script is self-contained and can be
run directly with ``python examples/communications/<section>/<script>.py``.

Sections
--------

- **decibels** -- the Bell System transmission unit.
- **propagation** -- Hertz's waves, the Friis transmission formula and the Telstar link budget.
- **antennas** -- Friis's effective aperture and Carter's half-wave dipole.
- **noise** -- Johnson-Nyquist noise and Friis's cascaded noise factor.
- **modulation** -- Fessenden's AM, Carson's suppressed carrier, Armstrong's FM and Crosby's noise analysis of FM against AM.
- **information** -- Nyquist's signalling rate and raised-cosine pulses, Hartley's law, Shannon's entropy and capacity, and Huffman's code.
- **detection** -- Kotel'nikov's optimum BPSK receiver, quadrature (QPSK) modulation, North's matched filter, Cahn's M-PSK and square QAM.
- **digitization** -- Reeves's pulse-code modulation and delta modulation.
- **coding** -- the Hamming(7,4) code, Elias's convolutional codes, the cyclic redundancy check and Viterbi decoding.
- **fading** -- Rayleigh's random phasors and diversity combining.
- **multiplexing** -- DFT-based OFDM with a cyclic prefix.

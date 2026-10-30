Examples
========

.. include:: /_generated/nav/signal_processing.rst

This gallery walks through every public feature of
``engineeringkit.signal_processing``, one breakthrough per script: Fourier
series, the FFT, the DCT and the Hilbert transform, spectral estimation
and autoregressive models, sampling, quantization and rate conversion,
digital filter design, optimal and adaptive filters, and wavelets.

:doc:`Units, assumptions and numerical validation </guides/signal_processing>`
gives the validity limits. Each script is self-contained and can be run
directly with ``python examples/signal_processing/<section>/<script>.py``.

Sections
--------

- **spectral** -- Schuster's periodogram, the Wiener-Khinchin theorem, Gabor's short-time transform, Blackman and Tukey's windows, Goertzel's algorithm, Welch's averaged periodogram, Yule's autoregressive model, Durbin's recursion and the cepstrum.
- **fourier** -- Gibbs's phenomenon, the Cooley-Tukey FFT, Hilbert's transform and the analytic signal, Stockham's fast convolution and the DCT.
- **wavelets** -- Haar's wavelet.
- **filters** -- Butterworth's filter, Hurewicz's z-transform, Wiener's optimal filter, the LMS adaptive filter, Kaiser's FIR design, second-order sections, Chebyshev and elliptic filters, bilinear design with prewarping and Parks-McClellan equiripple design.
- **sampling** -- Bennett's quantization noise, Shannon's sampling theorem and Schafer and Rabiner's interpolation and decimation.

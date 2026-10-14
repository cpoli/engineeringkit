Breakthroughs in Signal Processing
==================================

.. include:: /_generated/nav/signal_processing.rst

These landmarks trace a conceptual development from a physical law to a useful
engineering model. The order is a learning path, not a dated chronology or a
claim that any idea had a single inventor. SI units apply unless documented
otherwise.

Fourier: a waveform has a frequency description
-----------------------------------------------

Decomposing a signal into complex exponentials makes filtering and spectral analysis tractable. A finite record and a window alter the displayed spectrum.

.. math::

    X_k=\sum_{n=0}^{N-1}x_n e^{-j2\pi kn/N}

Use :py:obj:`engineeringkit.signal_processing.dft` to explore the model.
FFT bins are a finite-record representation; leakage and normalization must be considered.

Sampling: information can be lost before filtering
--------------------------------------------------

Sampling creates indistinguishable frequencies separated by integer multiples of the sample rate. Band-limited reconstruction requires an adequate rate and antialias filtering before digitization.

.. math::

    f_{Nyquist}=f_s/2

Use :py:obj:`engineeringkit.signal_processing.alias_frequency` to explore the model.
Finite sinc reconstruction is truncated; samples cannot reveal an already aliased original frequency.

Digital filtering: system memory is explicit
--------------------------------------------

Convolution connects an input to an impulse response. FIR and IIR implementations trade finite memory against recursive structure and phase behavior.

.. math::

    y[n]=\sum_k h[k]x[n-k]

Use :py:obj:`engineeringkit.signal_processing.filter_signal` to explore the model.
These helpers use causal forward filtering with default initial state, so startup transients and phase delay remain.

Explore and verify
------------------

:doc:`/guides/signal_processing` records units, conventions and independent checks.
Compare the gallery output against the equations before extending a model.

.. minigallery:: ../../examples/signal_processing/plot_signal_processing_models.py ../../examples/signal_processing/plot_signal_processing_benchmark.py

Further reading
---------------

`University or agency teaching reference <https://ocw.mit.edu/courses/res-6-007-signals-and-systems-spring-2011/pages/lecture-notes/>`_
provides background beyond the implemented subset. The equations and examples
here are concise original explanations, not reproductions of that reference.

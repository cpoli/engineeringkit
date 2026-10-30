Breakthroughs in Communications
===============================


.. include:: /_generated/nav/communications.rst

.. epigraph::

   "The fundamental problem of communication is that of reproducing at one
   point either exactly or approximately a message selected at another point."
   -- Claude E. Shannon, 1948

Electrical communication grew from Hertz's laboratory sparks to a theory
precise enough to say how many bits per second any channel can carry.
Telephone and radio engineers at Bell Labs and elsewhere supplied most of the
tools in this chronology: modulation, the decibel, the noise figure and the
link budget. Shannon's information theory then turned the field from an
art into a science with a hard limit. The entries trace the ideas behind
:mod:`engineeringkit.communications`. Each links to the implementing function
and to a gallery example devoted to it.

.. contents:: Timeline
   :local:
   :depth: 1

1880 -- Rayleigh's Sum of Random Phasors
----------------------------------------

Lord Rayleigh asked, for optics and acoustics, what amplitude results from
adding a large number of vibrations of the same pitch and arbitrary phase.
He found that the resultant's amplitude has the density
:math:`(2r/\Omega)e^{-r^2/\Omega}` and its intensity is exponentially
distributed, with mean equal to the sum of the components' intensities.
The same law describes a radio signal reaching a receiver by many scattered
paths. Its envelope fades deeply and often, and the average bit-error rate
of BPSK falls only as :math:`1/4\bar\gamma` instead of exponentially.

*Implementation:* :func:`engineeringkit.communications.rayleigh_pdf` and
:func:`engineeringkit.communications.ber_bpsk_rayleigh`. The example adds
50 random phasors, recovers the Rayleigh density, and measures the extra
power a fading link needs.

*References:* Lord Rayleigh, "On the resultant of a large number of
vibrations of the same pitch and of arbitrary phase," *Phil. Mag.* (5) 10,
no. 60 (1880) 73-78.

.. minigallery:: ../../examples/communications/fading/plot_01_rayleigh_random_phasors.py

1887-1888 -- Hertz Generates and Measures Radio Waves
-----------------------------------------------------

James Clerk Maxwell's equations (1865) predicted electromagnetic waves
travelling at the speed of light, but no one had produced or detected them.
Heinrich Hertz drove a spark-gap dipole at a few hundred megahertz and
detected the radiation with a resonant loop across the room. A metal sheet
reflected the waves into a standing pattern, and the spacing of the nodes,
half a wavelength, multiplied by the oscillation frequency gave the speed of
light. Hertz thought the waves had no practical use. Within a decade
Marconi was sending telegraph messages with them.

*Implementation:* :func:`engineeringkit.communications.wavelength`. The
example reconstructs the standing-wave measurement and maps wavelengths
across the radio spectrum.

*References:* H. Hertz, "Über Strahlen elektrischer Kraft," *Annalen der
Physik* 36 (1888) 769-783.

.. minigallery:: ../../examples/communications/propagation/plot_01_hertz_wavelength.py

1906 -- Fessenden's Amplitude Modulation
----------------------------------------

Spark transmitters could send only on-off Morse code. Reginald Fessenden
used a high-frequency alternator to produce a continuous carrier and varied
its amplitude with a microphone. On Christmas Eve 1906 he is said to have
broadcast voice and violin music from Brant Rock, Massachusetts, to ships at
sea. Conventional AM, :math:`A_c[1 + \mu m(t)]\cos 2\pi f_c t`, keeps the
carrier so that a simple envelope detector (later a crystal set) can
recover the message, at the cost of spending most of the power on the
carrier.

*Implementation:* :func:`engineeringkit.communications.am_signal`. The
example shows envelope detection succeeding at :math:`\mu = 0.7`, failing
under overmodulation, and the one-third limit on sideband power.

*References:* R. A. Fessenden, "Wireless telephony," *Trans. AIEE* 27
(1908) 553-629.

.. minigallery:: ../../examples/communications/modulation/plot_01_fessenden_amplitude_modulation.py

1915 -- Carson's Suppressed-Carrier and Single-Sideband Modulation
------------------------------------------------------------------

John Renshaw Carson of AT&T analysed modulation mathematically and saw that
the AM carrier conveys no information: the message lives entirely in the
two sidebands, each a mirror of the other. His 1915 patent application
suppressed the carrier, and then one sideband, halving the bandwidth and
saving most of the power. That made carrier-frequency multiplexing of many
telephone calls on one long-distance line practical. The receiver must
regenerate a coherent carrier, and any phase error scales the recovered
message by :math:`\cos\phi`.

*Implementation:* :func:`engineeringkit.communications.dsb_sc`. The example
demodulates coherently with phase errors of 0°, 45° and 90° and shows the
missing carrier line.

*References:* J. R. Carson, "Method and means for signaling with high
frequency waves," U.S. Patent 1,449,382 (filed 1915, issued 1923).

.. minigallery:: ../../examples/communications/modulation/plot_02_carson_suppressed_carrier.py

1924-1929 -- The Decibel
------------------------

Telephone engineers once measured loss in "miles of standard cable", a unit
that depended on frequency. In 1924 the Bell System adopted the
*transmission unit*, ten times the common logarithm of a power ratio, and in
1929 W. H. Martin announced its new name, the decibel, after Alexander
Graham Bell. Because logarithms turn products into sums, the gains and
losses along a chain of cables, repeaters and antennas can be budgeted by
simple addition.

*Implementation:* :func:`engineeringkit.communications.db`,
:func:`engineeringkit.communications.from_db`,
:func:`engineeringkit.communications.dbm` and
:func:`engineeringkit.communications.from_dbm`. The example checks that the
decibel sum along a repeater chain equals the linear product.

*References:* W. H. Martin, "DeciBel -- the name for the transmission unit,"
*Bell Syst. Tech. J.* 8 (1929) 1-2.

.. minigallery:: ../../examples/communications/decibels/plot_01_bell_decibel.py

1928 -- Nyquist's Signalling Rate
---------------------------------

How fast can telegraph pulses be sent through a line of limited bandwidth?
Harry Nyquist showed that a channel of bandwidth :math:`B` passes at most
:math:`2B` independent pulses per second. At that rate suitably shaped
pulses, in the ideal case :math:`\operatorname{sinc}(2Bt)`, vanish at every
other pulse's sampling instant, so there is no intersymbol interference.
The same paper contains the seed of the sampling theorem and the Nyquist
criterion for pulse shaping still used in every modem.

*Implementation:* :func:`engineeringkit.communications.nyquist_symbol_rate`
and :func:`engineeringkit.communications.symbol_rate`, and the
raised-cosine family that meets the criterion with excess bandwidth,
:func:`engineeringkit.communications.raised_cosine_pulse`,
:func:`engineeringkit.communications.raised_cosine_spectrum`,
:func:`engineeringkit.communications.root_raised_cosine_pulse` and
:func:`engineeringkit.communications.raised_cosine_bandwidth`, plotted with
:func:`engineeringkit.communications.visualizers.plot_eye_diagram`. The
first example sends BPSK sinc pulses at and above :math:`2B` and measures
the interference; the second checks the folded-spectrum criterion and opens
the eye of a root-raised-cosine link.

*References:* H. Nyquist, "Certain topics in telegraph transmission
theory," *Trans. AIEE* 47 (1928) 617-644.

.. minigallery:: ../../examples/communications/information/plot_01_nyquist_signalling_rate.py ../../examples/communications/information/plot_05_nyquist_raised_cosine.py

1928 -- Hartley's Law
---------------------

Ralph Hartley stripped "information" of meaning and defined it by choice. A
sequence of :math:`n` symbols, each drawn from :math:`M` distinguishable
values, singles out one of :math:`M^n` possibilities and so carries
:math:`n\log M` units. Combined with Nyquist's :math:`2B` symbols per
second, the line rate is :math:`2B\log_2 M` bits per second. Hartley
recognised that noise limits how many levels can be distinguished, but it
took Shannon to say exactly how.

*Implementation:* :func:`engineeringkit.communications.hartley_capacity`.
The example sends multilevel symbols through Gaussian noise and shows the
error rate rising as the levels crowd together.

*References:* R. V. L. Hartley, "Transmission of information," *Bell Syst.
Tech. J.* 7 (1928) 535-563.

.. minigallery:: ../../examples/communications/information/plot_02_hartley_law.py

1928 -- Johnson-Nyquist Thermal Noise
-------------------------------------

John B. Johnson at Bell Labs found a random voltage across every resistor,
proportional to its resistance, its absolute temperature and the measurement
bandwidth, and independent of its material. In the same issue of *Physical
Review*, Nyquist derived the result from thermodynamics and the equipartition
of energy among the modes of a transmission line. The available noise
power, :math:`kTB` or -174 dBm/Hz at 290 K, is the floor against which
every receiver's sensitivity is measured.

*Implementation:* :func:`engineeringkit.communications.thermal_noise_power`.
The example synthesises resistor noise, recovers :math:`kTB` and the
open-circuit voltage :math:`\sqrt{4kTRB}`, and checks that the spectrum is
flat.

*References:* J. B. Johnson, "Thermal agitation of electricity in
conductors," *Phys. Rev.* 32 (1928) 97-109; H. Nyquist, "Thermal agitation
of electric charge in conductors," *Phys. Rev.* 32 (1928) 110-113.

.. minigallery:: ../../examples/communications/noise/plot_01_johnson_nyquist_noise.py

1932 -- Carter's Half-Wave Dipole
---------------------------------

The half-wave dipole is the reference antenna of radio engineering. With a
sinusoidal current along the wire its far field varies as
:math:`\cos(\tfrac\pi2\cos\theta)/\sin\theta`, and integrating the
radiated power gives a directivity of 1.64 (2.15 dBi) and a radiation
resistance of about 73 ohms. Philip S. Carter of RCA put such calculations
into circuit form. His induced-EMF method gave the self-impedance of a
dipole and the mutual impedance between dipoles, from which arrays could be
designed like networks. The 73-ohm value itself was computed in the 1920s
from the far-field power. Its first appearance could be checked only in
textbooks, not in a primary source, so this entry credits Carter only with
the circuit method.

*Implementation:* :func:`engineeringkit.communications.half_wave_dipole_pattern`,
:func:`engineeringkit.communications.half_wave_dipole_directivity`,
:func:`engineeringkit.communications.half_wave_dipole_resistance` and
:func:`engineeringkit.communications.directivity`, plotted with
:func:`engineeringkit.communications.visualizers.plot_antenna_pattern`. The
example integrates the pattern for the directivity, radiation resistance and
78° beamwidth.

*References:* P. S. Carter, "Circuit relations in radiating systems and
applications to antenna problems," *Proc. IRE* 20, no. 6 (1932) 1004-1041.
The issue and date are checked against later citations; the page range
could not be confirmed.

.. minigallery:: ../../examples/communications/antennas/plot_02_carter_half_wave_dipole.py

1933-1936 -- Armstrong's Wideband FM
------------------------------------

The prevailing view, argued by Carson in 1922, was that frequency
modulation offered no bandwidth saving over AM and so no advantage. Edwin
Howard Armstrong turned the argument around. He used a frequency deviation
much *larger* than the audio bandwidth, and a limiter and discriminator in
the receiver, to trade bandwidth for a dramatic suppression of noise and
static. The spectrum of tone-modulated FM is a comb of sidebands with
Bessel-function amplitudes. Carson's own rule, :math:`B \approx 2(\Delta f
+ W)`, gives the bandwidth that holds about 98% of the power.

*Implementation:* :func:`engineeringkit.communications.fm_signal` and
:func:`engineeringkit.communications.carson_bandwidth`. The example matches
measured spectral lines to :math:`|J_n(\beta)|` and sizes broadcast FM.

*References:* E. H. Armstrong, "A method of reducing disturbances in radio
signaling by a system of frequency modulation," *Proc. IRE* 24 (1936)
689-740; J. R. Carson, "Notes on the theory of modulation," *Proc. IRE* 10
(1922) 57-64.

.. minigallery:: ../../examples/communications/modulation/plot_03_armstrong_frequency_modulation.py

1937 -- Crosby's Noise Analysis of FM and AM
--------------------------------------------

Armstrong's demonstrations showed that wideband FM was quieter than AM, but
not by how much, or when. Murray G. Crosby of RCA worked it out. Measured
against the SNR :math:`\gamma = P_R/N_0W` of baseband transmission with the
same power, an envelope-detected AM receiver delivers at most
:math:`\gamma/3`, because two thirds of its power sits in the carrier. A
limiter-discriminator delivers :math:`3\beta^2 P_m \gamma`, so each doubling
of the deviation buys 6 dB. Crosby also found the threshold: when the
carrier-to-noise ratio in the FM band falls below about 10 dB, noise spikes
appear and the advantage collapses.

*Implementation:* :func:`engineeringkit.communications.baseband_snr`,
:func:`engineeringkit.communications.am_output_snr`,
:func:`engineeringkit.communications.fm_output_snr` and
:func:`engineeringkit.communications.fm_threshold_snr`. The example simulates
an envelope detector and a discriminator over a range of SNR and finds the
threshold.

*References:* M. G. Crosby, "Frequency modulation noise characteristics,"
*Proc. IRE* 25, no. 4 (1937) 472-514.

.. minigallery:: ../../examples/communications/modulation/plot_04_crosby_fm_noise.py

1937-1938 -- Reeves's Pulse-Code Modulation
-------------------------------------------

Alec H. Reeves of International Telephone and Telegraph's Paris laboratory
proposed sampling speech, rounding each sample to one of :math:`2^N`
levels, and sending the level as :math:`N` on-off pulses. Rounding adds a
fixed noise, about :math:`6.02N + 1.76` dB below a full-scale sine. In
return a repeater has only to decide whether each pulse is present, so noise
no longer accumulates along a chain of repeaters. Vacuum tubes made PCM
impractical in 1938. It entered the telephone network with the transistor
and the Bell System's T1 carrier of 1962, which used mu-law companding to
keep quiet speech above the quantization noise. The 1937 date of the idea
comes from later accounts. The French filing of 1938 is checked against
citations of the patent, not the patent itself.

*Implementation:* :func:`engineeringkit.communications.pcm_encode`,
:func:`engineeringkit.communications.pcm_decode`,
:func:`engineeringkit.communications.pcm_sqnr_db`,
:func:`engineeringkit.communications.mu_law_compress` and
:func:`engineeringkit.communications.mu_law_expand`. The example measures 6 dB
per bit, compares PCM and analog repeater chains, and shows companding at
low levels.

*References:* A. H. Reeves, French Patent 852,183 (filed 1938, published
1939); U.S. Patent 2,272,070 (1942); B. M. Oliver, J. R. Pierce and C. E.
Shannon, "The philosophy of PCM," *Proc. IRE* 36 (1948) 1324-1331.

.. minigallery:: ../../examples/communications/digitization/plot_01_reeves_pcm.py

1943 -- North's Matched Filter
------------------------------

Dwight O. North at RCA Laboratories asked which receiving filter makes a
known radar pulse stand out best above white noise. His wartime report
answered: the filter whose impulse response is the pulse reversed in time.
Sampled at the end of the pulse its output signal-to-noise ratio is
:math:`2E/N_0`. That depends only on the pulse energy, not its shape, and by
the Cauchy-Schwarz inequality no linear filter does better. The report
circulated privately and was printed only in 1963. The same receiver
appears in Kotel'nikov's thesis and in every digital demodulator.

*Implementation:* :func:`engineeringkit.communications.matched_filter` and
:func:`engineeringkit.communications.matched_filter_snr`. The example finds
a chirp buried in noise and compares the matched filter with random
filters.

*References:* D. O. North, "An analysis of the factors which determine
signal/noise discrimination in pulsed-carrier systems," RCA Laboratories
Report PTR-6C (25 June 1943); reprinted in *Proc. IEEE* 51 (1963)
1016-1027.

.. minigallery:: ../../examples/communications/detection/plot_03_north_matched_filter.py

1944 -- Friis's Noise Figure
----------------------------

Harald Friis gave receiver engineers a single number to compare designs.
The *noise factor* :math:`F` is the ratio of the signal-to-noise ratio at a
device's input to that at its output, for a 290 K source. For stages in
cascade, :math:`F = F_1 + (F_2 - 1)/G_1 + (F_3 - 1)/(G_1 G_2) + \cdots`.
Each later stage's noise is divided by the gain ahead of it, so a
high-gain, low-noise first amplifier sets the noise figure of the whole
receiver.

*Implementation:* :func:`engineeringkit.communications.cascaded_noise_factor`
and :func:`engineeringkit.communications.noise_figure`. The example checks
the formula against a sample-by-sample simulation and compares orderings of
the same three stages.

*References:* H. T. Friis, "Noise figures of radio receivers," *Proc. IRE*
32 (1944) 419-422.

.. minigallery:: ../../examples/communications/noise/plot_02_friis_cascaded_noise.py

1946 -- The Friis Transmission Formula
--------------------------------------

Friis's two-page note wrote the received power between two antennas in
free space as :math:`P_r = P_t G_t G_r (\lambda/4\pi d)^2`. The formula
combines spherical spreading of the transmitted power with the effective
aperture :math:`G\lambda^2/4\pi` of the receiving antenna. It replaced a
collection of ad hoc field-strength formulas, and it is the propagation term
of every radio link budget.

*Implementation:* :func:`engineeringkit.communications.friis_received_power`,
:func:`engineeringkit.communications.free_space_path_loss` and
:func:`engineeringkit.communications.free_space_path_loss_db`, and in
Friis's own area form :func:`engineeringkit.communications.effective_aperture`
and :func:`engineeringkit.communications.aperture_gain`. The first example
verifies the spreading-times-aperture interpretation and the 6 dB loss per
doubling of distance; the second writes a satellite link with effective
areas and shows how dish gain grows with frequency.

*References:* H. T. Friis, "A note on a simple transmission formula," *Proc.
IRE* 34 (1946) 254-256.

.. minigallery:: ../../examples/communications/propagation/plot_02_friis_transmission.py ../../examples/communications/antennas/plot_01_friis_effective_aperture.py

1946-1952 -- Delta Modulation
-----------------------------

Delta modulation replaces PCM's multilevel coder with one comparator. The
encoder keeps a staircase copy of the signal and sends one bit per sample,
1 to step up by :math:`\Delta` and 0 to step down. The idea appears in a
French patent of 1946 by Maurice Deloraine and colleagues at ITT's Paris
laboratory. Frank de Jager of Philips Research gave the first full analysis
in English in 1952. The step size must balance two errors: *slope
overload* when the signal changes faster than :math:`\Delta f_s`, and
*granular noise* when it is flat. The co-inventors of the patent could be
found only in secondary sources, so they are not named here.

*Implementation:* :func:`engineeringkit.communications.delta_modulate`,
:func:`engineeringkit.communications.delta_demodulate` and
:func:`engineeringkit.communications.delta_slope_limit`. The example shows
overload and tracking, and finds the step that minimises the error.

*References:* F. de Jager, "Delta modulation, a method of P.C.M.
transmission using the 1-unit code," *Philips Res. Rep.* 7 (1952) 442-466.

.. minigallery:: ../../examples/communications/digitization/plot_02_delta_modulation.py

1947 -- Kotel'nikov's Theory of Optimum Noise Immunity
------------------------------------------------------

Vladimir Kotel'nikov's doctoral thesis asked which receiver minimises the
probability of error when one of several known signals arrives in white
Gaussian noise. The optimum decides for the signal closest to the received
waveform, which amounts to correlating with each candidate (a matched
filter). For antipodal binary signals, BPSK, the bit-error rate is
:math:`Q(\sqrt{2E_b/N_0})`, the benchmark against which every digital
modulation is still measured. The thesis was published in English in 1959.

*Implementation:* :func:`engineeringkit.communications.bpsk`,
:func:`engineeringkit.communications.ber_bpsk_awgn` and
:func:`engineeringkit.communications.q_function`, plotted with
:func:`engineeringkit.communications.visualizers.plot_ber_curve`. The
example confirms the curve by Monte Carlo simulation.

*References:* V. A. Kotel'nikov, *The Theory of Optimum Noise Immunity*
(1947; English translation, New York: McGraw-Hill, 1959).

.. minigallery:: ../../examples/communications/detection/plot_01_kotelnikov_bpsk_detection.py

1948 -- Shannon's Entropy and the Binary Symmetric Channel
----------------------------------------------------------

"A Mathematical Theory of Communication" defined the information of a
source as its entropy, :math:`H = -\sum p\log_2 p` bits per symbol, the
minimum average length of any lossless encoding. It defined the capacity of
a channel as the largest mutual information between input and output. And
it proved the startling *noisy-channel coding theorem*: below capacity,
codes exist that make the error probability as small as desired. For a
binary symmetric channel that flips each bit with probability :math:`p`,
the capacity is :math:`1 - H(p)`.

*Implementation:* :func:`engineeringkit.communications.entropy_binary` and
:func:`engineeringkit.communications.binary_symmetric_capacity`. The example
measures the empirical entropy of biased bit streams and maximises mutual
information over input distributions.

*References:* C. E. Shannon, "A mathematical theory of communication,"
*Bell Syst. Tech. J.* 27 (1948) 379-423, 623-656.

.. minigallery:: ../../examples/communications/information/plot_03_shannon_entropy.py

1949 -- The Shannon-Hartley Capacity Theorem
--------------------------------------------

In "Communication in the presence of noise", Shannon used a geometric
argument about spheres in signal space to give the capacity of a
band-limited channel with white Gaussian noise: :math:`C = B\log_2(1 + S/N)`
bits per second. The formula sets the exchange rate between bandwidth and
power. As bandwidth grows without limit, the energy per bit needed for
reliable communication falls to :math:`\ln 2`, or -1.59 dB, the *Shannon
limit* that turbo and LDPC codes approached half a century later.

*Implementation:* :func:`engineeringkit.communications.shannon_capacity` and
:func:`engineeringkit.communications.spectral_efficiency`. The example plots
the capacity curve and the minimum :math:`E_b/N_0` for each spectral
efficiency.

*References:* C. E. Shannon, "Communication in the presence of noise,"
*Proc. IRE* 37 (1949) 10-21.

.. minigallery:: ../../examples/communications/information/plot_04_shannon_hartley_capacity.py

1950 -- Hamming's Error-Correcting Codes
----------------------------------------

Richard Hamming's relay computer at Bell Labs detected errors over the
weekend but could only abort the job. Hamming designed codes that *correct*
them. In the (7,4) code, parity bits at positions 1, 2 and 4 each check the
positions whose binary index contains that bit. When a single bit flips,
the pattern of failed checks, the *syndrome*, spells out its position in
binary. Hamming also introduced the distance between codewords: a minimum
distance of three is exactly what single-error correction requires.

*Implementation:* :func:`engineeringkit.communications.hamming74_encode`,
:func:`engineeringkit.communications.hamming74_decode` (returning a
:class:`~engineeringkit.communications.HammingDecodeResult`) and
:func:`engineeringkit.communications.hamming_distance`. The example corrects
every single-bit error and compares coded and uncoded block error rates.

*References:* R. W. Hamming, "Error detecting and error correcting codes,"
*Bell Syst. Tech. J.* 29 (1950) 147-160.

.. minigallery:: ../../examples/communications/coding/plot_01_hamming_code.py

1952 -- Huffman's Minimum-Redundancy Code
-----------------------------------------

Offered a term paper in place of Robert Fano's final exam at MIT, David A.
Huffman was asked to find the most efficient binary code for a source with
known probabilities. Shannon and Fano had split the symbols top-down, which
is not always optimal. Huffman worked bottom-up. He merged the two least
likely symbols into one and repeated until one remained. The resulting
prefix code has the least average length of any symbol-by-symbol code,
within one bit of the entropy, and equal to it when every probability is a
power of one half. Huffman codes still sit inside JPEG, MP3 and DEFLATE.

*Implementation:* :func:`engineeringkit.communications.huffman_code`
(returning a :class:`~engineeringkit.communications.HuffmanCode` that
encodes and decodes) and :func:`engineeringkit.communications.entropy`. The
example codes the letters of English and a dyadic source.

*References:* D. A. Huffman, "A method for the construction of
minimum-redundancy codes," *Proc. IRE* 40 (1952) 1098-1101.

.. minigallery:: ../../examples/communications/information/plot_06_huffman_code.py

1954-1959 -- Kahn and Brennan's Diversity Combining
---------------------------------------------------

Two antennas a few wavelengths apart seldom fade at the same moment.
Short-wave receivers had long switched to whichever antenna was strongest.
Leonard R. Kahn's 1954 "ratio squarer" instead weighted each branch by its
own signal strength before adding them. With weights equal to the
conjugate channel gains this is maximal-ratio combining, whose output SNR is
the sum of the branch SNRs. D. G. Brennan's 1959 review analysed selection,
equal-gain and maximal-ratio combining together and showed how little
separates them. With :math:`L` independent branches the error rate falls as
:math:`\bar\gamma^{-L}`.

*Implementation:* :func:`engineeringkit.communications.maximal_ratio_combine`,
:func:`engineeringkit.communications.ber_bpsk_mrc_rayleigh` and
:func:`engineeringkit.communications.selection_outage_probability`. The
example simulates one, two and four branches and measures the diversity
order and the selection outage.

*References:* L. R. Kahn, "Ratio squarer," *Proc. IRE* 42, no. 11 (1954)
1704; D. G. Brennan, "Linear diversity combining techniques," *Proc. IRE*
47 (1959) 1075-1102.

.. minigallery:: ../../examples/communications/fading/plot_02_diversity_combining.py

1955 -- Elias's Convolutional Codes
-----------------------------------

Peter Elias's "Coding for noisy channels" showed that a code need not cut
the data into blocks. A shift register slides along the bit stream, and each
input bit produces :math:`n` output bits, sums modulo two of the register's
taps. Elias proved that these convolutional codes are at least as good as
block codes. Their error-correcting power is the *free distance*, the least
weight of a path that leaves the all-zero state of the encoder and returns
to it. It is 5 for the constraint-length-3 code with generators 7 and 5
(octal), and 10 for the constraint-length-7 code of deep-space and Wi-Fi
links.

*Implementation:* :func:`engineeringkit.communications.convolutional_encode`
and :func:`engineeringkit.communications.convolutional_free_distance`. The
example encodes with the (7, 5) code, checks linearity, and computes the
free distance of the best rate-1/2 codes.

*References:* P. Elias, "Coding for noisy channels," *IRE Convention
Record* 3, part 4 (1955) 37-46.

.. minigallery:: ../../examples/communications/coding/plot_02_elias_convolutional_codes.py

1959 -- Cahn's M-ary Phase-Shift Keying
---------------------------------------

Charles R. Cahn asked what happens when a carrier's phase takes :math:`M`
values instead of two, so that each symbol carries :math:`\log_2 M` bits.
The points sit on a circle and crowd together as :math:`M` grows: the
distance between neighbours is :math:`2\sin(\pi/M)`, and beyond four
phases each doubling of :math:`M` costs about 6 dB. Cahn considered both
coherent detection and detection by comparison with the previous phase. He
concluded that multiphase signalling trades bandwidth for signal-to-noise
ratio more efficiently than multilevel amplitude signalling. With Gray
labels, neighbouring phases differ in one bit.

*Implementation:* :func:`engineeringkit.communications.psk_constellation`,
:func:`engineeringkit.communications.ser_mpsk_awgn` (Craig's exact integral)
and :func:`engineeringkit.communications.minimum_distance_detect`. The
example simulates 2-, 4-, 8- and 16-PSK.

*References:* C. R. Cahn, "Performance of digital phase-modulation
communication systems," *IRE Trans. Commun. Syst.* 7 (1959) 3-6; J. W.
Craig, "A new, simple and exact result for calculating the probability of
error for two-dimensional signal constellations," *IEEE MILCOM* (1991)
571-575.

.. minigallery:: ../../examples/communications/detection/plot_04_cahn_mpsk.py

1960 -- Quadrature Modulation
-----------------------------

A cosine and a sine carrier at the same frequency are orthogonal, so each
can be modulated independently and separated again by coherent detection.
Charles Cahn's analysis of combined phase and amplitude modulation put
quadrature signalling on a quantitative footing for data transmission.
Quadrature phase-shift keying sends two bits per symbol in the bandwidth of
one, with the same bit-error rate per unit :math:`E_b/N_0` as BPSK. The same
in-phase/quadrature idea underlies every QAM modem.

*Implementation:* :func:`engineeringkit.communications.qpsk` and
:func:`engineeringkit.communications.ber_qpsk_awgn`, plotted with
:func:`engineeringkit.communications.visualizers.plot_constellation`. The
example simulates a noisy QPSK constellation and separates the two rails of
a passband waveform.

*References:* C. R. Cahn, "Combined digital phase and amplitude modulation
communication systems," *IRE Trans. Commun. Syst.* 8 (1960) 150-155.

.. minigallery:: ../../examples/communications/detection/plot_02_cahn_quadrature_qpsk.py

1961 -- Peterson and Brown's Cyclic Redundancy Check
----------------------------------------------------

W. Wesley Peterson and D. T. Brown showed how to detect errors with a few
shift-register stages and exclusive-OR gates. The message is treated as a
polynomial over GF(2), multiplied by :math:`x^r` and divided by a generator
of degree :math:`r`, and the remainder is sent after it. The receiver
divides again and declares an error if anything is left. Every burst of
:math:`r` or fewer errors is caught, and a random error pattern passes with
probability :math:`2^{-r}`. Cyclic redundancy checks now guard Ethernet
frames, disk sectors and USB packets.

*Implementation:* :func:`engineeringkit.communications.crc_remainder`,
:func:`engineeringkit.communications.crc_encode`,
:func:`engineeringkit.communications.crc_check` and the CRC-16-CCITT
generator :data:`engineeringkit.communications.CRC16_CCITT`. The example
reproduces the catalogue check value and tests every short burst.

*References:* W. W. Peterson and D. T. Brown, "Cyclic codes for error
detection," *Proc. IRE* 49, no. 1 (1961) 228-235.

.. minigallery:: ../../examples/communications/coding/plot_03_peterson_brown_crc.py

1962 -- Campopiano and Glazer's Square QAM
------------------------------------------

Phase-shift keying leaves the inside of the circle empty. C. N. Campopiano
and B. G. Glazer modulated the amplitudes of two carriers in quadrature
independently, so the signal points fill a square grid. Each rail is a
:math:`\sqrt M`-level amplitude decision, and the symbol-error rate follows
exactly from one rail's error probability. At 16 points and above, square
QAM needs several decibels less power than PSK for the same error rate, and
it has carried data in voiceband modems, cable, Wi-Fi and cellular systems
ever since.

*Implementation:* :func:`engineeringkit.communications.qam_constellation`
and :func:`engineeringkit.communications.ser_mqam_awgn`. The example
simulates 4-, 16- and 64-QAM and measures the 4 dB advantage of 16-QAM over
16-PSK.

*References:* C. N. Campopiano and B. G. Glazer, "A coherent digital
amplitude and phase modulation scheme," *IRE Trans. Commun. Syst.* 10,
no. 1 (1962) 90-95.

.. minigallery:: ../../examples/communications/detection/plot_05_campopiano_glazer_qam.py

1962 -- Telstar and the Satellite Link Budget
---------------------------------------------

Telstar 1, built by Bell Labs and launched in July 1962, relayed the first
live transatlantic television between Andover, Maine and stations in Europe.
Its 4 GHz downlink transmitted only a few watts, and the signal was received
by a large horn antenna with a maser preamplifier cooled to a few kelvin.
Such systems were designed with a *link budget*: transmitter power, antenna
gains, free-space and other losses added in decibels, then compared with
the receiver noise floor :math:`kTB` to give the signal-to-noise ratio.

*Implementation:* :func:`engineeringkit.communications.link_budget_db`
together with the propagation and noise functions. The example budgets a
representative Telstar downlink over a pass.

*References:* A. C. Dickieson et al., "The Telstar experiment," *Bell Syst.
Tech. J.* 42 (1963) 739-1907 (special issue).

.. minigallery:: ../../examples/communications/propagation/plot_03_telstar_link_budget.py

1966-1971 -- Orthogonal Frequency-Division Multiplexing
-------------------------------------------------------

Robert Chang showed that many data subcarriers can overlap in frequency yet
remain separable if they are spaced by exactly the symbol rate. Stephen
Weinstein and Paul Ebert generated all the subcarriers with a single
inverse discrete Fourier transform and recovered them with a forward DFT.
They also prefixed each block with a copy of its tail, a *cyclic prefix*
that makes a multipath channel's convolution circular. Each subcarrier then
sees a single complex gain, removed by one division. OFDM became practical
with fast DSP chips and now carries Wi-Fi, LTE, 5G and digital
broadcasting.

*Implementation:* :func:`engineeringkit.communications.ofdm_modulate` and
:func:`engineeringkit.communications.ofdm_demodulate`. The example passes
QPSK-loaded OFDM through a multipath channel with and without a cyclic
prefix.

*References:* R. W. Chang, "Synthesis of band-limited orthogonal signals for
multichannel data transmission," *Bell Syst. Tech. J.* 45 (1966) 1775-1796;
S. B. Weinstein and P. M. Ebert, "Data transmission by frequency-division
multiplexing using the discrete Fourier transform," *IEEE Trans. Commun.
Technol.* 19 (1971) 628-634.

.. minigallery:: ../../examples/communications/multiplexing/plot_01_ofdm_cyclic_prefix.py

1967 -- Viterbi's Decoding Algorithm
------------------------------------

Decoding a convolutional code by brute force means comparing the received
sequence with every possible codeword. Andrew J. Viterbi saw that two paths
entering the same encoder state at the same time share the same future, so
only the better of them, the *survivor*, need be kept. The work then grows
linearly with the message length. Viterbi introduced the algorithm to prove
error bounds. J. K. Omura and G. D. Forney soon showed that it is
maximum-likelihood decoding. It went on to decode deep-space telemetry, GSM
telephones and dial-up modems. Using the received values before slicing
(soft decisions) gains about 2 dB over hard bits.

*Implementation:* :func:`engineeringkit.communications.viterbi_decode`
(returning a :class:`~engineeringkit.communications.ViterbiDecodeResult`).
The example checks the survivor metric against exhaustive search and
measures the coding gain of the constraint-length-7 code.

*References:* A. J. Viterbi, "Error bounds for convolutional codes and an
asymptotically optimum decoding algorithm," *IEEE Trans. Inf. Theory* 13
(1967) 260-269.

.. minigallery:: ../../examples/communications/coding/plot_04_viterbi_algorithm.py

See Also
--------

- :doc:`/api/communications`
- :doc:`/guides/communications`
- :doc:`/history/signal_processing_breakthroughs`
- :doc:`/history/electrical_engineering_breakthroughs`

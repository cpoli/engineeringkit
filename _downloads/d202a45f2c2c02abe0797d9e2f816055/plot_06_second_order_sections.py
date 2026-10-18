r"""
Jackson's cascaded second-order sections
========================================

A digital filter can be computed from its transfer function
:math:`B(z)/A(z)` in one "direct form", but the roots of a high-degree
polynomial are extremely sensitive to rounding in its coefficients. When
the poles crowd near :math:`z = 1`, as they do for a narrow low-pass
filter, double precision is not enough. Kaiser warned of this in 1965, and
Leland Jackson (1970) analysed how roundoff noise and overflow depend on
how a filter is split into a cascade of second-order sections, each
holding one pair of poles. Cascades of "biquads" have been the standard
realization ever since
(:func:`~engineeringkit.signal_processing.butterworth_lowpass` with
``output="sos"``, :func:`~engineeringkit.signal_processing.filter_signal`).
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.signal_processing import butterworth_lowpass, filter_signal, frequency_response
from engineeringkit.signal_processing.visualizers import plot_magnitude_response

# %%
# Twelfth order, cutoff 1% of the sampling rate
# ---------------------------------------------

fs, fc, order = 1000.0, 10.0, 12
b, a = butterworth_lowpass(order, fc, fs)
sos = butterworth_lowpass(order, fc, fs, output="sos")
print(f"{sos.shape[0]} sections of [b0 b1 b2 1 a1 a2]; direct form has {a.size} denominator coefficients")
print(f"direct-form coefficients span {abs(a).min():.1e} to {abs(a).max():.1e}")

# %%
# Rounding moves the poles out of the unit circle
# -----------------------------------------------
#
# The exact poles lie on a circle of radius below 1. Recomputed from the
# rounded polynomial ``a``, some land outside it: the direct form is
# unstable. Each section's quadratic keeps its own two poles accurately.

poles_direct = np.roots(a)
poles_sections = np.concatenate([np.roots(section[3:]) for section in sos])
print(f"largest pole radius: direct form {abs(poles_direct).max():.4f}, sections {abs(poles_sections).max():.4f}")
assert abs(poles_direct).max() > 1 and abs(poles_sections).max() < 1

# %%
# Response and step response
# --------------------------

_, H_dc_sections = frequency_response(sos, None, worN=[0.0, fc], fs=fs)
_, H_dc_direct = frequency_response(b, a, worN=[0.0, fc], fs=fs)
print(f"gain at DC and at fc: sections {np.round(abs(H_dc_sections), 4)}, direct form {np.round(abs(H_dc_direct), 4)}")
assert np.allclose(abs(H_dc_sections), [1, 1 / np.sqrt(2)])

step = np.ones(3000)
y_sections = filter_signal(sos, None, step)
y_direct = filter_signal(b, a, step)
print(f"step response after 3 s: sections {y_sections[-1]:.6f}, direct form {y_direct[-1]:.3e}")
assert abs(y_sections[-1] - 1) < 1e-6 and abs(y_direct[-1]) > 1e6

# %%
# Low orders are safe either way
# ------------------------------

b4, a4 = butterworth_lowpass(4, fc, fs)
sos4 = butterworth_lowpass(4, fc, fs, output="sos")
assert np.allclose(filter_signal(b4, a4, step), filter_signal(sos4, None, step))

fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
plot_magnitude_response(sos, None, fs=fs, ax=axes[0], label="sections")
plot_magnitude_response(b, a, fs=fs, ax=axes[0], label="direct form")
axes[0].set(xlim=(0, 50), ylim=(-150, 60), title=f"Order {order}, $f_c$ = {fc:.0f} Hz")
theta = np.linspace(0, 2 * np.pi, 361)
axes[1].plot(np.cos(theta), np.sin(theta), color="gray", lw=0.8)
axes[1].plot(poles_sections.real, poles_sections.imag, "x", label="sections")
axes[1].plot(poles_direct.real, poles_direct.imag, "+", label="direct form")
axes[1].set(xlim=(0.85, 1.15), ylim=(-0.15, 0.15), aspect="equal", title="Poles near z = 1")
axes[1].legend()
axes[2].plot(y_sections, label="sections")
axes[2].plot(np.clip(y_direct, -2, 2), label="direct form (clipped)")
axes[2].set(xlabel="sample", ylabel="step response", ylim=(-2.2, 2.2))
axes[2].legend()
fig.tight_layout()

plt.show()

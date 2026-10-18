r"""
Bode's gain and phase margins
=============================

Hendrik Bode (1940, 1945) turned Nyquist's criterion into design numbers
read off a logarithmic plot of gain and phase
(:func:`~engineeringkit.control_systems.bode`). The *gain margin* is how
much the loop gain can grow before the phase-crossover point reaches
:math:`-1`; the *phase margin* is how much extra lag the loop tolerates at
the gain crossover (:func:`~engineeringkit.control_systems.stability_margins`).
For :math:`L = 2/(s(s+1)(s+2))` the phase crosses :math:`-180^\circ` at
:math:`\omega = \sqrt2`, where :math:`|L| = 1/3`, so the gain margin is 3.
A loop can also cross :math:`-180^\circ` more than once. Then the margin that
matters is the crossover nearest instability, which need not be the first.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from engineeringkit.control_systems import feedback, poles, stability_margins, transfer_function
from engineeringkit.control_systems.visualizers import plot_bode

loop = transfer_function([2], [1, 3, 2, 0])
m = stability_margins(loop)
print({k: round(m[k], 4) for k in ("gain_margin", "phase_margin_deg", "gain_cross_frequency", "phase_cross_frequency")})

# %%
# Independent checks
# ------------------

wc = brentq(lambda w: 2 / (w * np.hypot(w, 1) * np.hypot(w, 2)) - 1, 0.1, 2)
pm = 180 - np.degrees(np.pi / 2 + np.arctan(wc) + np.arctan(wc / 2))
assert np.isclose(m["gain_margin"], 3, rtol=1e-4) and np.isclose(m["phase_cross_frequency"], np.sqrt(2), rtol=1e-4)
assert np.isclose(m["gain_cross_frequency"], wc, rtol=1e-4) and np.isclose(m["phase_margin_deg"], pm, atol=0.01)
# Multiplying the gain by the gain margin puts closed-loop poles on the axis.
edge = poles(feedback(transfer_function([2 * m["gain_margin"]], [1, 3, 2, 0])))
assert np.isclose(np.abs(edge.imag).max(), np.sqrt(2), rtol=1e-3) and abs(edge.real.max()) < 1e-3

ax_mag, ax_phase = plot_bode(loop, np.logspace(-2, 2, 500))
ax_mag.axhline(0, color="k", lw=0.8)
ax_phase.axhline(-180, color="k", lw=0.8)
ax_mag.vlines(m["phase_cross_frequency"], -20 * np.log10(m["gain_margin"]), 0, color="C3", lw=2, label=f"GM = {m['gain_margin']:.2f}")
pm_label = f"PM = {m['phase_margin_deg']:.1f} deg"
ax_phase.vlines(m["gain_cross_frequency"], -180, -180 + m["phase_margin_deg"], color="C2", lw=2, label=pm_label)
ax_mag.legend()
ax_phase.legend()

# %%
# A conditionally stable loop
# ---------------------------
#
# :math:`L = 5(s+1)^2/(s^3(0.1s+1)^2)` has three integrators, so its phase
# starts at :math:`-270^\circ`. The lead of the double zero lifts it above
# :math:`-180^\circ` and the lag of the double pole pulls it back, so the
# phase crosses twice. The closed loop is stable, but *lowering* the gain by
# the factor at the first crossover destabilises it as surely as raising it
# by the factor at the second. The worst margin is the one closer to 1 on a
# logarithmic scale, here the upper one.

den = np.polymul([1, 0, 0, 0], np.polymul([0.1, 1], [0.1, 1]))
conditional = transfer_function(5 * np.polymul([1, 1], [1, 1]), den)
mc = stability_margins(conditional)
for wp, gm in zip(mc["phase_cross_frequencies"], mc["gain_margins"]):
    print(f"phase crossover {wp:.3f} rad/s: gain margin {gm:.3f}")
print(f"worst: {mc['gain_margin']:.3f} at {mc['phase_cross_frequency']:.3f} rad/s")
assert np.all(poles(feedback(conditional)).real < 0)
assert mc["gain_margins"].size == 2 and mc["gain_margins"][0] < 1 < mc["gain_margins"][1]
assert mc["gain_margin"] == mc["gain_margins"][1]
# Each margin, applied to the gain, puts a closed-loop pole pair on the imaginary axis.
for wp, gm in zip(mc["phase_cross_frequencies"], mc["gain_margins"]):
    edge = poles(feedback(transfer_function(5 * gm * np.polymul([1, 1], [1, 1]), den)))
    j = np.argmin(np.abs(edge - 1j * wp))
    assert abs(edge[j].real) < 1e-3 and np.isclose(edge[j].imag, wp, rtol=1e-3)

ax_mag, ax_phase = plot_bode(conditional, np.logspace(-1, 2.5, 500))
ax_mag.axhline(0, color="k", lw=0.8)
ax_phase.axhline(-180, color="k", lw=0.8)
for wp, gm in zip(mc["phase_cross_frequencies"], mc["gain_margins"]):
    ax_mag.vlines(wp, -20 * np.log10(gm), 0, color="C3", lw=2)
    ax_mag.annotate(f"GM = {gm:.2f}", (wp, -20 * np.log10(gm)), textcoords="offset points", xytext=(5, 0))

plt.show()

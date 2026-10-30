r"""
MIL-P-1629's failure mode and effects analysis (1949)
=====================================================

The US military procedure MIL-P-1629, *Procedures for Performing a Failure
Mode, Effects and Criticality Analysis* (9 November 1949), asked designers
to go through a system part by part, list each way the part can fail,
trace the effect of each failure and classify it by its consequence for
the mission and for safety. NASA used the method on Apollo, and the
automotive industry later made it routine, ranking each mode by a *risk
priority number*, severity times occurrence times detection, each scored
from 1 to 10 (:func:`~engineeringkit.reliability.risk_priority_number`,
:func:`~engineeringkit.reliability.fmea_ranking`). Where failure rates are
known the military standard's criticality number
:math:`C_m = \beta\alpha\lambda_p t`
(:func:`~engineeringkit.reliability.failure_mode_criticality`) ranks
modes by expected losses instead.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.reliability import failure_mode_criticality, fmea_ranking

# %%
# An FMEA of a centrifugal pump
# -----------------------------
#
# (severity, occurrence, detection) for each failure mode.

modes = {
    "mechanical seal leak": (6, 6, 3),
    "bearing seizure": (8, 3, 5),
    "impeller erosion": (5, 4, 7),
    "motor winding short": (9, 2, 4),
    "coupling misalignment": (4, 5, 5),
    "cavitation": (5, 5, 6),
    "casing crack": (10, 1, 6),
}
fmea = fmea_ranking(modes)
for name, s, o, d, rpn in zip(fmea.modes, fmea.severity, fmea.occurrence, fmea.detection, fmea.rpn):
    print(f"{name:22s} S={s:2d} O={o:2d} D={d:2d} RPN={rpn:4d}")
fig, ax = fmea.plot()
fig.set_size_inches(7, 4)
fig.tight_layout()

# %%
# What the RPN hides
# ------------------
#
# The casing crack has the lowest RPN but a severity of 10. Because the
# scores are ordinal, a high severity is reviewed whatever its product.

critical = [name for name, s in zip(fmea.modes, fmea.severity) if s >= 9]
print("modes reviewed for severity alone:", critical)

# %%
# Criticality numbers from failure rates
# --------------------------------------
#
# The seal of a pump with :math:`\lambda_p` = 50 failures per million hours
# fails by leaking (60% of failures, of which 20% stop the mission) or by
# sudden rupture (40%, all of which stop it); the mission is 2000 h.

leak, rupture = failure_mode_criticality(np.array([0.6, 0.4]), np.array([0.2, 1.0]), 50e-6, 2000.0)
print(f"criticality: leak {leak:.3f}, rupture {rupture:.3f} expected losses per mission")

assert fmea.modes[0] == "cavitation" and fmea.rpn[0] == 150
assert "casing crack" in critical and fmea.modes[-1] == "casing crack"
assert np.isclose(leak + rupture, 50e-6 * 2000 * (0.6 * 0.2 + 0.4))

plt.show()

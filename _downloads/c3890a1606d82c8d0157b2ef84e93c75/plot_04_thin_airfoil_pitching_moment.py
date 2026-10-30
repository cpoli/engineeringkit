r"""
Thin-airfoil pitching moment and the aerodynamic centre (1922-1926)
===================================================================

Munk's and Glauert's thin-airfoil theory predicts more than lift. The
pitching moment about the quarter chord,
:math:`c_{m,c/4} = \tfrac{\pi}{4}(A_2 - A_1)`
(:func:`~engineeringkit.aerodynamics.naca4_moment_coefficient`), depends
only on the camber line and not on the angle of attack. The quarter chord is
therefore the *aerodynamic centre*. A measured section's centre sits at
:math:`x_{ac}/c = \tfrac14 - m_0/a_0`
(:func:`~engineeringkit.aerodynamics.aerodynamic_center`), where :math:`m_0`
is its residual moment slope. This example reproduces Anderson's
:math:`c_{m,c/4} = -0.053` for the NACA 2412 and finds the aerodynamic
centre from panel-method moments.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.aerodynamics import aerodynamic_center, naca4_airfoil, naca4_moment_coefficient, panel_method

cm = naca4_moment_coefficient(0.02, 0.4)
print(f"NACA 2412 thin-airfoil c_m,c/4 = {cm:.4f}")
assert np.isclose(cm, -0.053, atol=5e-4)

# %%
# Moment about different reference points
# ---------------------------------------
#
# Moving the reference point from the quarter chord to :math:`x` adds
# :math:`c_l(x/c - \tfrac14)`. Only about the aerodynamic centre does the
# moment stay flat as the angle of attack changes.

alpha = np.radians(np.linspace(-4, 8, 7))
results = [panel_method(naca4_airfoil("2412", 121), a) for a in alpha]
cl = np.array([r.cl_pressure for r in results])
cm_quarter = np.array([r.cm_quarter_chord for r in results])
a0 = np.polyfit(alpha, cl, 1)[0]
m0 = np.polyfit(alpha, cm_quarter, 1)[0]
x_ac = float(aerodynamic_center(a0, m0))
print(f"panel method: a0 = {a0:.3f} per rad, m0 = {m0:.4f} per rad, x_ac/c = {x_ac:.4f}")
assert np.isclose(x_ac, 0.25, atol=0.02)

fig, ax = plt.subplots(figsize=(7, 4))
for x_ref in (0.0, 0.25, x_ac, 0.5):
    ax.plot(np.degrees(alpha), cm_quarter + cl * (x_ref - 0.25), "o-", ms=3, label=f"about x/c = {x_ref:.3f}")
ax.axhline(cm, color="k", ls="--", lw=0.8, label="thin airfoil, c/4")
ax.set(xlabel="angle of attack [deg]", ylabel="$c_m$", title="NACA 2412: the moment is flat about the aerodynamic centre")
ax.legend(fontsize=8)
fig.tight_layout()

flat = cm_quarter + cl * (x_ac - 0.25)
assert np.ptp(flat) < 0.1 * np.ptp(cm_quarter - cl * 0.25)

plt.show()

r"""
Lewis's gear-tooth bending formula (1892)
=========================================

Wilfred Lewis, of William Sellers and Co. in Philadelphia, gave the first
rational strength calculation for gear teeth in a talk to the Engineers'
Club of Philadelphia on 15 October 1892. He treated the tooth as a
cantilever with the whole load at its tip and noticed that a parabola
drawn from the load point, tangent to the tooth flanks at the root, is a
beam of uniform strength. The root section where they touch is therefore
the weakest, and the stress there is
:math:`\sigma = W_t/(F m Y)`
(:func:`~engineeringkit.mechanical_engineering.lewis_bending_stress`), with
a form factor :math:`Y` that depends only on the tooth shape and so on the
number of teeth
(:func:`~engineeringkit.mechanical_engineering.lewis_form_factor`). Carl
Barth's factor raises the stress for the dynamic load at speed
(:func:`~engineeringkit.mechanical_engineering.barth_velocity_factor`). The
AGMA bending equation still starts from Lewis's.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.mechanical_engineering import barth_velocity_factor, lewis_bending_stress, lewis_form_factor

# %%
# The parabola of uniform strength
# --------------------------------
#
# A cantilever of width :math:`F` with thickness :math:`t(y) = 2\sqrt{c\,y}`
# at distance :math:`y` below the load has bending stress
# :math:`6Wy/(F t^2) = 3W/(2Fc)` at every section. If the tooth's root
# thickness :math:`t` is reached at depth :math:`l`, Lewis's
# :math:`x = t^2/4l` and :math:`Y = 2x/3m` give the same stress.

W, F, module = 1500.0, 0.02, 0.003
y = np.linspace(1e-4, 0.006, 50)
c = 1.2e-3
thickness = 2 * np.sqrt(c * y)
stress = 6 * W * y / (F * thickness**2)
assert np.allclose(stress, stress[0])
t_root, depth = 0.0055, 0.0068
Y = 2 * (t_root**2 / (4 * depth)) / (3 * module)
print(f"form factor of this tooth: Y = {Y:.3f}")
assert np.isclose(lewis_bending_stress(W, F, module, Y), 6 * W * depth / (F * t_root**2))

# %%
# Sizing a pinion
# ---------------
#
# A 20-tooth, 3 mm module cut pinion transmits 5 kW at 1500 rpm with an
# allowable bending stress of 80 MPa. Find the face width.

teeth, power, rpm, allowable = 20, 5e3, 1500.0, 80e6
pitch_diameter = teeth * module
V = np.pi * pitch_diameter * rpm / 60
Wt = power / V
Kv = barth_velocity_factor(V, "cut")
face = Kv * Wt / (module * lewis_form_factor(teeth) * allowable)
print(f"V = {V:.2f} m/s, Wt = {Wt:.0f} N, Kv = {Kv:.3f}, face width = {face * 1e3:.1f} mm")
assert np.isclose(lewis_bending_stress(Wt, face, module, lewis_form_factor(teeth), Kv), allowable)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
n = np.arange(12, 401)
ax1.semilogx(n, lewis_form_factor(n))
ax1.set(xlabel="number of teeth", ylabel="Lewis form factor Y", title="Bigger gears have stronger teeth")
v = np.linspace(0, 25, 200)
for profile in ("cast", "cut", "hobbed", "ground"):
    ax2.plot(v, barth_velocity_factor(v, profile), label=profile)
ax2.set(xlabel="pitch-line velocity [m/s]", ylabel="Kv", title="Barth's dynamic factor")
ax2.legend()
fig.tight_layout()

plt.show()

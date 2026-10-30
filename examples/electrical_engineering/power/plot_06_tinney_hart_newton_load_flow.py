r"""
Tinney and Hart's Newton-Raphson load flow (1967)
=================================================

Gauss-Seidel load flow needs little memory, but on large networks it
converges slowly or not at all. Newton's method converges quadratically,
yet its Jacobian is as large as the network squared. William Tinney and
Clifford Hart of the Bonneville Power Administration made it practical in
1967 by keeping the Jacobian sparse and ordering the buses so that its
factorization stays sparse too. Newton-Raphson has been the standard power
flow method ever since. In polar form the unknowns are the bus angles and
the load-bus voltage magnitudes; a generator (PV) bus holds its voltage
magnitude and lets its reactive power follow
(:func:`~engineeringkit.electrical_engineering.newton_raphson_load_flow`).
This example solves a five-bus network with a generator bus, shows the
quadratic convergence against Gauss-Seidel
(:func:`~engineeringkit.electrical_engineering.gauss_seidel_load_flow`) and
checks the power balance.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.electrical_engineering import admittance_matrix, gauss_seidel_load_flow, newton_raphson_load_flow

# %%
# The network
# -----------
#
# Bus 0 is the slack, bus 2 a generator holding 1.03 pu while it delivers
# 0.6 pu; buses 1, 3 and 4 are loads.

lines = [
    (0, 1, 0.02 + 0.06j, 0.06j),
    (0, 2, 0.08 + 0.24j, 0.05j),
    (1, 2, 0.06 + 0.18j, 0.04j),
    (1, 3, 0.06 + 0.18j, 0.04j),
    (1, 4, 0.04 + 0.12j, 0.03j),
    (2, 3, 0.01 + 0.03j, 0.02j),
    (3, 4, 0.08 + 0.24j, 0.05j),
]
Y = admittance_matrix(5, lines)
p_load = [0.45, -0.6, 0.4, 0.6]
q_load = [0.15, 0.0, 0.05, 0.1]
r = newton_raphson_load_flow(Y, p_load, q_load, pv_buses={2: 1.03}, tol=1e-12)
for bus, (v, s) in enumerate(zip(r.voltages, r.power)):
    print(f"bus {bus}: |V| = {abs(v):.4f} pu at {np.degrees(np.angle(v)):7.3f} deg, injection {s.real:+.4f} {s.imag:+.4f}j pu")
assert r.converged and np.isclose(abs(r.voltages[2]), 1.03)
assert np.allclose(r.power[1:].real, -np.array(p_load))
assert np.allclose(r.power[[1, 3, 4]].imag, -np.array(q_load)[[0, 2, 3]])
# The slack supplies the net load plus the losses.
assert np.isclose(r.power[0].real, sum(p_load) + r.losses.real)
print(f"losses {r.losses.real:.4f} pu, generator reactive output {r.power[2].imag:.4f} pu")

# %%
# Quadratic convergence
# ---------------------
#
# Fix bus 2's reactive output at the value Newton found and the same case
# becomes all-PQ, which Gauss-Seidel can solve; Newton needs a handful of
# steps for what takes Gauss-Seidel dozens of sweeps.

q_pq = list(q_load)
q_pq[1] = -r.power[2].imag
gs = gauss_seidel_load_flow(Y, p_load, q_pq, tol=1e-10, max_iter=2000)
nr = newton_raphson_load_flow(Y, p_load, q_pq, tol=1e-10)
assert np.allclose(gs.voltages, r.voltages, atol=1e-8) and np.allclose(nr.voltages, r.voltages, atol=1e-10)
print(f"Newton-Raphson: {nr.iterations} steps; Gauss-Seidel: {gs.iterations} sweeps")
print("Newton mismatches:", " ".join(f"{m:.1e}" for m in nr.history))
assert nr.iterations <= 6 < gs.iterations
# Each mismatch is of the order of the square of the previous one, until
# rounding error (about 1e-14) is reached.
assert all(b < 10 * a**2 for a, b in zip(nr.history[1:], nr.history[2:]) if b > 1e-12)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
axes[0].semilogy(nr.history, "o-", label="Newton-Raphson (power mismatch)")
axes[0].semilogy(gs.history, label="Gauss-Seidel (voltage update)")
axes[0].set(xlabel="iteration", ylabel="convergence measure [pu]", title="Convergence")
axes[0].legend(fontsize="small")
r.plot(ax=axes[1])
axes[1].set_title("Bus voltages")
fig.tight_layout()

plt.show()

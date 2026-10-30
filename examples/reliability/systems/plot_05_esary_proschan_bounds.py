r"""
Esary and Proschan's path and cut bounds (1963)
===============================================

Most block diagrams are neither series nor parallel; the classic example is
the bridge. James Esary and Frank Proschan (1963) described any coherent
system by its *minimal path sets*, the smallest groups of components whose
working keeps it working (:func:`~engineeringkit.reliability.minimal_path_sets`),
and its *minimal cut sets*, the smallest groups whose failure stops it
(:func:`~engineeringkit.reliability.minimal_cut_sets`). The system is a
parallel arrangement of its paths, or a series arrangement of its cuts.
Treating the paths or the cuts as independent gives bounds on the exact
reliability (:func:`~engineeringkit.reliability.esary_proschan_bounds`,
:func:`~engineeringkit.reliability.network_reliability`), and the cut
bound is very close for reliable components.
"""

# %%
import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.reliability import (
    esary_proschan_bounds,
    minimal_cut_sets,
    minimal_path_sets,
    network_reliability,
    parallel_reliability,
    series_reliability,
)

# %%
# The bridge network
# ------------------
#
# Five components join source s to sink t; C bridges the two branches.

bridge = {"A": ("s", 1), "B": ("s", 2), "C": (1, 2), "D": (1, "t"), "E": (2, "t")}
paths, cuts = minimal_path_sets(bridge, "s", "t"), minimal_cut_sets(bridge, "s", "t")
print("minimal path sets:", paths)
print("minimal cut sets: ", cuts)

# %%
# Exact reliability and the two bounds
# ------------------------------------
#
# With identical components the exact result is
# :math:`2p^2 + 2p^3 - 5p^4 + 2p^5`. Removing C leaves two series pairs in
# parallel; shorting it leaves two parallel pairs in series. The bridge lies
# between the two.

p = np.linspace(0.01, 0.99, 99)
exact = np.array([network_reliability(bridge, dict.fromkeys(bridge, q), "s", "t") for q in p])
bounds = np.array([esary_proschan_bounds(bridge, dict.fromkeys(bridge, q), "s", "t") for q in p])
without_c = np.array([parallel_reliability(series_reliability(q, q), series_reliability(q, q)) for q in p])
shorted_c = np.array([series_reliability(parallel_reliability(q, q), parallel_reliability(q, q)) for q in p])

fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
axes[0].plot(p, exact, "k", label="exact")
axes[0].plot(p, bounds[:, 0], "--", label="cut-set lower bound")
axes[0].plot(p, bounds[:, 1], ":", label="path-set upper bound")
axes[0].set(xlabel="component reliability", ylabel="system reliability", title="Bridge network")
axes[0].legend(fontsize=8)
axes[1].semilogy(p, exact - bounds[:, 0] + 1e-18, label="exact - lower")
axes[1].semilogy(p, bounds[:, 1] - exact + 1e-18, label="upper - exact")
axes[1].set(xlabel="component reliability", ylabel="error of the bound", ylim=(1e-8, 1))
axes[1].legend(fontsize=8)
fig.tight_layout()

i = np.argmin(np.abs(p - 0.95))
print(f"p = 0.95: exact {exact[i]:.6f}, bounds [{bounds[i, 0]:.6f}, {bounds[i, 1]:.6f}]")

assert len(paths) == 4 and len(cuts) == 4
assert np.allclose(exact, 2 * p**2 + 2 * p**3 - 5 * p**4 + 2 * p**5)
assert np.all(bounds[:, 0] <= exact + 1e-12) and np.all(exact <= bounds[:, 1] + 1e-12)
assert np.all(without_c <= exact + 1e-12) and np.all(exact <= shorted_c + 1e-12)
assert exact[i] - bounds[i, 0] < 1e-4

plt.show()

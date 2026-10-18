r"""
Birnbaum, Esary and Saunders's k-out-of-n systems (1961)
========================================================

Birnbaum, Esary and Saunders (1961) put system reliability on a general
footing: a *structure function* says which combinations of working
components keep the system working. The :math:`k`-out-of-:math:`n` system,
which works if at least :math:`k` of :math:`n` identical components do,
spans the range from parallel (:math:`k = 1`) to series (:math:`k = n`), and
its reliability is a binomial tail
(:func:`~engineeringkit.reliability.k_of_n_reliability`). The example checks
it by enumerating component states and recovers the mean life of each
architecture.
"""

# %%
import itertools

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

from engineeringkit.reliability import exponential_reliability, k_of_n_reliability

# %%
# Enumerate the structure function
# --------------------------------
#
# Four engines, aircraft safe with any :math:`k` working, each with
# reliability 0.95 for the flight.

n, r = 4, 0.95
for k in range(1, n + 1):
    states = itertools.product([0, 1], repeat=n)
    exact = sum(np.prod([r if s else 1 - r for s in state]) for state in states if sum(state) >= k)
    print(f"{k}-out-of-{n}: binomial {k_of_n_reliability(n, k, r):.6f}, enumeration {exact:.6f}")
    assert np.isclose(k_of_n_reliability(n, k, r), exact)

# %%
# Mean life of each architecture
# ------------------------------
#
# With exponential components the :math:`k`-out-of-:math:`n` MTTF is
# :math:`\frac{1}{\lambda}\sum_{i=k}^{n} 1/i`. Time unit: hours.

lam = 1e-3
t = np.linspace(0, 4000, 300)
fig, ax = plt.subplots(figsize=(6.5, 3.8))
for k in range(1, n + 1):
    R = [k_of_n_reliability(n, k, float(x)) for x in exponential_reliability(t, lam)]
    mttf = quad(lambda s, k=k: k_of_n_reliability(n, k, float(exponential_reliability(s, lam))), 0, np.inf)[0]
    exact = sum(1 / i for i in range(k, n + 1)) / lam
    ax.plot(t, R, label=f"{k}-out-of-{n}, MTTF {mttf:.0f} h")
    assert np.isclose(mttf, exact)
ax.set(xlabel="time [h]", ylabel="system reliability")
ax.legend()
fig.tight_layout()

plt.show()

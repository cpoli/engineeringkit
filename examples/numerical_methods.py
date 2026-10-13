"""Compare explicit integrators with the analytic exponential decay solution."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from engineeringkit.integrators import integrate_ode


def main():
    output = Path("gallery")
    output.mkdir(exist_ok=True)
    t = np.linspace(0, 4, 21)
    fig, ax = plt.subplots()
    ax.plot(t, np.exp(-t), "k--", label="analytic")
    for method in ("euler", "rk4", "RK45"):
        result = integrate_ode(lambda time, y: -y, [1.0], t, method)
        print(result.summary())
        ax.plot(result.t, result.history[:, 0], label=method)
    ax.set(xlabel="time [s]", ylabel="normalized state")
    ax.legend()
    fig.savefig(output / "integrators.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()

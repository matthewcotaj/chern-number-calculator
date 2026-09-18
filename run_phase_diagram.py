"""
QWZ topological phase diagram.

    python run_phase_diagram.py

Sweeps the mass u and computes the Chern number at each value, producing the
staircase that jumps at u = -2, 0, 2. Saved to figures/phase_diagram.png.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from topo.models import brillouin_grid, qwz_hamiltonian_grid
from topo.berry import chern_number


def main():
    u_values = np.arange(-4.0, 4.0 + 0.1, 0.1) + 0.005
    n_grid = 60
    KX, KY = brillouin_grid(n_grid)

    chern = np.array([round(chern_number(qwz_hamiltonian_grid(KX, KY, u)))
                      for u in u_values])

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.step(u_values, chern, where="mid", lw=2)
    ax.plot(u_values, chern, "o", ms=3, alpha=0.4)
    for uc in (-2, 0, 2):
        ax.axvline(uc, color="crimson", ls="--", lw=1, alpha=0.7)
    ax.set_xlabel("mass parameter  u")
    ax.set_ylabel("Chern number  C")
    ax.set_title("QWZ topological phase diagram\n"
                 "(dashed red lines = gap-closing transitions at u = -2, 0, 2)")
    ax.set_yticks([-1, 0, 1])
    ax.grid(alpha=0.25)
    fig.tight_layout()
    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/phase_diagram.png", dpi=150)
    print("Saved figures/phase_diagram.png")

    print("\nPhases found:")
    print("  u < -2      : C =", chern[0])
    print("  -2 < u < 0  : C =", chern[np.argmin(np.abs(u_values + 1))])
    print("  0 < u < 2   : C =", chern[np.argmin(np.abs(u_values - 1))])
    print("  u > 2       : C =", chern[-1])


if __name__ == "__main__":
    main()

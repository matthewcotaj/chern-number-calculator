"""
QWZ Chern number and Berry curvature map.

    python run_chern.py

Prints the Chern number for several u (each should be very close to an integer)
and saves a Berry-curvature heatmap over the Brillouin zone to
figures/berry_curvature.png.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from topo.models import brillouin_grid, qwz_hamiltonian_grid
from topo.berry import chern_number, berry_curvature_field


def main():
    u_values = [-3.0, -1.0, 1.0, 3.0]
    n_grid = 80
    KX, KY = brillouin_grid(n_grid)

    print(f"Chern numbers (n_grid = {n_grid}):")
    print("-" * 34)
    for u in u_values:
        C = chern_number(qwz_hamiltonian_grid(KX, KY, u))
        print(f"  u = {u:+.0f}    C = {C:+.4f}   (rounds to {round(C):+d})")
    print()

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    extent = [0, 2 * np.pi, 0, 2 * np.pi]
    area = (2 * np.pi / n_grid) ** 2
    for ax, u in zip(axes.flat, u_values):
        H = qwz_hamiltonian_grid(KX, KY, u)
        F = berry_curvature_field(H)
        C = F.sum() / (2 * np.pi)
        im = ax.imshow((F / area).T, origin="lower", extent=extent,
                       cmap="RdBu_r", aspect="equal")
        ax.set_title(f"u = {u:+.0f}   ->   C = {round(C):+d}")
        ax.set_xlabel(r"$k_x$")
        ax.set_ylabel(r"$k_y$")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Berry curvature")

    fig.suptitle("QWZ Berry curvature over the Brillouin zone", fontsize=13)
    fig.tight_layout()
    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/berry_curvature.png", dpi=150)
    print("Saved figures/berry_curvature.png")


if __name__ == "__main__":
    main()

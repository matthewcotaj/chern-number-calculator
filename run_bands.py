"""
QWZ band structure. Watch the gap open and close as u changes.

    python run_bands.py

Draws the two bands along Gamma -> X -> M -> Gamma for several values of the mass
parameter u, saved to figures/bands.png. The gap closes at u = -2, 0, 2.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from topo.models import qwz_energies, qwz_path


def main():
    kxs, kys, ticks, labels = qwz_path(n_per_segment=200)
    idx = np.arange(len(kxs))
    u_values = [-3.0, -2.0, -1.0, 1.0]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    for ax, u in zip(axes.flat, u_values):
        lower, upper = [], []
        for kx, ky in zip(kxs, kys):
            el, eu = qwz_energies(kx, ky, u)
            lower.append(el)
            upper.append(eu)
        ax.plot(idx, lower, lw=2)
        ax.plot(idx, upper, lw=2)
        ax.axhline(0.0, color="gray", lw=0.6, ls="--")
        ax.set_title(f"u = {u:+.0f}")
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Energy")
        ax.grid(alpha=0.2)

    fig.suptitle("QWZ band structure (gap closes as u approaches -2)", fontsize=12)
    fig.tight_layout()
    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/bands.png", dpi=150)
    print("Saved figures/bands.png")


if __name__ == "__main__":
    main()

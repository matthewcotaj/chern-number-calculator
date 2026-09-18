"""
Haldane model: bands, Chern numbers, and the famous 2D phase diagram.

    python run_haldane.py

Produces two figures:

  figures/haldane_bands.png
    Band structure along Gamma -> K -> M -> K' -> Gamma for three masses. Notice
    that at the transition the gap closes at only ONE Dirac point (K or K'),
    which is the microscopic origin of the topological change.

  figures/haldane_phase_diagram.png
    Haldane's original phase diagram: the Chern number over the (phi, m/t2)
    plane. Two topological lobes (C = +1 and C = -1) separated from the trivial
    region (C = 0) by the curves m/t2 = +/- 3*sqrt(3)*sin(phi).

Also prints Chern numbers for a few representative points.
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from topo.models import (honeycomb_grid, haldane_hamiltonian_grid,
                         haldane_energies, haldane_path)
from topo.berry import chern_number


T2 = 0.15          # next-nearest-neighbor hopping strength
PHI = np.pi / 2    # the time-reversal-breaking phase used for the band plots


def plot_bands():
    kxs, kys, ticks, labels = haldane_path(n_per_segment=200)
    idx = np.arange(len(kxs))
    masses = [0.0, 3 * np.sqrt(3) * T2 * np.sin(PHI), 1.4]  # topo, transition, trivial
    titles = ["m = 0  (topological)",
              "m = m_c  (gap closes at K)",
              "m = 1.4  (trivial)"]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
    for ax, m, title in zip(axes, masses, titles):
        lower, upper = [], []
        for kx, ky in zip(kxs, kys):
            el, eu = haldane_energies(kx, ky, m=m, t2=T2, phi=PHI)
            lower.append(el)
            upper.append(eu)
        ax.plot(idx, lower, lw=2, color="#2f6fb0")
        ax.plot(idx, upper, lw=2, color="#c0562f")
        ax.set_title(title, fontsize=11)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Energy")
    fig.suptitle("Haldane band structure  (t2 = 0.15, phi = pi/2)", fontsize=13)
    fig.tight_layout()
    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/haldane_bands.png", dpi=150)
    print("Saved figures/haldane_bands.png")


def plot_phase_diagram():
    n_phi, n_m = 121, 121
    phis = np.linspace(-np.pi, np.pi, n_phi)
    m_over_t2 = np.linspace(-6.0, 6.0, n_m)

    KX, KY = honeycomb_grid(24)
    C = np.zeros((n_m, n_phi))
    for a, phi in enumerate(phis):
        for b, mr in enumerate(m_over_t2):
            H = haldane_hamiltonian_grid(KX, KY, m=mr * T2, t2=T2, phi=phi)
            C[b, a] = round(chern_number(H))

    fig, ax = plt.subplots(figsize=(9, 6.5))
    im = ax.imshow(C, origin="lower", aspect="auto",
                   extent=[-np.pi, np.pi, -6, 6],
                   cmap="coolwarm", vmin=-1.5, vmax=1.5)

    # Analytic phase boundaries m/t2 = +/- 3 sqrt(3) sin(phi).
    fine = np.linspace(-np.pi, np.pi, 400)
    ax.plot(fine, 3 * np.sqrt(3) * np.sin(fine), "k--", lw=1.2)
    ax.plot(fine, -3 * np.sqrt(3) * np.sin(fine), "k--", lw=1.2)

    ax.set_xlabel(r"phase  $\phi$")
    ax.set_ylabel(r"$m / t_2$")
    ax.set_title("Haldane phase diagram\nChern number over the (phi, m/t2) plane")
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    cbar = fig.colorbar(im, ax=ax, ticks=[-1, 0, 1], fraction=0.046, pad=0.04)
    cbar.set_label("Chern number")
    ax.text(-np.pi / 2, 0, "C = +1", ha="center", fontsize=11, fontweight="bold")
    ax.text(np.pi / 2, 0, "C = -1", ha="center", fontsize=11, fontweight="bold")
    ax.text(0, 4.3, "C = 0", ha="center", fontsize=11, fontweight="bold")
    ax.text(0, -4.6, "C = 0", ha="center", fontsize=11, fontweight="bold")

    fig.tight_layout()
    fig.savefig("figures/haldane_phase_diagram.png", dpi=150)
    print("Saved figures/haldane_phase_diagram.png")


def print_representative():
    KX, KY = honeycomb_grid(72)
    print("\nRepresentative Chern numbers (t2 = 0.15):")
    for phi, m, tag in [(np.pi / 2, 0.0, "phi=+pi/2, m=0  "),
                        (-np.pi / 2, 0.0, "phi=-pi/2, m=0  "),
                        (np.pi / 2, 1.4, "phi=+pi/2, m=1.4")]:
        C = chern_number(haldane_hamiltonian_grid(KX, KY, m=m, t2=0.15, phi=phi))
        print(f"  {tag} ->  C = {round(C):+d}")


def main():
    plot_bands()
    plot_phase_diagram()
    print_representative()


if __name__ == "__main__":
    main()

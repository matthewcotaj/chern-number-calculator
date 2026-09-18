"""
Bulk-boundary correspondence: edge states of the Haldane ribbon.

    python run_edge_states.py

Produces:

  figures/edge_bands.png
    Ribbon band structures side by side. Left: a topological ribbon (bulk C = -1)
    whose gap is bridged by edge bands, colored by which edge they live on. Right:
    a trivial ribbon (C = 0) with a clean, empty gap. The count of gap-crossing
    branches equals 2*|C|, and each branch is confirmed to be edge-localized by
    its color.

  figures/edge_wavefunction.png
    The probability profile |psi|^2 across the width of the ribbon for an in-gap
    state, showing it is pinned to one edge and decays into the insulating bulk.

Also prints, for several parameter choices, the bulk Chern number next to the
number of gap-crossing edge branches, i.e. the correspondence itself.
"""

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from topo.models import honeycomb_grid, haldane_hamiltonian_grid
from topo.berry import chern_number
from topo.ribbon import ribbon_bands, ribbon_hamiltonian, edge_mode_crossings, _ribbon_sites


T2 = 0.15
PHI = np.pi / 2
WIDTH = 24


def _colored_bands(ax, ks, E, pol, title):
    """Draw each band as a line colored by edge polarization (blue/red = edges)."""
    kk = ks / np.pi
    for n in range(E.shape[1]):
        pts = np.array([kk, E[:, n]]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        c = 0.5 * (pol[:-1, n] + pol[1:, n])
        lc = LineCollection(segs, cmap="coolwarm", norm=plt.Normalize(-1, 1))
        lc.set_array(c)
        lc.set_linewidth(1.6)
        ax.add_collection(lc)
    ax.set_xlim(0, 2)
    ax.set_ylim(E.min() * 1.05, E.max() * 1.05)
    ax.set_xlabel(r"momentum $k / \pi$")
    ax.set_title(title)
    ax.grid(alpha=0.2)


def plot_edge_bands():
    ks_t, E_t, pol_t = ribbon_bands(WIDTH, m=0.0, t2=T2, phi=PHI, n_k=241)
    ks_v, E_v, pol_v = ribbon_bands(WIDTH, m=1.4, t2=T2, phi=PHI, n_k=241)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)
    _colored_bands(axes[0], ks_t, E_t, pol_t,
                   "Topological ribbon (bulk C = -1)\nedge bands cross the gap")
    _colored_bands(axes[1], ks_v, E_v, pol_v,
                   "Trivial ribbon (C = 0)\nclean gap, no edge bands")
    axes[0].set_ylabel("energy")

    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(-1, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.04, pad=0.03)
    cbar.set_label("edge localization  (one edge  <->  other edge)")

    fig.suptitle("Bulk-boundary correspondence in the Haldane model", fontsize=13)
    os.makedirs("figures", exist_ok=True)
    fig.savefig("figures/edge_bands.png", dpi=150, bbox_inches="tight")
    print("Saved figures/edge_bands.png")


def plot_edge_wavefunction():
    # Pick a momentum where an edge mode sits inside the gap, then grab the
    # in-gap eigenstate and show its weight per row across the ribbon width.
    k = 2.0
    H = ribbon_hamiltonian(k, WIDTH, m=0.0, t2=T2, phi=PHI)
    evals, evecs = np.linalg.eigh(H)
    _, _, row = _ribbon_sites(WIDTH)

    # The state closest to mid-gap (E = 0 here) is an edge state.
    idx = np.argmin(np.abs(evals))
    psi = evecs[:, idx]
    weight = np.abs(psi) ** 2
    per_row = np.array([weight[row == l].sum() for l in range(WIDTH)])

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(np.arange(WIDTH), per_row, color="#2f6fb0")
    ax.set_xlabel("row across the ribbon width (edge -> bulk -> edge)")
    ax.set_ylabel(r"$|\psi|^2$ on that row")
    ax.set_title(f"An in-gap state at k = {k:.1f} is localized on one edge\n"
                 f"(energy E = {evals[idx]:+.3f}, inside the bulk gap)")
    ax.grid(alpha=0.2, axis="y")
    fig.tight_layout()
    fig.savefig("figures/edge_wavefunction.png", dpi=150)
    print("Saved figures/edge_wavefunction.png")


def print_correspondence():
    KX, KY = honeycomb_grid(60)
    print("\nBulk-boundary correspondence (width = {} rows):".format(WIDTH))
    print("  bulk Chern C   vs   edge branches crossing the gap")
    print("  " + "-" * 48)
    for m in [0.0, 0.4, 1.0, 1.4]:
        C = round(chern_number(haldane_hamiltonian_grid(KX, KY, m=m, t2=T2, phi=PHI)))
        crossings = edge_mode_crossings(WIDTH, m=m, t2=T2, phi=PHI)
        print(f"   m = {m:+.1f}:   C = {C:+d}   ->   {crossings} crossings"
              f"   (2|C| = {2 * abs(C)})")


def main():
    plot_edge_bands()
    plot_edge_wavefunction()
    print_correspondence()


if __name__ == "__main__":
    main()

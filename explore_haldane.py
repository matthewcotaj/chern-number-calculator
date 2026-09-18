"""
Interactive explorer for the Haldane model.

    python explore_haldane.py

Three sliders control the physics directly: m (staggered sublattice mass),
phi (the next-nearest-neighbor phase that breaks time-reversal symmetry), and
t2 (next-nearest-neighbor hopping strength). Dragging them updates the band
structure, the Berry curvature map, your position on the (phi, m/t2) phase
diagram, and the Chern number banner, all live.

The dashed curves on the phase diagram are the exact boundaries
m/t2 = +/- 3*sqrt(3)*sin(phi); cross one and watch the Chern number jump.

Needs a display. If no window appears, use run_haldane.py instead, which
saves the same figures as images.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from topo.models import (honeycomb_grid, haldane_hamiltonian_grid,
                         haldane_energies, haldane_path, haldane_dirac_masses)
from topo.berry import berry_curvature_field


# Appearance
BG = "#f4f4f2"
INK = "#22252a"
LOWER_COLOR = "#2f6fb0"
UPPER_COLOR = "#c0562f"
ACCENT = "#7a3fb0"
TRIVIAL_BG = "#dfe3e8"
TOPO_BG = "#e9dcf5"

plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": "white",
    "axes.edgecolor": "#bcbcc0", "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "font.size": 11, "axes.titlesize": 12,
})

PATH_KX, PATH_KY, PATH_TICKS, PATH_LABELS = haldane_path(n_per_segment=140)
PATH_INDEX = np.arange(len(PATH_KX))

# Precompute the (phi, m/t2) phase-diagram background once, since it doesn't
# depend on the sliders.
PD_PHI = np.linspace(-np.pi, np.pi, 90)
PD_MR = np.linspace(-6.0, 6.0, 90)
_PD_KX, _PD_KY = honeycomb_grid(22)
_T2_REF = 0.15
PD_C = np.zeros((len(PD_MR), len(PD_PHI)))
for _a, _phi in enumerate(PD_PHI):
    for _b, _mr in enumerate(PD_MR):
        _H = haldane_hamiltonian_grid(_PD_KX, _PD_KY, m=_mr * _T2_REF, t2=_T2_REF, phi=_phi)
        PD_C[_b, _a] = round(berry_curvature_field(_H).sum() / (2 * np.pi))


def band_curves(m, t2, phi):
    lower = np.empty_like(PATH_INDEX, dtype=float)
    upper = np.empty_like(PATH_INDEX, dtype=float)
    for n, (kx, ky) in enumerate(zip(PATH_KX, PATH_KY)):
        el, eu = haldane_energies(kx, ky, m=m, t2=t2, phi=phi)
        lower[n], upper[n] = el, eu
    return lower, upper


fig = plt.figure(figsize=(13.5, 8.8))
fig.canvas.manager.set_window_title("Haldane Model Explorer")

gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0],
                      left=0.07, right=0.97, top=0.80, bottom=0.24,
                      hspace=0.45, wspace=0.26)
ax_bands = fig.add_subplot(gs[0, 0])
ax_berry = fig.add_subplot(gs[0, 1])
ax_phase = fig.add_subplot(gs[1, :])

banner = fig.text(0.5, 0.905, "", ha="center", va="center", fontsize=17, fontweight="bold")
banner_box = dict(boxstyle="round,pad=0.6", linewidth=0)
sub_banner = fig.text(0.5, 0.855, "", ha="center", va="center", fontsize=11)

ax_m = fig.add_axes([0.22, 0.125, 0.58, 0.03])
ax_phi = fig.add_axes([0.22, 0.080, 0.58, 0.03])
ax_t2 = fig.add_axes([0.22, 0.035, 0.58, 0.03])

slider_m = Slider(ax_m, "mass  m", -1.6, 1.6, valinit=0.0, valstep=0.05, color=ACCENT)
slider_phi = Slider(ax_phi, "phase  phi", -np.pi, np.pi, valinit=np.pi / 2,
                    valstep=np.pi / 60, color=ACCENT)
slider_t2 = Slider(ax_t2, "hopping  t2", 0.02, 0.30, valinit=0.15, valstep=0.01,
                   color="#6b7280")

# Static phase-diagram background.
ax_phase.imshow(PD_C, origin="lower", aspect="auto", extent=[-np.pi, np.pi, -6, 6],
                cmap="coolwarm", vmin=-1.5, vmax=1.5)
_fine = np.linspace(-np.pi, np.pi, 400)
ax_phase.plot(_fine, 3 * np.sqrt(3) * np.sin(_fine), "k--", lw=1)
ax_phase.plot(_fine, -3 * np.sqrt(3) * np.sin(_fine), "k--", lw=1)
ax_phase.set_xlabel(r"phase  $\phi$")
ax_phase.set_ylabel(r"$m / t_2$")
ax_phase.set_title("Phase diagram (your position marked); dashed = exact boundaries")
ax_phase.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
ax_phase.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
phase_dot, = ax_phase.plot([], [], "o", color="black", ms=11,
                           markerfacecolor="yellow", markeredgewidth=2, zorder=6)


def update(_=None):
    """Redraw everything that depends on the sliders."""
    m = slider_m.val
    phi = slider_phi.val
    t2 = slider_t2.val
    n_grid = 45

    H = haldane_hamiltonian_grid(*honeycomb_grid(n_grid), m=m, t2=t2, phi=phi)
    field = berry_curvature_field(H)
    chern = int(round(field.sum() / (2 * np.pi)))

    # Bands
    lower, upper = band_curves(m, t2, phi)
    ax_bands.clear()
    ax_bands.plot(PATH_INDEX, lower, lw=2.2, color=LOWER_COLOR, label="lower band")
    ax_bands.plot(PATH_INDEX, upper, lw=2.2, color=UPPER_COLOR, label="upper band")
    ax_bands.fill_between(PATH_INDEX, lower, upper, color="#00000008")
    ax_bands.set_xticks(PATH_TICKS)
    ax_bands.set_xticklabels(PATH_LABELS)
    ax_bands.set_ylabel("energy")
    ax_bands.set_title("Band structure")
    ax_bands.legend(loc="upper right", fontsize=9, framealpha=0.9)
    mK, mKp = haldane_dirac_masses(m=m, t2=t2, phi=phi)
    ax_bands.text(0.03, 0.04,
                  f"Dirac masses:  m_K = {mK:+.2f},  m_K' = {mKp:+.2f}",
                  transform=ax_bands.transAxes, fontsize=8.5, color="#555")

    # Berry curvature (plotted over the reduced-coordinate cell that tiles the BZ)
    ax_berry.clear()
    vmax = max(np.abs(field).max(), 1e-9)
    ax_berry.imshow(field.T, origin="lower", extent=[0, 1, 0, 1],
                    cmap="RdBu_r", aspect="equal", vmin=-vmax, vmax=vmax)
    ax_berry.set_xlabel(r"$k_1$ (reduced)")
    ax_berry.set_ylabel(r"$k_2$ (reduced)")
    ax_berry.set_title("Berry curvature over the Brillouin zone")

    # Phase-diagram marker
    phase_dot.set_data([phi], [m / t2])

    # Banner
    if chern == 0:
        banner.set_text(f"TRIVIAL     C = {chern:+d}")
        banner_box["facecolor"] = TRIVIAL_BG
        sub_banner.set_text("ordinary insulator, no protected edge states")
    else:
        banner.set_text(f"TOPOLOGICAL     C = {chern:+d}")
        banner_box["facecolor"] = TOPO_BG
        sub_banner.set_text(f"Chern insulator: |C| = {abs(chern)} chiral edge mode(s)")
    banner.set_bbox(banner_box)

    fig.canvas.draw_idle()


_p = berry_curvature_field(haldane_hamiltonian_grid(*honeycomb_grid(45)))
_vmax = max(np.abs(_p).max(), 1e-9)
_im0 = ax_berry.imshow(np.zeros((2, 2)), origin="lower", extent=[0, 1, 0, 1],
                       cmap="RdBu_r", vmin=-_vmax, vmax=_vmax)
cbar = fig.colorbar(_im0, ax=ax_berry, fraction=0.046, pad=0.04)
cbar.set_label("curvature")

slider_m.on_changed(update)
slider_phi.on_changed(update)
slider_t2.on_changed(update)
update()

if __name__ == "__main__":
    plt.show()

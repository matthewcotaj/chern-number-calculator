"""
Interactive explorer for the QWZ topological model.

    python explore.py

Drag the mass slider and the bands, Berry curvature map, and Chern number all
update live, along with a marker on the phase diagram below showing where you
are on the u axis.

This needs a display -- if no window opens (headless/remote shell), use
run_bands.py / run_chern.py / run_phase_diagram.py instead, which save the
same figures as images.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from topo.models import qwz_energies, qwz_path, brillouin_grid, qwz_hamiltonian_grid
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
    "figure.facecolor": BG,
    "axes.facecolor": "white",
    "axes.edgecolor": "#bcbcc0",
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "font.size": 11,
    "axes.titlesize": 12,
})


# Fixed pieces that don't depend on the sliders, so compute them once.
PATH_KX, PATH_KY, PATH_TICKS, PATH_LABELS = qwz_path(n_per_segment=160)
PATH_INDEX = np.arange(len(PATH_KX))

# The phase-diagram staircase is independent of the slider too.
PHASE_U = np.arange(-4.0, 4.0 + 1e-9, 0.1) + 0.005
_PHASE_KX, _PHASE_KY = brillouin_grid(40)
PHASE_C = np.array([
    round(berry_curvature_field(qwz_hamiltonian_grid(_PHASE_KX, _PHASE_KY, u)).sum()
          / (2 * np.pi))
    for u in PHASE_U
])


def band_curves(u):
    lower = np.empty_like(PATH_INDEX, dtype=float)
    upper = np.empty_like(PATH_INDEX, dtype=float)
    for n, (kx, ky) in enumerate(zip(PATH_KX, PATH_KY)):
        el, eu = qwz_energies(kx, ky, u)
        lower[n], upper[n] = el, eu
    return lower, upper


fig = plt.figure(figsize=(13.5, 8.6))
fig.canvas.manager.set_window_title("QWZ Topological Explorer")

gs = fig.add_gridspec(
    2, 2,
    height_ratios=[1.0, 0.85],
    left=0.08, right=0.96, top=0.80, bottom=0.20,
    hspace=0.42, wspace=0.24,
)
ax_bands = fig.add_subplot(gs[0, 0])
ax_berry = fig.add_subplot(gs[0, 1])
ax_phase = fig.add_subplot(gs[1, :])

# Status banner across the top.
banner = fig.text(0.5, 0.90, "", ha="center", va="center",
                  fontsize=17, fontweight="bold")
banner_box = dict(boxstyle="round,pad=0.6", linewidth=0)
sub_banner = fig.text(0.5, 0.855, "", ha="center", va="center", fontsize=11)

# Slider axes.
ax_u = fig.add_axes([0.20, 0.085, 0.60, 0.03])
ax_grid = fig.add_axes([0.20, 0.035, 0.60, 0.03])

slider_u = Slider(ax_u, "mass  u", -4.0, 4.0, valinit=1.0, valstep=0.1,
                  color=ACCENT)
slider_grid = Slider(ax_grid, "grid  N", 20, 100, valinit=45, valstep=5,
                     color="#6b7280")


# Static parts of the phase diagram.
ax_phase.step(PHASE_U, PHASE_C, where="mid", lw=2.2, color=ACCENT)
for uc in (-2, 0, 2):
    ax_phase.axvline(uc, color="#c94b4b", ls="--", lw=1, alpha=0.6)
ax_phase.set_yticks([-1, 0, 1])
ax_phase.set_xlabel("mass parameter  u")
ax_phase.set_ylabel("Chern number  C")
ax_phase.set_title("Phase diagram  (dashed = gap-closing transitions)")
ax_phase.grid(alpha=0.25)
phase_marker = ax_phase.axvline(slider_u.valinit, color=INK, lw=2)
phase_dot, = ax_phase.plot([], [], "o", color=INK, ms=9, zorder=5)


def update(_=None):
    """Redraw everything that depends on the sliders."""
    u = slider_u.val
    n_grid = int(slider_grid.val)

    # Berry curvature + Chern number, computed once and reused below.
    field = berry_curvature_field(qwz_hamiltonian_grid(*brillouin_grid(n_grid), u))
    chern = int(round(field.sum() / (2 * np.pi)))
    area = (2 * np.pi / n_grid) ** 2
    density = field / area

    # Bands
    lower, upper = band_curves(u)
    ax_bands.clear()
    ax_bands.plot(PATH_INDEX, lower, lw=2.2, color=LOWER_COLOR, label="lower band")
    ax_bands.plot(PATH_INDEX, upper, lw=2.2, color=UPPER_COLOR, label="upper band")
    ax_bands.axhline(0, color="gray", lw=0.6, ls="--")
    ax_bands.fill_between(PATH_INDEX, lower, upper, color="#00000008")
    ax_bands.set_xticks(PATH_TICKS)
    ax_bands.set_xticklabels(PATH_LABELS)
    ax_bands.set_ylabel("energy")
    ax_bands.set_title("Band structure")
    ax_bands.set_ylim(-4.4, 4.4)
    ax_bands.legend(loc="upper right", fontsize=9, framealpha=0.9)
    gap = float(np.min(upper - lower))
    ax_bands.text(0.03, 0.04, f"min gap = {gap:.2f}", transform=ax_bands.transAxes,
                  fontsize=9, color="#555")

    # Berry curvature heatmap
    ax_berry.clear()
    vmax = max(np.abs(density).max(), 1e-9)
    im = ax_berry.imshow(density.T, origin="lower",
                         extent=[0, 2 * np.pi, 0, 2 * np.pi],
                         cmap="RdBu_r", aspect="equal", vmin=-vmax, vmax=vmax)
    ax_berry.set_xlabel(r"$k_x$")
    ax_berry.set_ylabel(r"$k_y$")
    ax_berry.set_title("Berry curvature over the Brillouin zone")

    # Phase-diagram marker
    phase_marker.set_xdata([u, u])
    phase_dot.set_data([u], [chern])

    # Status banner
    if chern == 0:
        banner.set_text(f"TRIVIAL     C = {chern:+d}")
        banner_box["facecolor"] = TRIVIAL_BG
        sub_banner.set_text("insulator with no protected edge states")
    else:
        banner.set_text(f"TOPOLOGICAL     C = {chern:+d}")
        banner_box["facecolor"] = TOPO_BG
        sub_banner.set_text(f"|C| = {abs(chern)} protected chiral edge state(s)")
    banner.set_bbox(banner_box)

    fig.canvas.draw_idle()


_probe = berry_curvature_field(
    qwz_hamiltonian_grid(*brillouin_grid(int(slider_grid.valinit)), slider_u.valinit))
_area = (2 * np.pi / int(slider_grid.valinit)) ** 2
_vmax = max(np.abs(_probe / _area).max(), 1e-9)
_im0 = ax_berry.imshow(np.zeros((2, 2)), origin="lower",
                       extent=[0, 2 * np.pi, 0, 2 * np.pi],
                       cmap="RdBu_r", vmin=-_vmax, vmax=_vmax)
cbar = fig.colorbar(_im0, ax=ax_berry, fraction=0.046, pad=0.04)
cbar.set_label("curvature")

slider_u.on_changed(update)
slider_grid.on_changed(update)
update()

if __name__ == "__main__":
    plt.show()

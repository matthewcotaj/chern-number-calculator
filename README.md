# Topological Band Structure & Chern Number Calculator

A from-scratch implementation of the Chern number, the integer invariant that classifies topological insulators. It's built on two lattice models sharing one topology engine:

- **Qi-Wu-Zhang (QWZ)** — the simplest Chern insulator, on a square lattice
- **Haldane model** — the honeycomb-lattice model that founded the field

The core idea: tuning a continuous parameter (a mass term, or a hopping phase) produces a quantized integer response, the Chern number, and that integer can only change when the band gap closes. Every result here is self-verifying since the Chern number has to come out as an exact integer, so if the code is wrong you'll know immediately.

![QWZ phase diagram](figures/qwz_phase_diagram.png) ![Haldane phase diagram](figures/haldane_phase_diagram.png)

## What you need

- Python 3.9+, plus `numpy` and `matplotlib`
- Linux, macOS, or Windows (WSL works fine)

Background-wise: linear algebra (NumPy handles the diagonalization, but you should know what an eigenvalue is), multivariable calculus (gradient, curl, surface integrals conceptually), and basic quantum mechanics (Hermitian operators, energy eigenstates). No condensed matter background assumed — I explain the physics in the docstrings and in section 5 below.

If you've never used a terminal before, read `GETTING_STARTED.md` first.

## Setup

```bash
cd chern-calculator
python -m venv .venv                 # optional
source .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Or skip the venv and just run `pip install numpy matplotlib`.

## Running it

Run everything from the project root (use `python3` if `python` doesn't work on your machine).

```bash
python test_models.py            # sanity check — should print all PASS

# QWZ (square lattice)
python run_bands.py              # band structure -> figures/bands.png
python run_chern.py              # Chern number + Berry curvature map
python run_phase_diagram.py      # C vs u staircase
python explore.py                # interactive mass slider

# Haldane (honeycomb lattice)
python run_haldane.py            # bands + 2D phase diagram
python explore_haldane.py        # interactive sliders (m, phi, t2)
python run_edge_states.py        # edge modes / bulk-boundary correspondence
```

Figures land in `figures/`.

## Project layout

```
chern-calculator/
├── README.md
├── GETTING_STARTED.md
├── requirements.txt
│
├── topo/
│   ├── __init__.py
│   ├── models.py       # QWZ and Haldane Hamiltonians
│   ├── berry.py         # model-agnostic Chern number engine (FHS method)
│   └── ribbon.py        # Haldane ribbon geometry for edge states
│
├── run_bands.py
├── run_chern.py
├── run_phase_diagram.py
├── explore.py
│
├── run_haldane.py
├── explore_haldane.py
├── run_edge_states.py
│
└── test_models.py
```

Best order to read the code: `topo/models.py` first (the physics), then `topo/berry.py` (the topology). `berry.py` doesn't import either model directly — it just takes a Hamiltonian evaluated on a grid and returns a Chern number. Both models plug into the same engine, which is the part I'm most happy with.

## Suggested order to go through this

1. Run `test_models.py`, make sure everything passes.
2. Read the QWZ part of `models.py`, run `run_bands.py`, figure out why the gap closes at u = -2, 0, 2.
3. Read `berry.py` (there's a docstring about the gauge problem that's worth reading closely), run `run_chern.py`, check the curvature heatmap.
4. Run `run_phase_diagram.py` — this is the main QWZ result.
5. Play with `explore.py`, drag the slider across u = -2, 0, 2 and watch the Chern number jump.
6. Read the Haldane section of `models.py`, run `run_haldane.py`. Notice the gap only closes at one Dirac point at the transition, not both.
7. Play with `explore_haldane.py` across the phase boundaries.
8. Run `run_edge_states.py` — this is the physical payoff. A topological ribbon has conducting edge channels bridging the gap, and the number of them equals the bulk Chern number. Compare to the trivial phase, where the gap stays clean.

## The physics/math

Each Bloch Hamiltonian here is a 2x2 Hermitian matrix, H(k) = eps(k) I + d(k)·sigma. The eps·I term just shifts energies, so all the topology lives in the vector field d(k).

- **QWZ:** d = (sin kx, sin ky, u + cos kx + cos ky). Gap closes at u = -2, 0, 2.
- **Haldane:** nearest-neighbor hopping t1 gives graphene's Dirac cones. A complex next-nearest-neighbor hopping t2·e^(i·phi) breaks time-reversal symmetry, and a staggered mass m breaks inversion symmetry. The gaps at the two Dirac points are m_K = m − 3√3·t2·sin(phi) and m_K' = m + 3√3·t2·sin(phi). Opposite signs → topological, same sign → trivial. Phase boundary: m = ±3√3·t2·sin(phi).

Each eigenstate has a Berry connection A(k) = i⟨u|∇_k|u⟩ (a gradient), a Berry curvature Ω = ∇×A (a curl), and a Chern number C = (1/2π)∫Ω d²k over the Brillouin zone (a surface integral). C comes out as an integer, always.

You can't just differentiate the eigenstates directly, though — a numerical eigensolver returns eigenvectors with arbitrary phases, which wrecks any naive gradient. The Fukui-Hatsugai-Suzuki method gets around this by rewriting the integral in terms of overlaps between neighboring k-points, so the phases cancel and the result is automatically an integer. That's what `berry.py` implements.

The Chern number is a bulk quantity, but it determines edge behavior: a finite ribbon of a Chern insulator has exactly |C| conducting edge channels per edge, one-way, with direction set by the sign of C. `topo/ribbon.py` builds the Haldane model on a finite strip, and `run_edge_states.py` checks that the number of gap-crossing edge bands is 2|C| (one set per edge), and confirms each one is actually localized at an edge.

One thing I had to get right: the honeycomb Brillouin zone is a hexagon, not a square like QWZ's, so the Haldane calculation samples the actual hexagonal reciprocal lattice (`honeycomb_grid`). Sampling a square grid instead gives a wrong Chern number even with a correct Hamiltonian — this tripped me up early on.

## Possible extensions

- Convergence study: Chern number error vs. grid size
- Anomalous Hall conductivity, sigma_xy = C·e²/h — ties the integer to something measurable
- Add the BHZ model (a Z2 topological insulator) to the same engine, moving from Chern to time-reversal-invariant topology
- Add disorder to the ribbon and check that edge channels survive while bulk states localize
- An interactive edge-state visualizer

## References

- Qi, Wu, Zhang, *Phys. Rev. B* 74, 085308 (2006) — the QWZ model
- Haldane, *Phys. Rev. Lett.* 61, 2015 (1988) — the honeycomb model
- Fukui, Hatsugai, Suzuki, *J. Phys. Soc. Jpn.* 74, 1674 (2005) — the discrete Chern number method used in `berry.py`
- Bernevig & Hughes, *Topological Insulators and Topological Superconductors* (Princeton, 2013)

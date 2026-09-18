"""
The Haldane model on a ribbon: where the bulk-boundary correspondence shows up.

The idea
--------
The Chern number is a property of the infinite 2D crystal (the "bulk"). Its
physical punchline is about *edges*: a strip of a Chern insulator must carry
exactly |C| one-way ("chiral") conducting channels along each edge, and their
direction is set by the sign of C. Conducting edges wrapped around an insulating
interior is the whole reason these materials are interesting.

To see it, we cut the 2D lattice into a ribbon: infinite (periodic) along one
lattice direction, but only a finite number of rows W wide in the other. Momentum
k along the periodic direction is still a good quantum number, so we get a band
structure E_n(k). Inside the bulk energy gap, extra bands appear that cross the
gap; each such crossing branch is an edge channel, and we confirm it is real by
checking that its wavefunction is localized on an edge row.

Construction
------------
We build the ribbon's Bloch Hamiltonian directly from the tight-binding rules,
reusing the same geometry as the bulk model (so signs stay consistent):

  * on-site energy  +m on sublattice A, -m on sublattice B,
  * nearest-neighbor hopping t1 between A and B (bond length 1),
  * next-nearest-neighbor hopping t2 with Haldane's phase: +phi or -phi depending
    on the hop direction and sublattice (this is the time-reversal-breaking term).

The ribbon unit cell holds 2*W sites (an A and a B in each of the W rows), so
H(k) is a (2W x 2W) matrix. We diagonalize it at each k.
"""

import numpy as np

from topo.models import DELTA, NU

_SQRT3 = np.sqrt(3.0)

# Lattice vectors for the ribbon. A1 is the periodic direction (momentum k runs
# along it); A2 is the finite "stacking" direction (the ribbon is W cells of A2
# wide). Both are next-nearest-neighbor (same-sublattice) lattice vectors, so
# they are consistent with the bulk geometry in models.py.
A1 = np.array([1.5, -_SQRT3 / 2])   # periodic
A2 = np.array([1.5, _SQRT3 / 2])    # finite
DELTA_INTRA = DELTA[0]              # the A -> B bond taken inside one cell

_TOL = 1e-6


def _ribbon_sites(width):
    """
    Positions and sublattice labels of the 2*W sites in one ribbon unit cell.

    Returns
    -------
    pos : ndarray, shape (2W, 2)      2D position of each site
    sub : ndarray, shape (2W,)        +1 on sublattice A, -1 on sublattice B
    row : ndarray, shape (2W,)        row index l = 0..W-1 (which edge a site is near)
    Index convention: site 2*l is the A of row l, site 2*l+1 is the B of row l.
    """
    pos = np.zeros((2 * width, 2))
    sub = np.zeros(2 * width)
    row = np.zeros(2 * width, dtype=int)
    for l in range(width):
        a_pos = l * A2
        pos[2 * l] = a_pos
        pos[2 * l + 1] = a_pos + DELTA_INTRA
        sub[2 * l] = +1
        sub[2 * l + 1] = -1
        row[2 * l] = l
        row[2 * l + 1] = l
    return pos, sub, row


def ribbon_hamiltonian(k, width, m=0.0, t1=1.0, t2=0.15, phi=np.pi / 2):
    """
    The (2W x 2W) Bloch Hamiltonian of the Haldane ribbon at momentum k.

    k is dimensionless, in [0, 2*pi), measured along the periodic direction A1.
    """
    pos, sub, _ = _ribbon_sites(width)
    n_sites = 2 * width
    H = np.zeros((n_sites, n_sites), dtype=complex)

    # On-site staggered mass.
    for i in range(n_sites):
        H[i, i] += m * sub[i]

    # Hoppings: loop over site pairs and the neighboring cells n_off along A1.
    for i in range(n_sites):
        for n_off in (-1, 0, 1):
            bloch = np.exp(1j * k * n_off)
            for j in range(n_sites):
                disp = pos[i] - (pos[j] + n_off * A1)   # vector from j-image to i
                r = np.hypot(disp[0], disp[1])

                if abs(r - 1.0) < _TOL and sub[i] != sub[j]:
                    # Nearest-neighbor hop (A <-> B).
                    H[i, j] += t1 * bloch

                elif abs(r - _SQRT3) < _TOL and sub[i] == sub[j]:
                    # Next-nearest-neighbor hop with Haldane phase.
                    is_positive = any(np.hypot(*(disp - nu)) < _TOL for nu in NU)
                    dir_sign = 1.0 if is_positive else -1.0
                    sub_sign = 1.0 if sub[i] > 0 else -1.0
                    amp = t2 * np.exp(1j * dir_sign * sub_sign * phi)
                    H[i, j] += amp * bloch

    return H


def ribbon_bands(width, m=0.0, t1=1.0, t2=0.15, phi=np.pi / 2, n_k=241):
    """
    Band structure of the ribbon plus an edge-localization label for each state.

    Returns
    -------
    ks : ndarray, shape (n_k,)                  momentum samples in [0, 2pi)
    energies : ndarray, shape (n_k, 2W)         sorted eigenvalues at each k
    edge_pol : ndarray, shape (n_k, 2W)         edge polarization of each state,
        in [-1, 1]: near +1 if localized on the first rows (one edge), near -1 on
        the last rows (other edge), near 0 for bulk states.
    """
    ks = np.linspace(0.0, 2 * np.pi, n_k)
    _, _, row = _ribbon_sites(width)
    n_states = 2 * width
    energies = np.empty((n_k, n_states))
    edge_pol = np.empty((n_k, n_states))

    # Rows counted as "edge" = the outer quarter on each side.
    edge_cut = max(1, width // 4)
    bottom = row < edge_cut
    top = row >= (width - edge_cut)

    for a, k in enumerate(ks):
        H = ribbon_hamiltonian(k, width, m, t1, t2, phi)
        evals, evecs = np.linalg.eigh(H)
        energies[a] = evals
        weight = np.abs(evecs) ** 2                    # (site, state)
        w_bottom = weight[bottom, :].sum(axis=0)
        w_top = weight[top, :].sum(axis=0)
        edge_pol[a] = w_bottom - w_top

    return ks, energies, edge_pol


def _indirect_gap(m, t1, t2, phi, n=120):
    """
    Return (valence_max, conduction_min) of the 2D bulk bands: the edges of the
    indirect gap. Any energy strictly between them has no bulk states at any k,
    which makes it a safe reference for counting edge modes.
    """
    from topo.models import honeycomb_grid, haldane_energies
    KX, KY = honeycomb_grid(n)
    lower, upper = haldane_energies(KX, KY, m=m, t1=t1, t2=t2, phi=phi)
    return lower.max(), upper.min()


def edge_mode_crossings(width, m=0.0, t1=1.0, t2=0.15, phi=np.pi / 2,
                        n_k=481, e_ref=None):
    """
    Count how many bands cross a reference energy inside the bulk gap.

    For a Chern insulator this equals 2*|C| (|C| channels on each of the two
    edges). We count sign changes of (E_sorted[:, k] - e_ref) as k advances;
    bulk bands never enter the gap, so only edge channels contribute.

    If e_ref is None it is placed automatically off-center inside the indirect
    gap. Off-center matters: at a symmetric spectrum (e.g. phi = pi/2, m = 0) the
    two edges cross at exactly the same energy, and a reference sitting right on
    that energy would miss them.
    """
    if e_ref is None:
        v_max, c_min = _indirect_gap(m, t1, t2, phi)
        if c_min <= v_max:
            e_ref = 0.0                      # no clean gap; fall back
        else:
            e_ref = v_max + 0.6 * (c_min - v_max)

    ks = np.linspace(0.0, 2 * np.pi, n_k, endpoint=True)
    crossings = 0
    prev = None
    for k in ks:
        evals = np.linalg.eigvalsh(ribbon_hamiltonian(k, width, m, t1, t2, phi))
        s = np.sign(evals - e_ref)
        if prev is not None:
            crossings += int(np.count_nonzero(s != prev))
        prev = s
    return crossings

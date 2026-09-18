"""
The two lattice models: Qi-Wu-Zhang (QWZ) and Haldane.

Common structure
----------------
Both are two-band models, so every Bloch Hamiltonian is a 2x2 Hermitian matrix.
Any such matrix can be written as

    H(k) = eps(k) * I  +  d_x(k) * sigma_x  +  d_y(k) * sigma_y  +  d_z(k) * sigma_z

The eps(k) * I piece shifts both bands equally. It changes the energies but not
the eigenvectors, so it has zero effect on the topology (the Chern number). The
vector field d(k) = (d_x, d_y, d_z) is where all the topology lives.

To keep the Chern-number engine fast, each model provides a *grid builder* that
evaluates H on a whole mesh of k-points at once and returns an array of shape
(N, N, 2, 2). NumPy can then diagonalize the entire stack in a single call.
"""

import numpy as np

_SQRT3 = np.sqrt(3.0)

# ---------------------------------------------------------------------------
# Pauli matrices (used by the scalar convenience Hamiltonians below)
# ---------------------------------------------------------------------------
SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
SIGMA_Z = np.array([[1, 0], [0, -1]], dtype=complex)
IDENTITY = np.eye(2, dtype=complex)


def _assemble(eps, dx, dy, dz):
    """
    Build a stack of 2x2 Hamiltonians from component arrays.

    eps, dx, dy, dz can be scalars or arrays of the same shape (say (N, N)).
    Returns an array of shape (*shape, 2, 2):

        [[ eps + dz ,  dx - i dy ],
         [ dx + i dy,  eps - dz  ]]
    """
    eps = np.asarray(eps, dtype=float)
    dx = np.asarray(dx, dtype=float)
    dy = np.asarray(dy, dtype=float)
    dz = np.asarray(dz, dtype=float)

    H = np.empty(eps.shape + (2, 2), dtype=complex)
    H[..., 0, 0] = eps + dz
    H[..., 0, 1] = dx - 1j * dy
    H[..., 1, 0] = dx + 1j * dy
    H[..., 1, 1] = eps - dz
    return H


def brillouin_grid(n_grid):
    """
    A square mesh of k-points covering the square-lattice Brillouin zone
    [0, 2pi) x [0, 2pi). Used by the QWZ model.

    Returns two arrays KX, KY of shape (n_grid, n_grid) with 'ij' indexing, so
    KX[i, j], KY[i, j] is the (kx, ky) of grid point (i, j). The mesh is periodic
    (a torus), which is what the Chern engine needs.
    """
    k = np.linspace(0.0, 2 * np.pi, n_grid, endpoint=False)
    return np.meshgrid(k, k, indexing="ij")


# Reciprocal lattice vectors of the honeycomb lattice (bond length = 1). The
# honeycomb Brillouin zone is a hexagon, NOT the square [0, 2pi)^2, so we must
# sample it using these vectors or the Chern number comes out wrong.
HONEYCOMB_B1 = np.array([2 * np.pi / 3, 2 * np.pi / _SQRT3])
HONEYCOMB_B2 = np.array([2 * np.pi / 3, -2 * np.pi / _SQRT3])


def honeycomb_grid(n_grid):
    """
    A mesh covering exactly one honeycomb Brillouin zone.

    We lay a uniform grid over the *reduced* coordinates (k1, k2) in [0, 1) x
    [0, 1) and map each point to Cartesian momentum via

        k = k1 * b1 + k2 * b2

    with b1, b2 the reciprocal lattice vectors. This tiles one full BZ and is
    periodic on a torus (rolling off one edge lands on the opposite edge), which
    is exactly what the Fukui-Hatsugai-Suzuki method requires. Sampling the wrong
    region -- e.g. a plain square in (kx, ky) -- gives a wrong, non-unit Chern
    number even though the Hamiltonian is correct.

    Returns KX, KY of shape (n_grid, n_grid).
    """
    s = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    K1, K2 = np.meshgrid(s, s, indexing="ij")
    KX = K1 * HONEYCOMB_B1[0] + K2 * HONEYCOMB_B2[0]
    KY = K1 * HONEYCOMB_B1[1] + K2 * HONEYCOMB_B2[1]
    return KX, KY


# ===========================================================================
# Model 1: Qi-Wu-Zhang (square lattice)
# ===========================================================================
#
#   d_x = sin(kx)
#   d_y = sin(ky)
#   d_z = u + cos(kx) + cos(ky)      (u is the mass parameter)
#   eps = 0
#
# Topological transitions (gap closes) at u = -2, 0, 2.

def qwz_components(kx, ky, u):
    """Return (eps, d_x, d_y, d_z) for the QWZ model. Inputs may be arrays."""
    eps = np.zeros_like(np.asarray(kx, dtype=float))
    dx = np.sin(kx)
    dy = np.sin(ky)
    dz = u + np.cos(kx) + np.cos(ky)
    return eps, dx, dy, dz


def qwz_hamiltonian(kx, ky, u):
    """The 2x2 QWZ Bloch Hamiltonian at a single k-point."""
    _, dx, dy, dz = qwz_components(kx, ky, u)
    return dx * SIGMA_X + dy * SIGMA_Y + dz * SIGMA_Z


def qwz_hamiltonian_grid(KX, KY, u):
    """QWZ Hamiltonian evaluated on a whole grid: shape (N, N, 2, 2)."""
    eps, dx, dy, dz = qwz_components(KX, KY, u)
    return _assemble(eps, dx, dy, dz)


def qwz_energies(kx, ky, u):
    """Analytic bands E_pm = +/- |d| for the QWZ model."""
    _, dx, dy, dz = qwz_components(kx, ky, u)
    mag = np.sqrt(dx**2 + dy**2 + dz**2)
    return -mag, mag


# ===========================================================================
# Model 2: Haldane (honeycomb lattice)
# ===========================================================================
#
# The honeycomb lattice has two sublattices, A and B. Three ingredients:
#   t1  : ordinary nearest-neighbor hopping (A <-> B). Gives the graphene cones.
#   t2, phi : next-nearest-neighbor hopping (A<->A and B<->B) with a complex
#             phase +/- phi. The phase breaks time-reversal symmetry and is what
#             makes the model topological.
#   m   : a staggered on-site energy (+m on A, -m on B). Breaks inversion
#         symmetry and competes with the t2 term.
#
# Geometry (bond length = 1):
#   delta_i : the three vectors from an A site to its B neighbors.
#   nu_i    : the three next-nearest-neighbor vectors (same sublattice).
#
# Dirac points sit at K and K'. The gap at each is set by a "Dirac mass":
#   m_K  = m - 3*sqrt(3)*t2*sin(phi)
#   m_Kp = m + 3*sqrt(3)*t2*sin(phi)
# When these two masses have opposite signs the band is topological (C = +/-1);
# when they share a sign it is trivial (C = 0). So the phase boundaries are
#   m = +/- 3*sqrt(3)*t2*sin(phi).

# Nearest-neighbor vectors (A -> B), each of length 1.
DELTA = np.array([
    [0.5, _SQRT3 / 2.0],
    [0.5, -_SQRT3 / 2.0],
    [-1.0, 0.0],
])

# Next-nearest-neighbor vectors (same sublattice), each of length sqrt(3).
# Taken as delta_i - delta_j so the three are related by 120-degree rotations.
NU = np.array([
    [0.0, _SQRT3],       # delta_1 - delta_2
    [1.5, -_SQRT3 / 2],  # delta_2 - delta_3
    [-1.5, -_SQRT3 / 2],  # delta_3 - delta_1
])

# High-symmetry points of the honeycomb Brillouin zone, for band plots.
K_POINT = np.array([2 * np.pi / 3, 2 * np.pi / (3 * _SQRT3)])
KP_POINT = np.array([2 * np.pi / 3, -2 * np.pi / (3 * _SQRT3)])
M_POINT = np.array([2 * np.pi / 3, 0.0])


def haldane_components(kx, ky, m=0.0, t1=1.0, t2=0.15, phi=np.pi / 2):
    """
    Return (eps, d_x, d_y, d_z) for the Haldane model. Inputs may be arrays.

    d_x, d_y come from the nearest-neighbor (t1) term; d_z carries the mass m and
    the time-reversal-breaking t2*sin(phi) term; eps is the pure-identity part
    from t2*cos(phi) (it shifts energies but not topology).
    """
    kx = np.asarray(kx, dtype=float)
    ky = np.asarray(ky, dtype=float)

    # Nearest-neighbor sums: sum over the three delta vectors.
    cos_sum = np.zeros_like(kx)
    sin_sum = np.zeros_like(kx)
    for dx_, dy_ in DELTA:
        phase = kx * dx_ + ky * dy_
        cos_sum += np.cos(phase)
        sin_sum += np.sin(phase)

    dx = t1 * cos_sum
    dy = t1 * sin_sum

    # Next-nearest-neighbor sums: sum over the three nu vectors.
    nnn_cos = np.zeros_like(kx)
    nnn_sin = np.zeros_like(kx)
    for nx_, ny_ in NU:
        phase = kx * nx_ + ky * ny_
        nnn_cos += np.cos(phase)
        nnn_sin += np.sin(phase)

    eps = 2.0 * t2 * np.cos(phi) * nnn_cos
    dz = m - 2.0 * t2 * np.sin(phi) * nnn_sin
    return eps, dx, dy, dz


def haldane_hamiltonian(kx, ky, m=0.0, t1=1.0, t2=0.15, phi=np.pi / 2):
    """The 2x2 Haldane Bloch Hamiltonian at a single k-point."""
    eps, dx, dy, dz = haldane_components(kx, ky, m, t1, t2, phi)
    return eps * IDENTITY + dx * SIGMA_X + dy * SIGMA_Y + dz * SIGMA_Z


def haldane_hamiltonian_grid(KX, KY, m=0.0, t1=1.0, t2=0.15, phi=np.pi / 2):
    """Haldane Hamiltonian on a whole grid: shape (N, N, 2, 2)."""
    eps, dx, dy, dz = haldane_components(KX, KY, m, t1, t2, phi)
    return _assemble(eps, dx, dy, dz)


def haldane_energies(kx, ky, m=0.0, t1=1.0, t2=0.15, phi=np.pi / 2):
    """Bands for the Haldane model: eps +/- |d| (the eps shift is included)."""
    eps, dx, dy, dz = haldane_components(kx, ky, m, t1, t2, phi)
    mag = np.sqrt(dx**2 + dy**2 + dz**2)
    return eps - mag, eps + mag


def haldane_dirac_masses(m=0.0, t2=0.15, phi=np.pi / 2):
    """
    Return (m_K, m_Kp), the Dirac masses at K and K'.

    Their signs decide the phase: opposite signs -> topological (C = +/-1),
    same sign -> trivial (C = 0).
    """
    shift = 3.0 * _SQRT3 * t2 * np.sin(phi)
    return m - shift, m + shift


# ---------------------------------------------------------------------------
# Band-structure paths through high-symmetry points
# ---------------------------------------------------------------------------
def _path_through(corners, labels, n_per_segment):
    kxs, kys, ticks = [], [], [0]
    for a, b in zip(corners[:-1], corners[1:]):
        seg_x = np.linspace(a[0], b[0], n_per_segment, endpoint=False)
        seg_y = np.linspace(a[1], b[1], n_per_segment, endpoint=False)
        kxs.extend(seg_x)
        kys.extend(seg_y)
        ticks.append(len(kxs))
    kxs.append(corners[-1][0])
    kys.append(corners[-1][1])
    return np.array(kxs), np.array(kys), ticks, labels


def qwz_path(n_per_segment=160):
    """Square-lattice path Gamma -> X -> M -> Gamma."""
    corners = [(0, 0), (np.pi, 0), (np.pi, np.pi), (0, 0)]
    labels = [r"$\Gamma$", "X", "M", r"$\Gamma$"]
    return _path_through(corners, labels, n_per_segment)


def haldane_path(n_per_segment=160):
    """
    Honeycomb path Gamma -> K -> M -> K' -> Gamma.

    Passing through both K and K' is deliberate: at a topological transition the
    gap closes at only ONE of them, so seeing both makes the transition obvious.
    """
    corners = [(0, 0), tuple(K_POINT), tuple(M_POINT), tuple(KP_POINT), (0, 0)]
    labels = [r"$\Gamma$", "K", "M", "K'", r"$\Gamma$"]
    return _path_through(corners, labels, n_per_segment)

"""
Berry curvature and Chern number, by the Fukui-Hatsugai-Suzuki (FHS) method.

This engine is model-agnostic: it takes a stack of Bloch Hamiltonians already
evaluated on a grid (shape (N, N, 2, 2)) and returns the Chern number of the
lower band. It never needs to know whether those matrices came from the QWZ
model or the Haldane model. That is the whole idea: the topology machinery is
universal.

The gauge problem, and the fix
------------------------------
The honest definitions are a Berry connection A(k) = i <u|grad_k|u> (a gradient),
its curl Omega = curl A (the Berry curvature), and the Chern number
C = (1/2pi) * integral of Omega over the Brillouin zone (a surface integral),
which is guaranteed to be an integer. But a numerical eigensolver returns each
eigenvector |u(k)> with an arbitrary phase, so the gradient is meaningless
point to point. FHS rewrites the integral using overlaps between neighboring
k-points ("link variables"), and those arbitrary phases cancel exactly. The
result is an integer automatically.

Everything below is vectorized with NumPy: one batched diagonalization for the
whole grid, then array shifts (np.roll) to form the plaquette loops. That makes
sweeping a 2D phase diagram fast.
"""

import numpy as np


def lower_band_states(hamiltonian_grid):
    """
    Diagonalize a whole grid of 2x2 Hamiltonians at once and return the
    normalized lower-band eigenvector at every k-point.

    Parameters
    ----------
    hamiltonian_grid : ndarray, shape (N, N, 2, 2), complex

    Returns
    -------
    ndarray, shape (N, N, 2), complex
        states[i, j] is the lower-band eigenvector at grid point (i, j).
    """
    # np.linalg.eigh is batched: it diagonalizes every 2x2 block along the last
    # two axes in one call, returning eigenvalues ascending. Column 0 is the
    # lower band.
    _, evecs = np.linalg.eigh(hamiltonian_grid)
    return evecs[..., :, 0]


def _links(states, axis):
    """
    Link variables along one grid direction: U(k) = <u(k)|u(k+step)> normalized
    to unit modulus, for every k at once. np.roll implements the periodic
    "k+step" shift.
    """
    shifted = np.roll(states, -1, axis=axis)
    overlap = np.sum(np.conj(states) * shifted, axis=-1)  # <u(k)|u(k+step)>
    return overlap / np.abs(overlap)


def berry_curvature_field(hamiltonian_grid):
    """
    Per-plaquette Berry curvature ("field strength") on the grid.

    Returns
    -------
    ndarray, shape (N, N), real
        Summing this and dividing by 2pi gives the Chern number.
    """
    states = lower_band_states(hamiltonian_grid)

    U_x = _links(states, axis=0)   # links along kx
    U_y = _links(states, axis=1)   # links along ky

    # Wilson loop around each plaquette, built entirely with shifts:
    #   U_x(k) * U_y(k + x) * conj(U_x(k + y)) * conj(U_y(k))
    U_y_shift_x = np.roll(U_y, -1, axis=0)
    U_x_shift_y = np.roll(U_x, -1, axis=1)
    wilson = U_x * U_y_shift_x * np.conj(U_x_shift_y) * np.conj(U_y)

    # Field strength = argument of the loop, in (-pi, pi].
    return np.angle(wilson)


def chern_number(hamiltonian_grid):
    """
    The Chern number of the lower band, summed from the curvature field.

    Returns
    -------
    float
        Very close to an integer; round it for the exact value.
    """
    field = berry_curvature_field(hamiltonian_grid)
    return field.sum() / (2 * np.pi)

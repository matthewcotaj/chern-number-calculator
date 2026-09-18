"""
topo: a small toolkit for two-band topological lattice models.

models.py has the QWZ and Haldane Bloch Hamiltonians (energies, geometry, BZ
paths). berry.py has the Fukui-Hatsugai-Suzuki Chern number calculation, which
is model-agnostic -- it only ever sees a Hamiltonian evaluated on a grid, so
both models plug into the same engine.
"""

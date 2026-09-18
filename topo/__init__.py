"""
topo: a small toolkit for two-band topological lattice models.

It contains two models and one shared engine for the topology:

models
------
- Qi-Wu-Zhang (QWZ): the simplest Chern insulator, on a square lattice.
- Haldane: the honeycomb-lattice model that started the whole field.

engine
------
- berry: the Fukui-Hatsugai-Suzuki Chern-number calculation. It works on *any*
  two-band model, because it only ever sees the Hamiltonian evaluated on a grid.
  Building both models on top of the same engine is the point: the topology
  machinery is universal, only the Hamiltonian changes.

modules
-------
models : both Bloch Hamiltonians, their energies, geometry, and BZ paths
berry  : Berry curvature and Chern number (model-agnostic)
"""

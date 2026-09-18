"""
Tests: proof that both models are correct.

    python test_models.py

Checks:
  1. QWZ Chern number is the known integer in each phase.
  2. QWZ gap closes at u = -2, 0, 2.
  3. Haldane Chern number is +/-1 in the topological lobes and 0 outside,
     with the phase boundary at |m| = 3*sqrt(3)*t2*sin(phi).
  4. Haldane gap closes at exactly the predicted mass at one Dirac point.

If every line says PASS, the calculators work.
"""

import numpy as np

from topo.models import (brillouin_grid, honeycomb_grid, qwz_hamiltonian_grid,
                         qwz_energies, haldane_hamiltonian_grid,
                         haldane_energies, haldane_dirac_masses,
                         K_POINT, KP_POINT)
from topo.berry import chern_number
from topo.ribbon import edge_mode_crossings


def check(name, condition):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}")
    return condition


def test_qwz_chern():
    print("QWZ Chern number in each phase:")
    ok = True
    KX, KY = brillouin_grid(60)
    for u, expected in {-3.0: 0, -1.0: -1, 1.0: 1, 3.0: 0}.items():
        C = chern_number(qwz_hamiltonian_grid(KX, KY, u))
        ok &= check(f"u = {u:+.0f} gives C = {expected:+d} (got {C:+.3f})",
                    round(C) == expected)
    return ok


def test_qwz_gap():
    print("\nQWZ gap:")
    ok = True
    kv = np.linspace(0, 2 * np.pi, 120, endpoint=False)
    for u_c in (-2.0, 0.0, 2.0):
        gap = min(qwz_energies(kx, ky, u_c)[1] - qwz_energies(kx, ky, u_c)[0]
                  for kx in kv for ky in kv)
        ok &= check(f"gap ~ 0 at u = {u_c:+.0f} (min = {gap:.3f})", gap < 0.15)
    return ok


def test_haldane_chern():
    print("\nHaldane Chern number (t2 = 0.15, phi = pi/2):")
    ok = True
    KX, KY = honeycomb_grid(72)
    # boundary at |m| = 3 sqrt(3) t2 ~ 0.779
    cases = {-1.2: 0, -0.4: -1, 0.0: -1, 0.4: -1, 1.2: 0}
    for m, expected in cases.items():
        C = chern_number(haldane_hamiltonian_grid(KX, KY, m=m, t2=0.15, phi=np.pi / 2))
        ok &= check(f"m = {m:+.1f} gives |C| = {abs(expected)} (got {C:+.3f})",
                    round(C) == expected)
    # flipping phi flips the sign of C
    Cp = chern_number(haldane_hamiltonian_grid(KX, KY, m=0.0, t2=0.15, phi=np.pi / 2))
    Cm = chern_number(haldane_hamiltonian_grid(KX, KY, m=0.0, t2=0.15, phi=-np.pi / 2))
    ok &= check(f"phi -> -phi flips C ({Cp:+.2f} vs {Cm:+.2f})",
                round(Cp) == -round(Cm) and round(Cp) != 0)
    return ok


def test_haldane_gap_closing():
    print("\nHaldane gap closes at the predicted mass:")
    ok = True
    t2, phi = 0.15, np.pi / 2
    m_c = 3 * np.sqrt(3) * t2 * np.sin(phi)   # ~0.779
    # At m = +m_c the K point should be gapless; K' should stay open.
    gap_K = (haldane_energies(*K_POINT, m=m_c, t2=t2, phi=phi)[1]
             - haldane_energies(*K_POINT, m=m_c, t2=t2, phi=phi)[0])
    gap_Kp = (haldane_energies(*KP_POINT, m=m_c, t2=t2, phi=phi)[1]
              - haldane_energies(*KP_POINT, m=m_c, t2=t2, phi=phi)[0])
    ok &= check(f"gap at K ~ 0 when m = m_c = {m_c:.3f} (gap = {gap_K:.3f})", gap_K < 1e-6)
    ok &= check(f"gap at K' stays open there (gap = {gap_Kp:.3f})", gap_Kp > 0.5)
    # Dirac-mass bookkeeping matches the closed form.
    mK, mKp = haldane_dirac_masses(m=m_c, t2=t2, phi=phi)
    ok &= check(f"m_K ~ 0 at the boundary (m_K = {mK:.3e})", abs(mK) < 1e-9)
    return ok


def test_bulk_boundary():
    print("\nBulk-boundary correspondence (ribbon edge modes = 2|C|):")
    ok = True
    KX, KY = honeycomb_grid(60)
    for m in (0.0, 0.4, 1.4):
        C = round(chern_number(haldane_hamiltonian_grid(KX, KY, m=m, t2=0.15, phi=np.pi / 2)))
        crossings = edge_mode_crossings(width=20, m=m, t2=0.15, phi=np.pi / 2)
        ok &= check(f"m = {m:+.1f}: C = {C:+d}, edge crossings = {crossings} (want {2*abs(C)})",
                    crossings == 2 * abs(C))
    return ok


if __name__ == "__main__":
    all_ok = True
    all_ok &= test_qwz_chern()
    all_ok &= test_qwz_gap()
    all_ok &= test_haldane_chern()
    all_ok &= test_haldane_gap_closing()
    all_ok &= test_bulk_boundary()
    print("\n" + ("ALL TESTS PASSED" if all_ok else "SOME TESTS FAILED"))

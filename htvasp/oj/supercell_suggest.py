# -*- coding: utf-8 -*-
"""Suggest supercell size for OstravaJ based on structure + dist_cutoff.

Two competing constraints:
  1)  Uniqueness of J-pair vectors — need supercell to distinguish
      atom pairs at distances up to dist_cutoff.
  2)  Tractable search space — OJ's config generator searches 2^{N_atoms}
      spin configurations; N_atoms > 128 becomes expensive on a single node.

Empirical calibration (from dry-runs on BCC 2-atom, Fe at 2.866 A):

    supercell   atoms   result
    (4,4,2)      64    ✅  7x7 full-rank, solved
    (3,3,3)      54    ❌  only 4 configs for 7 unknowns

The heuristic:
  1.  n_i >= ceil(2 * dist_cutoff / a_min) — geometric bound
  2.  Reduce if search space > 2^{80} — OJ tractability
  3.  sum(n_i) >= 10 — enough independent configs
  4.  Never shrink below a validated minimum (4,4,2)
"""

import math


def suggest_supercell(struct, dist_cutoff=6.0, atoms_per_cell=None):
    """Suggest supercell size for OstravaJ.

    Args:
        struct: pymatgen Structure.
        dist_cutoff: Distance cutoff for J-pairs (Angstrom).
        atoms_per_cell: Number of atoms in the primitive cell.

    Returns:
        (nx, ny, nz), reason_string
    """
    if atoms_per_cell is None:
        atoms_per_cell = len(struct)

    a, b, c = struct.lattice.abc

    # Step 1: geometric bound
    n_req = int(math.ceil(2.0 * dist_cutoff / max(a, b, c)))
    nx = max(n_req, 4)
    ny = max(n_req, 4)
    nz = max(n_req, 2)

    def _natoms():
        return atoms_per_cell * nx * ny * nz

    # Step 2: cap search space to 2^80
    while 2 ** _natoms() > 2 ** 80 and nx > 2 and ny > 2 and nz > 1:
        if nz >= nx and nz >= ny and nz > max(2, min(nx, ny)):
            nz = max(nz - 1, 2)
        elif ny >= nx and ny > 2:
            ny -= 1
        elif nx > 2:
            nx -= 1
        else:
            break

    # Step 3: ensure sum(n_i) >= 10
    while nx + ny + nz < 10 and _natoms() <= 256:
        if nz <= nx and nz <= ny and nz < 6:
            nz += 1
        elif ny <= nx and ny < 6:
            ny += 1
        elif nx < 6:
            nx += 1
        else:
            break

    # Step 4: never below the smallest validated supercell
    # (4,4,2) confirmed via dry-run on BCC 2-atom Fe (64 atoms, 7x7 full rank)
    if nx < 4: nx = 4
    if ny < 4: ny = 4
    if nz < 2: nz = 2

    n_total = _natoms()
    reasons = [
        f"({nx},{ny},{nz}) = {n_total} atoms; "
        f"2^{n_total} ~ 10^{n_total * math.log10(2):.0f} search space",
    ]

    # J-pair coverage
    # OJ needs the supercell to distinguish atom pairs up to dist_cutoff.
    # The strict condition (unique images) is n_i*a_i >= 2*d_max.
    # The minimum condition (vector fits in box) is n_i*a_i > d_max.
    labels = [("a", a), ("b", b), ("c", c)]
    for (label, a_i), n_i in zip(labels, [nx, ny, nz]):
        extent = n_i * a_i
        if extent >= 2 * dist_cutoff:
            reasons.append(f"{label}: {extent:.1f}A >= 2*d_max=12.0A (optimal)")
        elif extent >= dist_cutoff:
            reasons.append(f"{label}: {extent:.1f}A in [{dist_cutoff:.0f}, 12.0)A (adequate)")
        else:
            reasons.append(f"{label}: {extent:.1f}A < d_max (may need larger)")

    if nx + ny + nz >= 10:
        reasons.append(f"sum(n_i)={nx}+{ny}+{nz}={nx+ny+nz} >= 10 OK")

    return (nx, ny, nz), "; ".join(reasons)

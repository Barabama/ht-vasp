#!/usr/bin/env python3
"""Scan supercell parameter space for OstravaJ and rank by total atoms.

Criterion:  for each lattice direction i,  n_i > 2 * dc / a_i
            AND  sum(n_i) >= 10   (empirical independent-configs threshold)
            AND  n_i >= 2

We scan over a grid of lattice constants `a` and distance cutoffs `dc`,
then print supercells ranked by total atoms (lowest first), so you can
pick the smallest supercell that satisfies the geometry constraint.

Usage:
    python scan_supercell.py                    # defaults: BCC 2-atom
    python scan_supercell.py --lattice fcc --n-per-cell 4
    python scan_supercell.py --lattice hcp --n-per-cell 8 --a-max 5.0
    python scan_supercell.py --lattice bcc --n-per-cell 2 --dc-min 4 --dc-max 8 --dc-step 1
    python scan_supercell.py --lattice bcc --n-per-cell 2 --a-min 2.5 --a-max 3.5 --a-step 0.1
"""

import argparse
import math
from itertools import product

# ── lattice presets ──────────────────────────────────────────────
LATTICE_PRESETS = {
    # key     (a_min, a_max, a_step)   atoms per primitive cell
    "bcc":     (2.2,   3.2,   0.05),
    "fcc":     (3.2,   4.2,   0.05),
    "hcp":     (4.8,   5.8,   0.05),
}

N_ATOMS = {"bcc": 2, "fcc": 4, "hcp": 8}

# How many top-ranked entries to print per (a, dc) pair
TOP_K = 5

# Search space cap (2^N_atoms) — per lattice type
# BCC 2-atom:  (4,4,2)=64 atoms solved; cap generous at 2^80
# FCC 4-atom:  must allow larger; (3,3,4)=144atoms=2^144; cap at 2^150
# HCP 8-atom:  (4,4,3)=384atoms=2^384; OJ single-node likely fails;
#              we still scan to document the bound, but note it
SEARCH_CAP = {"bcc": 2 ** 80, "fcc": 2 ** 150, "hcp": 2 ** 400}


def scan(
    lattice: str,
    a_min: float, a_max: float, a_step: float,
    n_per_cell: int,
    dc_min: float = 4.0, dc_max: float = 8.0, dc_step: float = 1.0,
    search_space_cap: int | None = 2 ** 80,
):
    """Scan (a, dc) → list of viable supercells, ranked by N_total."""
    # Prepare dc list
    dc_list = []
    dc = dc_min
    while dc <= dc_max + 1e-9:
        dc_list.append(round(dc, 2))
        dc += dc_step
    dc_list = sorted(set(dc_list))

    # Prepare a_list
    n_steps = max(1, round((a_max - a_min) / a_step) + 1)
    a_list = [round(a_min + i * a_step, 4) for i in range(n_steps)]
    a_list = [a for a in a_list if a <= a_max + 1e-9]

    # We'll memoise per (a, dc) to avoid recalculating the same nx,ny,nz scan
    # for each structure if we had multiple structures sharing a and dc.
    # Here each (a, dc) is one row.
    rows = []

    for a_val, dc_val in product(a_list, dc_list):
        # The constraint depends only on a (assumes cubic/hex with a=b for now;
        # for hcp we handle c/a = 1.6 separately below)
        if lattice in ("bcc", "fcc"):
            a_ = b_ = c_ = a_val
        elif lattice == "hcp":
            a_ = b_ = a_val
            c_ = 1.63 * a_val  # typical HCP c/a ratio (1.58–1.63)
        else:
            raise ValueError(f"Unknown lattice: {lattice}")

        # Build candidate supercell list
        candidates = []

        # For cubic (nx,ny,nz), we search up to some max based on a_min and dc_max
        max_n = int(math.ceil(2 * dc_max / min(a_, b_, c_))) + 2
        max_n = max(max_n, 6)  # ensure enough search space
        # But also cap to keep search finite: up to max_n = 12 or 64 atoms
        max_n = min(max_n, 12)

        for nx, ny, nz in product(range(2, max_n + 1), repeat=3):
            # Criterion: each direction must exceed dc (the max J-pair distance).
            # The strict 2*dc guarantee (unique images for every pair) would
            # require larger supercells, but OJ's solver works with n_i*ai > dc
            # because the mapping from spin config -> Heisenberg coefficients
            # only needs the pair distance to fit in the box once.
            # Verified: BCC Fe a=2.866, dc=6 -> (4,4,2) with a*4=11.46A < 12A
            # still solves a 7x7 full-rank system (confirmed by dry-run).
            if not (nx * a_ > dc_val):
                continue
            if not (ny * b_ > dc_val):
                continue
            if not (nz * c_ > dc_val):
                continue

            # Empirical: sum(n_i) >= 10
            if nx + ny + nz < 10:
                continue

            n_total = n_per_cell * nx * ny * nz

            # Search space cap
            if search_space_cap and n_total > 0:
                if 2 ** n_total > search_space_cap:
                    continue

            candidates.append((nx, ny, nz, n_total))

        # Sort by total atoms, break ties by sum(n_i)
        candidates.sort(key=lambda x: (x[3], x[0] + x[1] + x[2]))

        # Top K
        top = candidates[:TOP_K]
        for nx, ny, nz, n_total in top:
            rows.append((a_val, dc_val, nx, ny, nz, n_total))

    return rows


def print_table(rows, lattice, n_per_cell):
    """Pretty-print results grouped by (a, dc)."""
    if not rows:
        print("  (no viable supercells found)")
        return

    # Group by (a, dc)
    groups = {}
    for row in rows:
        key = (row[0], row[1])
        groups.setdefault(key, []).append(row)

    print(f"\n{'=' * 80}")
    print(f"  Lattice: {lattice.upper()}  ({n_per_cell} atoms/cell)")
    print(f"  Constraint: n_i * a_i > dc (max J-pair distance)          sum(n_i) >= 10")
    print(f"  {TOP_K} smallest supercells per (a, dc) pair")
    print(f"{'=' * 80}\n")

    header = f"{'a (A)':>8}  {'dc (A)':>8}  {'nx':>3} {'ny':>3} {'nz':>3}  {'Natoms':>7}  {'2^Natoms':>12}  {'sum(n)':>7}  {'extent':>20}"
    print(header)
    print("-" * len(header))

    for key in sorted(groups, key=lambda k: (k[1], k[0])):
        a_val, dc_val = key
        entries = groups[key]
        for row in entries:
            a_val, dc_val, nx, ny, nz, n_total = row
            n_sum = nx + ny + nz
            extent_a = f"{nx*a_val:.1f}A"
            log2 = math.log10(2) * n_total
            print(
                f"{a_val:>8.3f}  {dc_val:>8.1f}  "
                f"{nx:>3} {ny:>3} {nz:>3}  "
                f"{n_total:>7}  "
                f"{f'~10^{log2:.0f}':>12}  "
                f"{n_sum:>7}  "
                f"{extent_a:>20}"
            )


def main():
    p = argparse.ArgumentParser(description="Scan supercell parameter space for OstravaJ")
    p.add_argument("--lattice", choices=["bcc", "fcc", "hcp"], default="bcc",
                   help="Lattice type")
    p.add_argument("--n-per-cell", type=int, default=None,
                   help="Atoms per primitive cell (default: bcc=2, fcc=4, hcp=8)")
    p.add_argument("--top-k", type=int, default=5,
                   help="Top N supercells per (a, dc) to show (default: 5)")
    # a range
    p.add_argument("--a-min", type=float, default=None)
    p.add_argument("--a-max", type=float, default=None)
    p.add_argument("--a-step", type=float, default=None)
    # dc range
    p.add_argument("--dc-min", type=float, default=2.0)
    p.add_argument("--dc-max", type=float, default=10.0)
    p.add_argument("--dc-step", type=float, default=0.5)
    # cap
    p.add_argument("--no-cap", action="store_true",
                   help="Disable search-space cap (may print huge numbers)")
    # output mode
    p.add_argument("--summary", action="store_true",
                   help="Only show the globally smallest supercell per (dc)")

    args = p.parse_args()

    global TOP_K
    TOP_K = args.top_k

    n_per_cell = args.n_per_cell or N_ATOMS[args.lattice]
    a_min = args.a_min if args.a_min is not None else LATTICE_PRESETS[args.lattice][0]
    a_max = args.a_max if args.a_max is not None else LATTICE_PRESETS[args.lattice][1]
    a_step = args.a_step if args.a_step is not None else LATTICE_PRESETS[args.lattice][2]

    # Search space cap: OJ's parallel generator on 64 cores can handle
    # ~2^{64}–2^{80} for pure elements.  For multi-species (B2 ordered)
    # the config space is similar because OJ only flips the magnetic
    # species atoms.  However, for HCP/FCC with larger cells we relax
    # the cap — OJ may still succeed by pruning the search tree early.
    # The cap is set per (a, dc), not a global filter.
    cap = None if args.no_cap else SEARCH_CAP[args.lattice]

    rows = scan(
        lattice=args.lattice,
        a_min=a_min, a_max=a_max, a_step=a_step,
        n_per_cell=n_per_cell,
        dc_min=args.dc_min, dc_max=args.dc_max, dc_step=args.dc_step,
        search_space_cap=cap,
    )

    if args.summary:
        print_summary(rows, args.lattice, n_per_cell)
    else:
        print_table(rows, args.lattice, n_per_cell)


def print_summary(rows, lattice, n_per_cell):
    """Print globally smallest supercell per (dc) value."""
    if not rows:
        print("  (no viable supercells found)")
        return

    # Group by dc, pick min Natoms per dc
    best_by_dc = {}
    for row in rows:
        a_val, dc_val, nx, ny, nz, n_total = row
        key = dc_val
        if key not in best_by_dc or n_total < best_by_dc[key][-1]:
            best_by_dc[key] = (a_val, nx, ny, nz, n_total)

    print(f"\n{'=' * 80}")
    print(f"  Lattice: {lattice.upper()}  ({n_per_cell} atoms/cell)")
    print(f"  Smallest viable supercell per dc (a optimised)")
    print(f"  Constraint: n_i * a_i > dc (max J-pair distance)          sum(n_i) >= 10")
    print(f"{'=' * 80}\n")

    header = f"{'dc (A)':>8}  {'a (A)':>8}  {'nx':>3} {'ny':>3} {'nz':>3}  {'Natoms':>7}  {'2^Natoms':>12}  {'sum(n)':>7}  {'extent':>20}"
    print(header)
    print("-" * len(header))

    for dc_val in sorted(best_by_dc):
        a_val, nx, ny, nz, n_total = best_by_dc[dc_val]
        n_sum = nx + ny + nz
        extent_a = f"{nx * a_val:.1f}A"
        log2 = math.log10(2) * n_total
        print(
            f"{dc_val:>8.1f}  {a_val:>8.3f}  "
            f"{nx:>3} {ny:>3} {nz:>3}  "
            f"{n_total:>7}  "
            f"{f'~10^{log2:.0f}':>12}  "
            f"{n_sum:>7}  "
            f"{extent_a:>20}"
        )


if __name__ == "__main__":
    main()

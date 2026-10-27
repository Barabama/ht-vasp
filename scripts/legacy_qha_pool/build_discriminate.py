#!/usr/bin/env python
"""Discrimination experiment: why does the 2-atom production cell show a much
larger k-mesh error than the 16-atom phonon supercell at comparable reciprocal
sampling density?

The review rewrote the question so that no cause is assumed up front. The
candidate mechanisms it lists map onto concrete knobs here:

  (1) is the sampling density really equivalent?      -> matched-volume pair, both cells
  (2) BZ folding / irreducible-count differences      -> record NKPTS + IBZKPT
  (3) is ISMEAR=1's sigma->0 extrapolation fragile?   -> same grid, s1 vs sm5
  (4) does symmetry reduction give different sampling? -> even-mesh shift 0 vs 0.5
  (5) different Fermi-surface sampling                 -> record E-fermi, magnetization

Design (all at V/V0 = 1.0 unless stated, so the two cells are an EXACT
2x2x2 repeat pair -- same shape, same space group Pm-3m, same volume per
2-atom formula unit):

  A. a2_f4   2-atom, a = A0*(1+eps4) = frame 4 of the production QHA grid
             (V/V0 = 0.97872). Completes the frame-4 k-series that produced
             the 12.7 meV/atom anomaly.
  B. a2_v0   2-atom, a = A0 = 3.2246015525754168 A  (V/V0 = 1.0)
  C. a16_v0  16-atom, a = 2*A0  (V/V0 = 1.0, exactly 2x2x2 of B)

For every case: ISMEAR=1/SIGMA=0.2 (production) and ISMEAR=-5 (tetrahedron).
Even meshes are run at BOTH shift 0 0 0 and shift 0.5 0.5 0.5, because for an
even n the MP grid with shift 0 does not contain Gamma while shift 0.5 does --
so a "density-equivalent" comparison of even meshes is only meaningful once
that is controlled for.

Recipe: identical to the alignment test (build_align_inputs.incar_text):
ENCUT=520, ALGO=Normal, EDIFF=1e-6, ISYM=2, NCORE=8, NELM=120, NELMIN=6,
ISPIN=2, MAGMOM Al 1.0 / Ta 3.0, POTCAR PBE(52) standard variants.
The 16-atom cases use the SAME INCAR text, so EDIFF/ISYM/NCORE do not differ
between the cells -- only the cell and the mesh do.

NOT done deliberately: 2x1x1. It changes shape, symmetry, BZ and band
degeneracies, so it is not a clean "electron count" control; 2x2x2 (Pm-3m,
same shape) is. The review asked for exactly this substitution.

Writes only under --out/disc. Refuses to overwrite any existing INCAR.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from build_align_inputs import A0, FRAMES, STRAIN, b2_cell, write_case
from build_d0_inputs import b2 as b2_super          # 2x2x2 (or n x n x n) supercell

SMEARINGS = {"s1": (1, 0.2), "sm5": (-5, None)}

# frame 4 of the production QHA volume grid
EPS4 = float(np.linspace(STRAIN[0], STRAIN[1], FRAMES)[3])
A_F4 = A0 * (1.0 + EPS4)

# (kgrid, shift) node lists
A2_F4 = [(6, 0.0), (6, 0.5), (7, 0.0), (8, 0.0), (8, 0.5), (9, 0.0)]
A2_V0 = [(6, 0.0), (6, 0.5), (7, 0.0), (8, 0.0), (8, 0.5), (9, 0.0)]
A16_V0 = [(3, 0.0), (4, 0.0), (4, 0.5), (5, 0.0), (6, 0.0), (7, 0.0), (9, 0.0)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--blocks", default="a2_f4,a2_v0,a16_v0",
                    help="comma-separated subset to (re)build")
    args = ap.parse_args()
    root = args.out / "disc"
    root.mkdir(parents=True, exist_ok=True)
    want = set(args.blocks.split(","))

    cases: list[dict] = []

    def add(block: str, struct, nodes, cell_label: str, vol_scale: float):
        for k, sh in nodes:
            for tag, (ism, sig) in SMEARINGS.items():
                shs = f"s{int(sh * 2)}"          # 0 or 5, i.e. shift 0 or 0.5
                name = f"k{k:02d}_{shs}_{tag}"
                cases.append(write_case(
                    root / block / name, struct, k,
                    f"disc-{block}-k{k}-{shs}-{tag}", ism, sig,
                    {"block": block, "cell": cell_label, "frame": None,
                     "volume_angstrom3": float(struct.volume),
                     "volume_scale_vs_v0": vol_scale, "shift_kind": shs},
                    shift=(sh, sh, sh)))

    # A. 2-atom at frame 4 (the anomalous volume)
    if "a2_f4" in want:
        add("a2_f4", b2_cell(A_F4), A2_F4, "a2", (A_F4 / A0) ** 3)
    # B. 2-atom at V/V0 = 1
    if "a2_v0" in want:
        add("a2_v0", b2_cell(A0), A2_V0, "a2", 1.0)
    # C. 16-atom (2x2x2 of B) at V/V0 = 1
    if "a16_v0" in want:
        add("a16_v0", b2_super(2), A16_V0, "a16", 1.0)

    bad = [c for c in cases if not c["ok"]]
    for c in cases:
        mark = "OK " if c["ok"] else "BAD"
        print(f"  {mark} {c['dir'].split('disc/')[-1]:<22} natoms={c['natoms']:>3} "
              f"V={c['volume_angstrom3']:8.4f} k={c['kgrid']:>2} shift={c['shift']} "
              f"ISMEAR={c['ismear']:>2}")
        for p in c["problems"]:
            print(f"        !! {p}")
    (root / "validation.json").write_text(json.dumps(cases, indent=2) + "\n")
    print(f"\n{len(cases) - len(bad)}/{len(cases)} validated")
    if bad:
        raise SystemExit("refusing to submit: some inputs failed validation")


if __name__ == "__main__":
    main()

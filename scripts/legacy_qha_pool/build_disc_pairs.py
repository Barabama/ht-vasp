#!/usr/bin/env python
"""Add the parity-matched density pairs that decide the cell-size question.

The discrimination run left one clean parity-matched pair (2-atom k=8 <->
16-atom k=4, both even so both exclude Gamma) and one parity-MISMATCHED pair
(2-atom k=6 even/excludes Gamma <-> 16-atom k=3 odd/includes Gamma). The
matched pair agreed to 0.011 meV/atom; the mismatched pair differed by 4.4.
That is suggestive but rests on a single pair, so this adds two more matched
pairs:

    density 12 : 2-atom k=12 (even) <-> 16-atom k=6 (even)   [16@6 exists]
    density 16 : 2-atom k=16 (even) <-> 16-atom k=8 (even)

and 2-atom k=10 to fill the series between them. All at V/V0 = 1.0, both
smearings, shift 0.

Writes only under --out/disc, refusing to overwrite.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from build_align_inputs import A0, b2_cell, write_case
from build_d0_inputs import b2 as b2_super

SMEARINGS = {"s1": (1, 0.2), "sm5": (-5, None)}
A2_EXTRA = (10, 12, 16)          # even, shift 0 -> no Gamma
A16_EXTRA = (8,)                 # even, shift 0 -> no Gamma, pairs with 2-atom k=16


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    root = args.out / "disc"

    cases = []
    for k in A2_EXTRA:
        for tag, (ism, sig) in SMEARINGS.items():
            cases.append(write_case(
                root / "a2_v0" / f"k{k:02d}_s0_{tag}", b2_cell(A0), k,
                f"disc-a2_v0-k{k}-s0-{tag}", ism, sig,
                {"block": "a2_v0", "cell": "a2", "frame": None,
                 "volume_angstrom3": float(b2_cell(A0).volume),
                 "volume_scale_vs_v0": 1.0, "shift_kind": "s0"},
                shift=(0.0, 0.0, 0.0)))
    for k in A16_EXTRA:
        for tag, (ism, sig) in SMEARINGS.items():
            cases.append(write_case(
                root / "a16_v0" / f"k{k:02d}_s0_{tag}", b2_super(2), k,
                f"disc-a16_v0-k{k}-s0-{tag}", ism, sig,
                {"block": "a16_v0", "cell": "a16", "frame": None,
                 "volume_angstrom3": float(b2_super(2).volume),
                 "volume_scale_vs_v0": 1.0, "shift_kind": "s0"},
                shift=(0.0, 0.0, 0.0)))

    for c in cases:
        mark = "OK " if c["ok"] else "BAD"
        print(f"  {mark} {c['dir'].split('disc/')[-1]:<22} natoms={c['natoms']:>3} "
              f"k={c['kgrid']:>2} ISMEAR={c['ismear']:>2}")
        for p in c["problems"]:
            print(f"        !! {p}")
    bad = [c for c in cases if not c["ok"]]
    print(f"\n{len(cases) - len(bad)}/{len(cases)} validated")
    if bad:
        raise SystemExit("refusing to submit")


if __name__ == "__main__":
    main()

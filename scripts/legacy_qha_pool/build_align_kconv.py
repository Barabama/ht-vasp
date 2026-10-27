#!/usr/bin/env python
"""k-mesh convergence on the PRODUCTION cell (2-atom B2 AlTa conventional cell).

The alignment test found that E(k=11) - E(k=7) on the 2-atom cell is
-13.8 meV/atom (ISMEAR=1), i.e. an order of magnitude LARGER than the whole
smearing effect (<=1 meV/atom). That is not surprising once the cell is
accounted for: the earlier KGRID scan ran on the 2x2x2 (16-atom) supercell,
whose BZ is 8x smaller, so a mesh of n there samples 2n on the conventional
cell. In other words the KGRID scan's "k=11" was effectively 22x22x22 on the
cell the EOS is actually fitted on.

This block pins down the mesh the production E(V) statics need, at the middle
volume point, ISMEAR=1 (the intended production smearing).

Writes only under --out/kconv.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from build_align_inputs import (A0, FRAMES, KMAIN, STRAIN, b2_cell, write_case)

MESHES = (5, 9, 13, 15, 17)          # 7 and 11 already exist (kmesh/ and primary/)
MESHES_SM5 = (9, 13)                 # check whether dE(-5 -> s1) shrinks with mesh


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    out = args.out / "kconv"
    out.mkdir(parents=True, exist_ok=True)

    eps_l = np.linspace(STRAIN[0], STRAIN[1], FRAMES)
    idx = 4
    eps = float(eps_l[idx - 1])
    a = A0 * (1.0 + eps)
    st = b2_cell(a)
    v0 = A0 ** 3
    meta = {"block": "kconv", "frame": idx, "eps": eps, "a_angstrom": a,
            "volume_angstrom3": float(st.volume),
            "volume_scale_vs_v0": float(st.volume) / v0}

    cases = []
    for k in MESHES:
        cases.append(write_case(out / f"f{idx:02d}_k{k}_s1", st, k,
                                f"align-kconv-k{k}-s1", 1, 0.2, meta))
    for k in MESHES_SM5:
        cases.append(write_case(out / f"f{idx:02d}_k{k}_sm5", st, k,
                                f"align-kconv-k{k}-sm5", -5, None, meta))

    for c in cases:
        mark = "OK " if c["ok"] else "BAD"
        print(f"  {mark} {c['dir']:<44} V={c['volume_angstrom3']:.4f} "
              f"k={c['kgrid']:>2} ISMEAR={c['ismear']:>2}")
        for p in c["problems"]:
            print(f"        !! {p}")
    (out / "validation.json").write_text(json.dumps(cases, indent=2) + "\n")
    bad = [c for c in cases if not c["ok"]]
    print(f"\n{len(cases) - len(bad)}/{len(cases)} validated")
    if bad:
        raise SystemExit("refusing to submit")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""k=15 repeat of the 5-volume scan, for ISMEAR=1 and ISMEAR=-5 only.

Why: the alignment test's k-series on the 2-atom production cell showed the
residual k-mesh error at k=11 is ~1 meV/atom (k=11 vs k=15 = -1.11 meV/atom),
which is the same size as the whole dE(-5 -> s1) range it was measuring. So the
primary 5-volume verdict has to be re-taken at a denser mesh.

This also yields a second E(V) curve (k=15, ISMEAR=1) that can be fitted and
compared against the k=11 curve, which is the direct answer to "does the k-mesh
error distort the QHA EOS?"

ISMEAR=2 is dropped here: the sigma/s2 results already showed s2-s1 <= 0.2
meV/atom range, so s2 is not where the uncertainty is.

Writes only under --out/k15.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from build_align_inputs import A0, FRAMES, SELECTED, STRAIN, b2_cell, write_case

K15 = 15
TAGS = {"s1": (1, 0.2), "sm5": (-5, None)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    out = args.out / "k15"
    out.mkdir(parents=True, exist_ok=True)

    eps_l = np.linspace(STRAIN[0], STRAIN[1], FRAMES)
    v0 = A0 ** 3
    cases = []
    for idx in SELECTED:
        eps = float(eps_l[idx - 1])
        a = A0 * (1.0 + eps)
        st = b2_cell(a)
        meta = {"block": "k15", "frame": idx, "eps": eps, "a_angstrom": a,
                "volume_angstrom3": float(st.volume),
                "volume_scale_vs_v0": float(st.volume) / v0}
        for tag, (ism, sig) in TAGS.items():
            cases.append(write_case(out / f"f{idx:02d}_k{K15}_{tag}", st, K15,
                                    f"align-k{K15}-f{idx:02d}-{tag}", ism, sig, meta))

    for c in cases:
        mark = "OK " if c["ok"] else "BAD"
        print(f"  {mark} {c['dir']:<46} V={c['volume_angstrom3']:.4f} "
              f"(x{c['volume_scale_vs_v0']:.5f}) k={c['kgrid']} ISMEAR={c['ismear']:>2}")
        for p in c["problems"]:
            print(f"        !! {p}")
    (out / "validation.json").write_text(json.dumps(cases, indent=2) + "\n")
    bad = [c for c in cases if not c["ok"]]
    print(f"\n{len(cases) - len(bad)}/{len(cases)} validated")
    if bad:
        raise SystemExit("refusing to submit")


if __name__ == "__main__":
    main()

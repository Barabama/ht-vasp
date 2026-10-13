#!/usr/bin/env python
"""Stage 0A follow-ups, decided by what Stage 0A actually found.

Three questions the first pass could not settle:

1. Which Gamma-centred family converges faster?  Stage 0A compared them at
   different reciprocal spacings, so the comparison was inconclusive.  This adds
   denser points in BOTH families on the same AlTa cell:
       family A (even n, shift 0.5): n = 20, 24
       family B (odd  n, shift 0  ): n = 19, 21
   (n = 12 / 16 for A and 11..17 for B already exist.)

2. Is the SER-Al-1 dE(12 -> 16) = 6.26 meV/atom real convergence, or is the
   mesh still far from converged?  Adds n = 24, 28 (shift 0.5).

3. Is the Cr magnetism mesh-dependent?  Stage 0A screened spins only at n = 16.
   Adds ISPIN=1 and ISPIN=2/Cr-start at n = 12.

Writes only under --out/stage0 (refuses to overwrite an existing INCAR).
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import build_stage0 as S

A_ALTA = S.A_ALTA


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    root = args.out / "stage0"

    reps = S.build_structures()
    cases = []

    # ---- 1. denser points in both Gamma-centred families, AlTa cell --------
    st = reps["BCC-Al-Ta-2"]["struct"]
    for n, shift, fam in ((20, 0.5, "A/even+shift0.5"), (24, 0.5, "A/even+shift0.5"),
                          (19, 0.0, "B/odd+shift0"), (21, 0.0, "B/odd+shift0")):
        cases.append(S.write_case(
            root / "BCC-Al-Ta-2" / f"family{n:02d}_s{int(shift * 2)}_s1", st, n, shift,
            f"stage0-AlTa-family-{fam}-n{n}", 1, 0.2, S.magmom_for(st), 2,
            {"endpoint": "BCC-Al-Ta-2", "class": f"family probe {fam}",
             "a_angstrom": A_ALTA, "a_source": "same cell as Stage 0A",
             "density_target": None,
             "dk_angstrom_inv": 2 * math.pi / (A_ALTA * n),
             "volume_angstrom3": float(st.volume),
             "density_invariant_n3V": n ** 3 * float(st.volume)}))

    # ---- 2. SER-Al-1: push the mesh up ------------------------------------
    st = reps["SER-Al-1"]["struct"]
    for n in (24, 28):
        cases.append(S.write_case(
            root / "SER-Al-1" / f"dense_k{n:02d}_s1", st, n, 0.5,
            f"stage0-SER-Al-1-dense-k{n}", 1, 0.2, S.magmom_for(st), 2,
            {"endpoint": "SER-Al-1", "class": "density follow-up",
             "a_angstrom": reps["SER-Al-1"]["a"],
             "a_source": reps["SER-Al-1"]["source"],
             "density_target": None,
             "dk_angstrom_inv": 2 * math.pi / (reps["SER-Al-1"]["a"] * n),
             "volume_angstrom3": float(st.volume),
             "density_invariant_n3V": n ** 3 * float(st.volume)}))

    # ---- 3. Cr magnetism at the coarser mesh ------------------------------
    st = reps["BCC-Al-Cr-2"]["struct"]
    n = S.mesh_for(float(st.volume), 12)
    base = {"endpoint": "BCC-Al-Cr-2", "class": "magnetic screening @d12",
            "a_angstrom": reps["BCC-Al-Cr-2"]["a"],
            "a_source": reps["BCC-Al-Cr-2"]["source"],
            "density_target": 12,
            "dk_angstrom_inv": 2 * math.pi / (reps["BCC-Al-Cr-2"]["a"] * n),
            "volume_angstrom3": float(st.volume)}
    for tag, ispin, mm, desc in (
            ("ispin1", 1, "MAGMOM = 0.0 0.0", "ISPIN=1 (no spin)"),
            ("ispin2_cr", 2, S.magmom_for(st), "ISPIN=2, Al 1.0 / Cr 3.0")):
        cases.append(S.write_case(
            root / "BCC-Al-Cr-2" / f"mag_{tag}_k{n:02d}", st, n, 0.5,
            f"stage0-AlCr-{tag}-d12", 1, 0.2, mm, ispin,
            {**base, "magmom_policy": desc}))

    bad = [c for c in cases if not c["ok"]]
    for c in cases:
        mark = "OK " if c["ok"] else "BAD"
        print(f"  {mark} {c['dir'].split('stage0/')[-1]:<28} n={c['kgrid']:>2} "
              f"shift={c['shift_value']} V={c['volume_angstrom3']:8.3f} "
              f"n3V={c.get('density_invariant_n3V', 0):>8.0f} ISMEAR={c['ismear']}")
        for p in c["problems"]:
            print(f"        !! {p}")
    print(f"\n{len(cases) - len(bad)}/{len(cases)} pass the input contract")
    if bad:
        raise SystemExit("refusing to submit")


if __name__ == "__main__":
    main()

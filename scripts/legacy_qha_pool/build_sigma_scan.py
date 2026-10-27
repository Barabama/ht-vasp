#!/usr/bin/env python
"""SIGMA scan: which smearing parameter can serve both E(V) and the phonon forces?

Stage 0A left one question open (report §4.9.4): the residual k-mesh scatter in
AlTa is 1.6-3.5 meV/atom with ISMEAR=1/SIGMA=0.2, and tetrahedron gives only
0.39.  Tetrahedron is ruled out for the force path (9% force error, and it does
not respect the folding equivalence).  So the question becomes: does a SMALLER
SIGMA remove the extrapolation noise while keeping the SCF and the forces sane?

Minimal experiment, per the review:
  * 4 representatives -- the three that matter plus a stable baseline
        BCC-Al-Ta-2   finite-q soft mode, most sensitive in Stage 0A
        SER-Al-1      largest dE(12->16) in Stage 0A (6.26 meV/atom)
        BCC-Al-Cr-2   magnetic branch
        BCC-Ta-Ta-2   stable non-magnetic baseline
  * 5 volume points from the PRODUCTION QHA grid, so the result applies to the
    E(V) actually fitted:
        QhaWorker(linear_strain=(-0.05, 0.05), number_of_frames=8)
        frames 1, 3, 5, 6, 8  ->  eps = -0.0500, -0.0214, +0.0071, +0.0214, +0.0500
  * 5 smearings
        s1_005  ISMEAR=1 SIGMA=0.05   production candidate
        s1_010  ISMEAR=1 SIGMA=0.10   production candidate
        s1_020  ISMEAR=1 SIGMA=0.20   current production (reference)
        s2_010  ISMEAR=2 SIGMA=0.10   candidate B
        sm5     ISMEAR=-5             static audit only

Held fixed across every case (the review's requirement):
    same POSCAR per volume point, same explicit KPOINTS per endpoint, same
    POTCAR, same ISPIN/MAGMOM, same ENCUT, same ALGO.

The k mesh is FIXED per endpoint (computed at that endpoint's equilibrium
volume, d = 16 in AlTa-equivalent units, even n + shift 0.5).  It deliberately
does NOT scale with volume -- report §6.1 requires that volume changes must not
silently change the mesh, and a fixed mesh is also the conventional QHA choice.

AlCr starts from the SAME magnetic branch (Cr 3.0) at every volume and every
sigma: Stage 0A showed a zero starting moment converges to a non-magnetic
metastable solution 283.7 meV/atom higher (report §4.9.3).

Writes only under --out.  Refuses to overwrite an existing INCAR.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import build_stage0 as S

# production QHA volume grid
FRAMES = 8
STRAIN = (-0.05, 0.05)
VOLUME_FRAMES = (1, 3, 5, 6, 8)        # 1-indexed
DENSITY = 16                           # AlTa-equivalent k density

SMEARINGS = {
    "s1_005": (1, 0.05),
    "s1_010": (1, 0.10),
    "s1_020": (1, 0.20),
    "s2_010": (2, 0.10),
    "sm5": (-5, None),
}
REFERENCE = "s1_020"


def scaled(struct, eps: float):
    """Uniformly scale a cubic cell by (1 + eps)."""
    lat = struct.lattice.matrix * (1.0 + eps)
    from pymatgen.core import Lattice, Structure
    return Structure(Lattice(lat), [s.specie for s in struct],
                     [s.frac_coords for s in struct], coords_are_cartesian=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    root = args.out / "sigma"

    reps = S.build_structures()
    targets = ("BCC-Al-Ta-2", "SER-Al-1", "BCC-Al-Cr-2", "BCC-Ta-Ta-2")
    eps_l = np.linspace(STRAIN[0], STRAIN[1], FRAMES)

    cases = []
    for name in targets:
        spec = reps[name]
        base = spec["struct"]
        # mesh fixed at the equilibrium volume, then held for every volume point
        n = S.mesh_for(float(base.volume), DENSITY)
        for fr in VOLUME_FRAMES:
            eps = float(eps_l[fr - 1])
            st = scaled(base, eps)
            meta = {
                "endpoint": name, "class": spec["class"],
                "a_source": spec["source"],
                "a_angstrom": spec["a"] * (1.0 + eps),
                "a0_angstrom": spec["a"],
                "frame": fr, "eps": eps,
                "volume_angstrom3": float(st.volume),
                "volume_scale_vs_v0": float(st.volume) / float(base.volume),
                "density_target": DENSITY,
                "density_invariant_n3V": n ** 3 * float(base.volume),
                "kpoints_policy": ("fixed per endpoint at the equilibrium volume; "
                                   "does not scale with volume"),
                "block": "sigma",
            }
            for tag, (ism, sig) in SMEARINGS.items():
                cases.append(S.write_case(
                    root / name / f"f{fr:02d}_{tag}", st, n, 0.5,
                    f"sigma-{name}-f{fr:02d}-{tag}", ism, sig,
                    S.magmom_for(st), 2, meta))

    bad = [c for c in cases if not c["ok"]]
    report = {
        "schema": "legacy-qha-pool/sigma-scan/v1",
        "purpose": "choose the SIGMA that keeps E(V) smooth AND forces sane",
        "volume_grid": {
            "source": f"QhaWorker linear_strain={STRAIN} number_of_frames={FRAMES}",
            "frames_1indexed": list(VOLUME_FRAMES),
            "eps": [float(eps_l[f - 1]) for f in VOLUME_FRAMES],
        },
        "smearings": {k: {"ISMEAR": v[0], "SIGMA": v[1]} for k, v in SMEARINGS.items()},
        "reference_smearing": REFERENCE,
        "density_target": DENSITY,
        "kpoints_family": "even n, shift 0.5 0.5 0.5 (explicit folding-compatible)",
        "recipe": {"ENCUT": S.ENCUT, "EDIFF": S.EDIFF, "ISYM": S.ISYM,
                   "NCORE": S.NCORE, "NELM": S.NELM, "NELMIN": S.NELMIN,
                   "ALGO": S.ALGO, "LORBIT": S.LORBIT,
                   "potcar_functional": S.POTCAR_FUNCTIONAL,
                   "magmom_table": S.MAGMOM,
                   "alcr_magnetic_policy": "MAGMOM Al 1.0 / Cr 3.0 at every volume "
                                           "and every sigma (report 4.9.3)"},
        "n_cases": len(cases), "n_bad": len(bad), "cases": cases,
    }
    (root / "validation.json").write_text(json.dumps(report, indent=2) + "\n")

    print(f"{'endpoint':<16} {'n':>3} {'frame':>5} {'eps':>9} {'V/V0':>8} "
          f"{'V/A^3':>9}  smearings")
    print("-" * 88)
    shown = set()
    for c in cases:
        key = (c["endpoint"], c["frame"])
        if key in shown:
            continue
        shown.add(key)
        print(f"{c['endpoint']:<16} {c['kgrid']:>3} {c['frame']:>5} {c['eps']:>9.5f} "
              f"{c['volume_scale_vs_v0']:>8.5f} {c['volume_angstrom3']:>9.4f}  "
              f"{len(SMEARINGS)}")
    print(f"\n{len(cases) - len(bad)}/{len(cases)} cases pass the input contract")
    for c in bad:
        print(f"  !! {c['dir']}")
        for p in c["problems"]:
            print(f"       {p}")
    if bad:
        raise SystemExit("refusing to submit: some inputs failed validation")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Extract the soft-mode eigenvector from ALREADY-COMPUTED force constants.

Why this exists
---------------
The pool's dominant failure mode is a finite-q soft mode (report 4.2b / 7.6b).
Diagnosing it needs the eigenVECTOR, not the frequency: the frequency says how
deep the instability is, the eigenvector says which way to distort.  Frequency
alone cannot answer "does this point at a low-symmetry phase, a double well, or
a numerical artefact".

phonopy computes eigenvectors for free once the force constants exist -- they
are a property of the dynamical matrix, not of VASP.  So for any case whose
force constants were kept, the eigenvector can be recovered with ZERO new DFT.

Two honest limits, both of which the report states:

* This only works where the FORCE CONSTANTS survive.  It cannot be reconstructed
  from a thermal_properties YAML, so the legacy pool (which has only the YAML)
  cannot be retrofitted -- only cases we computed and kept the fc for.
* "Zero DFT" is not "free": it is a post-processing step that has to be wired
  into the pipeline BEFORE the thermal-properties stage, or the fc gets thrown
  away and the eigenvector with it.

What it writes, per case: the lowest-frequency q-point on a mesh, its frequency,
the branch index, and the mass-weighted eigenvector as phonopy defines it (plus
the atom ordering so the vector can be mapped back to the supercell), together
with the hashes needed to tie it back to the run.

Usage
-----
    python extract_soft_mode.py --fcdiff /tmp/fcdiff/stage0b \
                                --root /nfs_ssd/tmp/legacy-qha-pool/stage0b \
                                --json soft_modes.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import analyze_stage0b as A
import build_stage0b as B0
from htvasp.dfpt_qha.phonons import degeneracy, eigenvector

ENDPOINTS = ("BCC-Ta-Ta-2", "BCC-Al-Ta-2", "BCC-Al-Cr-2")
MESH = 12          # same mesh analyze_stage0b used for mesh_min / mesh_n_negative
TOL_THz = 1e-4     # what counts as imaginary, same as analyze_stage0b


def fc_hash(fc: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(fc, dtype=np.float64).tobytes()
                          ).hexdigest()[:16]


def lowest_mode(endpoint: str, fc: np.ndarray) -> dict:
    """Lowest-frequency mode on the mesh, with its eigenvector.

    THREE traps this function has to survive, all of which produce
    plausible-looking output when mishandled:

    1. BASIS.  phonopy's mesh runs on the PRIMITIVE cell, so `eigenvectors` has
       shape (nq, 3*n_prim, 3*n_prim): for the 2-atom B2 cell that is 6, NOT the
       16-atom supercell's 48.  Indexing a 6-vector with supercell symbols still
       looks reasonable (supercell atoms 0..7 are all Al, so "dominant atom = Al"
       comes out right by accident).  Same failure shape as build_stage0b's
       force/POSCAR ordering trap.

    2. ROW vs COLUMN.  `np.linalg.eigh` puts eigenvectors in the COLUMNS and
       phonopy stores that matrix unchanged, so the eigenvector of band j is
       `eigenvectors[iq][:, j]` -- NOT `[iq][j, :]`.  Both are unit vectors, so
       the wrong one is invisible in the norm and only shows up as a wrong
       direction.  `htvasp.dfpt_qha.phonons.eigenvector()` reads the column AND
       asserts |D v - lambda v| ~ 0, so this cannot pass silently.

    3. DEGENERACY.  On a degenerate branch `eigh` returns an ARBITRARY basis of
       the subspace, so a single eigenvector has no meaning.  The degeneracy is
       recorded here and the caller is expected to treat such a q as a subspace
       (see analyze_soft_mode.py).

    Also: this is the MINIMUM ON A SAMPLING MESH, which is NOT the instability
    wavevector -- the mesh may not even contain the symmetry point.  See
    analyze_soft_mode.py; do not use this q for mode following.
    """
    ph, _ = A.make_phonopy(endpoint)
    ph.force_constants = fc
    ph.run_mesh([MESH] * 3)
    d = ph.get_mesh_dict()
    freqs = np.array(d["frequencies"])           # (nq, nbands), nbands = 3*n_prim
    qpoints = np.array(d["qpoints"])
    iq, ib = np.unravel_index(np.argmin(freqs), freqs.shape)
    q = [float(x) for x in qpoints[iq]]

    pr = ph.primitive
    be = eigenvector(ph, q, int(ib))   # checked column read
    vec, f = be.vector, be.frequency_thz
    assert vec.shape[0] == 3 * len(pr), (vec.shape, len(pr))
    deg = degeneracy(ph, q, int(ib))
    per_atom = np.array([np.linalg.norm(vec[3 * a:3 * a + 3])
                         for a in range(len(pr))])
    return {
        "mesh": [MESH] * 3,
        "qpoint": q,
        "qpoint_distance_to_Gamma": float(np.linalg.norm(qpoints[iq])),
        "band_index": int(ib),
        "frequency_THz": float(f),
        "is_imaginary": bool(f < -TOL_THz),
        "n_imaginary_on_mesh": int((freqs < -TOL_THz).sum()),
        "n_modes_on_mesh": int(freqs.size),
        "degenerate_bands": deg,
        "is_degenerate": len(deg) > 1,
        "eigen_residual": be.residual,
        "row_reading_residual": be.row_reading_residual,
        "eigenvector_basis": ("primitive cell (3*n_prim components), NOT the "
                              "supercell; map to the supercell via p2s_map. "
                              "Read from the COLUMN of phonopy's eigenvector "
                              "matrix; the residual above is the proof."),
        "eigenvector_convention": ("phonopy mass-weighted eigenvector of the "
                                   "dynamical matrix at this q-point"),
        "primitive_symbols": [str(s) for s in pr.symbols],
        "primitive_scaled_positions": [[float(x) for x in p]
                                       for p in pr.scaled_positions],
        "p2s_map": [int(x) for x in pr.p2s_map],
        "eigenvector": [complex(z).real if abs(complex(z).imag) < 1e-10
                        else {"re": complex(z).real, "im": complex(z).imag}
                        for z in vec],
        "eigenvector_abs": [float(abs(z)) for z in vec],
        "per_primitive_atom_amplitude": [float(x) for x in per_atom],
        "dominant_primitive_atom": int(np.argmax(per_atom)),
        "dominant_primitive_symbol": str(pr.symbols[int(np.argmax(per_atom))]),
        "primitive_species_weight": {
            s: float(sum(per_atom[a] ** 2 for a in range(len(pr))
                         if str(pr.symbols[a]) == s))
            for s in dict.fromkeys(str(x) for x in pr.symbols)},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fcdiff", type=Path, action="append", required=True,
                    help="root holding <endpoint>/<smearing>/fc.npy.  Repeat it to "
                         "compare densities; the first is the reference.")
    ap.add_argument("--label", action="append", default=None,
                    help="label per --fcdiff, same order (default: directory name)")
    ap.add_argument("--root", type=Path, default=None,
                    help="VASP tree, only to record whether the case dirs survive")
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    labels = args.label or [p.name for p in args.fcdiff]
    if len(labels) != len(args.fcdiff):
        raise SystemExit("--label count must match --fcdiff count")

    out: dict = {}
    for lab, root_fc in zip(labels, args.fcdiff):
        for ep in ENDPOINTS:
            for tag in A.SMEARINGS:
                key = f"{lab}|{ep}|{tag}"
                f = root_fc / ep / tag / "fc.npy"
                if not f.is_file():
                    out[key] = {"error": "no fc.npy"}
                    continue
                fc = np.load(f)
                rec = lowest_mode(ep, fc)
                rec["fcdiff_label"] = lab
                rec["fc_file"] = str(f)
                rec["fc_sha256_16"] = fc_hash(fc)
                rec["fc_shape"] = list(fc.shape)
                rec["n_supercell_atoms"] = int(fc.shape[0])
                if args.root:
                    d = args.root / ep
                    rec["vasp_case_dirs_present"] = sorted(
                        p.name for p in d.glob(f"a0.01_d*_{tag}") if p.is_dir())
                out[key] = rec

    print(f"{'endpoint':<15} {'smearing':>8} {'q*':>22} {'f* /THz':>9} "
          f"{'band':>5} {'imag':>5} {'degen':>6} {'dominant':>10}")
    print("-" * 92)
    for lab in labels:
        print(f"\n########## {lab} ##########")
        for ep in ENDPOINTS:
            for tag in A.SMEARINGS:
                r = out[f"{lab}|{ep}|{tag}"]
                if "error" in r:
                    print(f"{ep:<15} {A.LABEL[tag]:>8}   {r['error']}")
                    continue
                q = r["qpoint"]
                print(f"{ep:<15} {A.LABEL[tag]:>8} "
                      f"({q[0]:>6.4f},{q[1]:>6.4f},{q[2]:>6.4f}) "
                      f"{r['frequency_THz']:>9.4f} {r['band_index']:>5} "
                      f"{str(r['is_imaginary']):>5} "
                      f"{len(r['degenerate_bands']):>6} "
                      f"{r['dominant_primitive_symbol']:>10}")

        print(f"\n=== 本征矢的物种权重（primitive 基，参照展宽 {A.LABEL[A.REF]}）===")
        for ep in ENDPOINTS:
            r = out[f"{lab}|{ep}|{A.REF}"]
            if "error" in r:
                continue
            w = r["primitive_species_weight"]
            tot = sum(w.values()) or 1.0
            s = "  ".join(f"{k}={100 * v / tot:.1f}%" for k, v in
                          sorted(w.items(), key=lambda x: -x[1]))
            print(f"  {ep:<15} f*={r['frequency_THz']:+.4f} THz @ "
                  f"{[round(x, 4) for x in r['qpoint']]}  {s}")

    if len(labels) > 1:
        print("\n=== 跨密度的本征矢一致性（同一个软模是否同样指向）===")
        print(f"{'endpoint':<15} {'smearing':>8} {'q* 是否相同':>12} "
              f"{'Δf* /THz':>10} {'物种权重是否相同':>18}")
        print("-" * 70)
        for ep in ENDPOINTS:
            for tag in A.SMEARINGS:
                rs = [out[f"{lab}|{ep}|{tag}"] for lab in labels]
                if any("error" in r for r in rs):
                    continue
                qs = [tuple(round(x, 6) for x in r["qpoint"]) for r in rs]
                same_q = len(set(qs)) == 1
                df = rs[-1]["frequency_THz"] - rs[0]["frequency_THz"]
                # compare the species weight ordering, not the raw numbers
                orders = []
                for r in rs:
                    w = r["primitive_species_weight"]
                    tot = sum(w.values()) or 1.0
                    orders.append(tuple(k for k, _ in
                                        sorted(w.items(), key=lambda x: -x[1] / tot)))
                print(f"{ep:<15} {A.LABEL[tag]:>8} {str(same_q):>12} "
                      f"{df:>10.4f} {str(len(set(orders)) == 1):>18}")

    if args.json:
        args.json.write_text(json.dumps(out, indent=2, default=str) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()

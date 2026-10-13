#!/usr/bin/env python
"""Stage 0B at the SAME reciprocal density as the E(V) path (d = 16, not 12).

Why this exists
---------------
Stage 0B (build_stage0b.py) ran its phonon supercells at mesh 6, i.e. density 12
in AlTa-equivalent units, while the E(V) work (Stage 0A, the SIGMA scan) ran the
2-atom cell at density 16.  The comment in build_stage0b.py claimed the two were
commensurate; that was wrong -- the density-16 partner of a 16-atom supercell is
mesh 8, not 6.

That matters because report 4.6.5/4.9.1 measured a ~3.5 meV/atom change in the
static energy of BCC-Al-Ta-2 between density 12 and 16, and report 4.11.7
measured a ~4% change in the N-point soft-mode depth between two meshes.  So
Stage 0B's phonon absolute values are at a coarser density than the E(V) they
would be combined with in a QHA.

This batch re-runs the same 40 cases at supercell mesh 8 (parent 16) so the two
can be compared at matched density.  Nothing is overwritten: the mesh-6 tree
stays under stage0b/, this one goes to stage0b-d16/.

Everything else is identical to build_stage0b -- same cells, same displacements,
same four smearings, same AlCr magnetic branch, explicit folding-compatible
KPOINTS (16-atom k=8 shift 0.5 <-> 2-atom k=16 shift 0.5).

Writes only under --out.  Refuses to overwrite an existing INCAR.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import build_stage0 as S0
import build_stage0b as B0

SUPER_MESH = 8          # density 16 partner of the 2-atom mesh 16
PARENT_MESH = 16
AMPLITUDE = 0.01        # 0.03 is answered by the mesh-6 batch (report 4.11.5)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    root = args.out / "stage0b-d16"

    reps = S0.build_structures()
    cases: list[dict] = []
    manifest: dict = {}

    for name in B0.TARGETS:
        base = reps[name]["struct"]
        # cross-check the density claim instead of asserting it in a comment
        n_ev = S0.mesh_for(float(base.volume), 16)
        assert n_ev == PARENT_MESH, (name, n_ev, PARENT_MESH)
        assert SUPER_MESH * 2 == PARENT_MESH
        ph, sup = B0.phonopy_supercell(base, B0.SUPERCELL)
        assert len(sup) == B0.SUPER_NATOMS
        fp = B0.cell_fingerprint(sup)
        disps = B0.displacements_for(ph, AMPLITUDE)
        manifest[name] = {
            "n_supercell_atoms": len(sup),
            "supercell_mesh": SUPER_MESH, "parent_2atom_mesh": PARENT_MESH,
            "folding_verified": True,
            "cell_fingerprint": fp,
            "n_displacements": len(disps),
            "displacements": disps,
            "density_note": (f"2-atom mesh {n_ev} (d=16) <-> 16-atom mesh "
                             f"{SUPER_MESH}; matches the E(V) density"),
        }
        for disp in disps:
            dstruct = B0.displaced(sup, disp["atom"], disp["vector"])
            for tag, (ism, sig) in B0.SMEARINGS.items():
                meta = {
                    "endpoint": name, "class": reps[name]["class"],
                    "block": "stage0b-d16",
                    "supercell": f"{B0.SUPERCELL}x{B0.SUPERCELL}x{B0.SUPERCELL}",
                    "n_supercell_atoms": len(sup),
                    "cell_fingerprint": fp,
                    "atom_order_source": "phonopy supercell order",
                    "amplitude_angstrom": AMPLITUDE,
                    "disp_index": disp["index"], "disp_atom": disp["atom"],
                    "disp_atom_species": str(sup[disp["atom"]].specie),
                    "disp_vector": disp["vector"],
                    "grid": [SUPER_MESH] * 3, "shift": [0.5] * 3,
                    "cell_atoms": 2,
                    "folding_parent_grid": [PARENT_MESH] * 3,
                    "folding_verified": True,
                    "a_supercell_angstrom": sup.lattice.a,
                    "volume_angstrom3": float(sup.volume),
                }
                cases.append(B0.write_case(
                    root / name / f"a{AMPLITUDE:g}_d{disp['index']:02d}_{tag}",
                    dstruct, SUPER_MESH, 0.5,
                    f"stage0b-d16-{name}-d{disp['index']}-{tag}",
                    ism, sig, S0.magmom_for(sup), meta))

    bad = [c for c in cases if not c["ok"]]
    report = {
        "schema": "legacy-qha-pool/stage0b-d16/v1",
        "purpose": ("re-run Stage 0B at the E(V) density (d=16 / supercell mesh 8) "
                    "because the mesh-6 batch was at d=12"),
        "smearings": {k: {"ISMEAR": v[0], "SIGMA": v[1]}
                      for k, v in B0.SMEARINGS.items()},
        "amplitude_angstrom": AMPLITUDE,
        "kpoints_policy": ("explicit folding-compatible: 16-atom k=8 shift 0.5 "
                           "<-> 2-atom k=16 shift 0.5"),
        "manifest": manifest, "n_cases": len(cases), "n_bad": len(bad),
        "cases": cases,
    }
    (root / "validation.json").write_text(json.dumps(report, indent=2) + "\n")

    print(f"{'endpoint':<15} {'#disp':>6} {'atoms':>6} {'mesh':>5} {'parent':>7} "
          f"{'cases':>6}")
    for name in B0.TARGETS:
        m = manifest[name]
        print(f"{name:<15} {m['n_displacements']:>6} {m['n_supercell_atoms']:>6} "
              f"{m['supercell_mesh']:>5} {m['parent_2atom_mesh']:>7} "
              f"{m['n_displacements'] * len(B0.SMEARINGS):>6}")
    print(f"\n{len(cases) - len(bad)}/{len(cases)} cases pass the input contract")
    for c in bad:
        print(f"  !! {c['dir']}")
        for p in c["problems"]:
            print(f"       {p}")
    if bad:
        raise SystemExit("refusing to submit")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Stage 0C, step C1+C2: FCC/L1_2 and SER inputs, plus the input-contract gate.

Two endpoints, deliberately the two classes Stage 0A/0B never touched:

    FCC-Zr-Al-4   L1_2, 4-atom conventional cell, Zr on 1a / Al on 3c
    SER-Al-1      fcc primitive, 1 atom

Each gets:
  * one equilibrium-volume static (C2's single reference volume) at the frozen
    recipe ISMEAR=1 / SIGMA=0.20, and
  * one phonon-supercell run per symmetry-reduced displacement (C3's object),
    so the supercell ordering can be checked before any real phonon job is
    submitted.

The folding-compatible KPOINTS rule generalises: for a supercell that is an
m x m x m repeat of the input cell, equal reciprocal spacing needs
k_super = k_cell / m, and the meshes sample the same folded points when the
Gamma parity matches (report 4.8.1).  For FCC the input cell is the 4-atom
conventional cell, so its 2x2x2 supercell holds 32 atoms -- larger than the
16-atom B2 supercell, and recorded as such.

Writes only under --out.  Refuses to overwrite an existing INCAR.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp.inputs import Potcar

import build_stage0 as S0
import build_stage0b as B0

SMEARING = (1, 0.20)          # frozen Stage 0 provisional recipe
DENSITY = 16
SUPERCELL = 2
TARGETS = ("FCC-Zr-Al-4", "SER-Al-1")


def write_case(d: Path, st: Structure, mesh: int, shift: float, system: str,
               ismear: int, sigma: float | None, magmom: str, meta: dict) -> dict:
    if (d / "INCAR").exists():
        raise SystemExit(f"refusing to overwrite existing case: {d}")
    d.mkdir(parents=True, exist_ok=True)
    Poscar(st).write_file(str(d / "POSCAR"))
    (d / "INCAR").write_text(S0.incar_text(system, ismear, sigma, magmom))
    (d / "KPOINTS").write_text(S0.kpoints_text(mesh, shift))
    species_order = list(dict.fromkeys(str(s.specie) for s in st))
    Potcar(species_order, functional=S0.POTCAR_FUNCTIONAL).write_file(str(d / "POTCAR"))
    from htvasp.dfpt_qha.gate import validate_dir
    gate = validate_dir(d)
    record = {**gate.model_dump(mode="json"), "kgrid": mesh, "shift_value": shift,
              **meta}
    (d / "case.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    root = args.out / "stage0c" / "inputs"

    reps = S0.build_structures()
    cases: list[dict] = []
    manifest: dict = {}

    for name in TARGETS:
        base = reps[name]["struct"]
        # --- equilibrium-volume static (C2) ------------------------------
        n = S0.mesh_for(float(base.volume), DENSITY)
        shift = 0.5 if n % 2 == 0 else 0.0
        manifest[name] = {
            "class": reps[name]["class"], "a_source": reps[name]["source"],
            "n_cell_atoms": len(base),
            "cell_fingerprint": B0.cell_fingerprint(base),
            "mesh": n, "shift": shift,
            "contains_gamma": (n % 2 == 1) if shift == 0 else (n % 2 == 0),
            "n3V": n ** 3 * float(base.volume),
            "supercell_n_atoms": None,
        }
        cases.append(write_case(
            root / name / "static_eq", base, n, shift,
            f"stage0c-{name}-static", SMEARING[0], SMEARING[1],
            S0.magmom_for(base),
            {"endpoint": name, "block": "stage0c-static",
             "class": reps[name]["class"], "a_source": reps[name]["source"],
             "n_cell_atoms": len(base), "volume_angstrom3": float(base.volume),
             "cell_fingerprint": B0.cell_fingerprint(base),
             "role": "C2 single reference volume, verifies VASP accepts the input"}))

        # --- phonon supercell displacements (C3's object) ----------------
        ph, sup = B0.phonopy_supercell(base, SUPERCELL)
        manifest[name]["supercell_n_atoms"] = len(sup)
        manifest[name]["supercell_mesh"] = n // 2 if n % 2 == 0 else None
        manifest[name]["supercell_shift"] = (
            [0.0] * 3 if (n // 2) % 2 else [0.5] * 3) if n % 2 == 0 else None
        manifest[name]["supercell_fingerprint"] = B0.cell_fingerprint(sup)
        amp = 0.01
        ph.generate_displacements(distance=amp, is_plusminus="auto", is_diagonal=True)
        disps = [{"index": i + 1, "atom": int(d[0]),
                  "species": str(sup[int(d[0])].specie),
                  "vector": [float(d[1]), float(d[2]), float(d[3])]}
                 for i, d in enumerate(ph.displacements)]
        manifest[name]["n_displacements"] = len(disps)
        manifest[name]["displacements"] = disps
        sm = n // 2 if n % 2 == 0 else n          # supercell mesh
        ss = 0.0 if sm % 2 else 0.5
        for disp in disps:
            dstruct = B0.displaced(sup, disp["atom"], disp["vector"])
            cases.append(write_case(
                root / name / f"phon_d{disp['index']:02d}", dstruct, sm, ss,
                f"stage0c-{name}-phon-d{disp['index']}", SMEARING[0], SMEARING[1],
                S0.magmom_for(sup),
                {"endpoint": name, "block": "stage0c-phonon",
                 "class": reps[name]["class"],
                 "n_supercell_atoms": len(sup),
                 "cell_fingerprint": B0.cell_fingerprint(sup),
                 "amplitude_angstrom": amp,
                 "disp_index": disp["index"], "disp_atom": disp["atom"],
                 "disp_atom_species": disp["species"],
                 "disp_vector": disp["vector"],
                 "grid": [sm] * 3, "shift": [ss] * 3,
                 "cell_atoms": len(base),
                 "folding_parent_grid": [n] * 3,
                 "folding_verified": n % 2 == 0,
                 "volume_angstrom3": float(sup.volume)}))

    bad = [c for c in cases if not c["ok"]]
    report = {
        "schema": "legacy-qha-pool/stage0c-inputs/v1",
        "purpose": "C1/C2 inputs for FCC/L1_2 and SER: static + phonon supercell",
        "recipe": {"ISMEAR": SMEARING[0], "SIGMA": SMEARING[1],
                   "ENCUT": S0.ENCUT, "EDIFF": S0.EDIFF, "ISYM": S0.ISYM,
                   "NCORE": S0.NCORE, "ALGO": S0.ALGO, "LORBIT": S0.LORBIT,
                   "potcar_functional": S0.POTCAR_FUNCTIONAL,
                   "displacement_angstrom": 0.01,
                   "status": "provisional_stage0 (tested on BCC/B2 representatives)"},
        "manifest": manifest, "n_cases": len(cases), "n_bad": len(bad),
        "cases": cases,
    }
    (root / "validation.json").write_text(json.dumps(report, indent=2) + "\n")

    for name in TARGETS:
        m = manifest[name]
        print(f"{name:<14} cell={m['n_cell_atoms']:>2} atoms  mesh={m['mesh']:>2} "
              f"shift={m['shift']} gamma={m['contains_gamma']}  "
              f"supercell={m['supercell_n_atoms']} atoms mesh={m['supercell_mesh']} "
              f"shift={m['supercell_shift']}  disp={m['n_displacements']}")
        for d in m["displacements"]:
            print(f"    d{d['index']}: atom {d['atom']:>2} ({d['species']:<3}) "
                  f"vec={[round(v, 5) for v in d['vector']]}")
    print(f"\n{len(cases) - len(bad)}/{len(cases)} cases pass the input contract")
    for c in bad:
        print(f"  !! {c['dir']}")
        for p in c["problems"]:
            print(f"       {p}")
    if bad:
        raise SystemExit("refusing to submit: some inputs failed validation")


if __name__ == "__main__":
    main()

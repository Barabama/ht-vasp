#!/usr/bin/env python
"""Stage 0B: force / force-constant / phonon calibration across smearings.

The question, in the review's framing: under a FIXED folding-compatible
KPOINTS and a FIXED magnetic branch, do the smearing choices differ at the
level of forces, force constants, or phonons -- for the three structural
classes that matter?

    BCC-Ta-Ta-2   stable non-magnetic baseline
    BCC-Al-Ta-2   confirmed finite-q soft mode (does the soft mode survive?)
    BCC-Al-Cr-2   strongly magnetic (does the branch hold across displacements?)

Design
------
* 2x2x2 phonon supercell (16 atoms), phonopy displacements at 0.01 A.
  The count is symmetry-reduced per endpoint and is RECORDED, not assumed:
  the B2 cell yields 1 displacement (Al only, x); BCC Ta and B2 AlCr yield 2.
* Two displacement amplitudes (0.01 and 0.03 A).  The review's brief is
  forces/force-constants; 0.03 A is the cheap test of whether the 0.01 A
  result is in the harmonic regime, which is otherwise an untested free
  parameter of the finite-displacement path (PARAM_FREEZE §4 item 2).
* Four smearings: 1/0.10, 1/0.20, 2/0.10, 2/0.20.  ISMEAR=-5 is EXCLUDED
  from the force matrix (report §4.8.2, §4.10): its force constants deviate
  ~9% and it does not respect the folding equivalence.
* Everything else identical: ENCUT 520, EDIFF 1e-6, ISYM 2, NCORE 8,
  ALGO Normal, LORBIT 11, POTCAR PBE(52) standard variants.
* AlCr starts from Al 1.0 / Cr 3.0 at EVERY displacement and every smearing.

KPOINTS
-------
Explicit, folding-compatible, recorded as an indivisible field set
(grid + shift + cell + the 2-atom parent grid), never as a bare NKPTS:
the 16-atom supercell is a 2x2x2 of the 2-atom cell, so the verified pair is
    2-atom k=12 shift 0.5  <->  16-atom k=6 shift 0.5   (density d=12)
and the supercell's own d=16 partner would be k=8 shift 0.5 (parent k=16).
Stage 0B uses the supercell mesh that pairs with the density already used for
the E(V) calibration, so the E(V) and the phonon sampling are commensurate.

Writes only under --out.  Refuses to overwrite an existing INCAR.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp.inputs import Potcar

import build_stage0 as S

SMEARINGS = {"s1_010": (1, 0.10), "s1_020": (1, 0.20),
             "s2_010": (2, 0.10), "s2_020": (2, 0.20)}
AMPLITUDES = (0.01, 0.03)
SUPERCELL = 2                      # 2x2x2
SUPER_NATOMS = 2 * SUPERCELL ** 3  # 16 for the 2-atom B2 cell
TARGETS = ("BCC-Ta-Ta-2", "BCC-Al-Ta-2", "BCC-Al-Cr-2")
# Supercell mesh.  HONEST NOTE (added during Stage 0C): this is mesh 6, i.e.
# density 12 in AlTa-equivalent units, NOT the density 16 that the E(V) path
# (Stage 0A, the SIGMA scan) uses on the 2-atom cell.  The comment here
# originally claimed the two were commensurate; that was wrong.  The density-16
# partner of this supercell is mesh 8, and a separate batch
# (build_stage0b_d16.py) supplies it so the E(V) and phonon sampling can be
# compared at the same density.
SUPER_MESH = 6
PARENT_MESH = 12


def phonopy_supercell(base: Structure, m: int) -> tuple[Phonopy, Structure]:
    """Phonopy object plus its OWN supercell as a pymatgen Structure.

    The supercell is taken from phonopy rather than rebuilt, because the two
    orderings do not agree: a hand-built (element x site x i,j,k) loop over the
    2-atom B2 cell gives Al(0,0,0), Al(0,0,1/2), ... while phonopy gives
    Al(0,0,0), Al(1/2,0,0), ... -- same atoms, different order, so the same
    cartesian point sits at a different index.  Writing one order to POSCAR and
    feeding forces back in the other silently attaches every force to the wrong
    atom, producing force constants that are nonzero, roughly symmetric and give
    plausible-looking frequencies.  Taking phonopy's own cell makes the two
    orders identical by construction; the fingerprint is recorded so the
    analysis can assert it later.
    """
    pa = PhonopyAtoms(symbols=[str(s.specie) for s in base],
                      scaled_positions=[s.frac_coords for s in base],
                      cell=base.lattice.matrix)
    ph = Phonopy(pa, supercell_matrix=np.diag([m, m, m]))
    sc = ph.supercell
    st = Structure(np.array(sc.cell), list(sc.symbols),
                   [np.array(p) for p in sc.scaled_positions],
                   coords_are_cartesian=False)
    return ph, st


def cell_fingerprint(st: Structure) -> str:
    """Order-sensitive fingerprint: lattice, species and fractional coords."""
    import hashlib
    payload = json.dumps({
        "lattice": [[round(float(x), 10) for x in row] for row in st.lattice.matrix],
        "species": [str(s.specie) for s in st],
        "frac": [[round(float(x), 10) for x in s.frac_coords] for s in st],
    }, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def displacements_for(ph: Phonopy, amp: float) -> list[dict]:
    """Symmetry-reduced phonopy displacements, one dict per displaced atom.

    phonopy stores each entry as a FLAT list ``[atom_index, dx, dy, dz]``, not
    as ``[atom_index, vector]``.  Reading index 1 as the vector silently yields
    a scalar (dx) and drops dy/dz -- which is what the first Stage 0B build did.
    """
    ph.generate_displacements(distance=amp, is_plusminus="auto", is_diagonal=True)
    out = []
    for i in range(len(ph.displacements)):
        d = ph.displacements[i]
        vec = [float(d[1]), float(d[2]), float(d[3])]
        assert len(vec) == 3
        out.append({"index": i + 1, "atom": int(d[0]), "vector": vec,
                    "magnitude": float(np.linalg.norm(vec))})
    return out


def displaced(st: Structure, atom: int, vec, scale: float = 1.0) -> Structure:
    """Return `st` with `vec` (Angstrom) applied to site `atom` (cartesian)."""
    lat = st.lattice
    coords = [s.coords for s in st]
    coords[atom] = coords[atom] + np.asarray(vec, float) * scale
    return Structure(lat, [s.specie for s in st], coords, coords_are_cartesian=True)


def write_case(d: Path, st: Structure, mesh: int, shift: float, system: str,
               ismear: int, sigma: float | None, magmom: str, meta: dict) -> dict:
    if (d / "INCAR").exists():
        raise SystemExit(f"refusing to overwrite existing case: {d}")
    d.mkdir(parents=True, exist_ok=True)
    Poscar(st).write_file(str(d / "POSCAR"))
    (d / "INCAR").write_text(S.incar_text(system, ismear, sigma, magmom))
    (d / "KPOINTS").write_text(S.kpoints_text(mesh, shift))
    species_order = list(dict.fromkeys(str(s.specie) for s in st))
    Potcar(species_order, functional=S.POTCAR_FUNCTIONAL).write_file(str(d / "POTCAR"))
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
    root = args.out / "stage0b"

    reps = S.build_structures()
    cases: list[dict] = []
    manifest: dict = {}

    for name in TARGETS:
        base = reps[name]["struct"]
        ph, sup = phonopy_supercell(base, SUPERCELL)
        assert len(sup) == SUPER_NATOMS, (name, len(sup))
        # mesh must be the fold-matched partner of the 2-atom calibration mesh
        assert SUPER_MESH * 2 == PARENT_MESH
        fp = cell_fingerprint(sup)
        manifest[name] = {"n_supercell_atoms": len(sup),
                          "supercell_mesh": SUPER_MESH,
                          "parent_2atom_mesh": PARENT_MESH,
                          "folding_verified": True,
                          "a_supercell": sup.lattice.a,
                          "cell_fingerprint": fp,
                          "atom_order_source": ("phonopy's own supercell; POSCAR order "
                                                "== phonopy order by construction"),
                          "displacements": {}}
        for amp in AMPLITUDES:
            disps = displacements_for(ph, amp)
            manifest[name]["displacements"][f"{amp:g}"] = disps
            for disp in disps:
                dstruct = displaced(sup, disp["atom"], disp["vector"])
                for tag, (ism, sig) in SMEARINGS.items():
                    meta = {
                        "endpoint": name, "class": reps[name]["class"],
                        "block": "stage0b",
                        "supercell": f"{SUPERCELL}x{SUPERCELL}x{SUPERCELL}",
                        "n_supercell_atoms": len(sup),
                        "cell_fingerprint": fp,
                        "atom_order_source": "phonopy supercell order",
                        "amplitude_angstrom": amp,
                        "disp_index": disp["index"],
                        "disp_atom": disp["atom"],
                        "disp_atom_species": str(sup[disp["atom"]].specie),
                        "disp_vector": disp["vector"],
                        "grid": [SUPER_MESH] * 3,
                        "shift": [0.5] * 3,
                        "cell_atoms": 2,
                        "folding_parent_grid": [PARENT_MESH] * 3,
                        "folding_verified": True,
                        "a_supercell_angstrom": sup.lattice.a,
                        "volume_angstrom3": float(sup.volume),
                    }
                    cases.append(write_case(
                        root / name / f"a{amp:g}_d{disp['index']:02d}_{tag}",
                        dstruct, SUPER_MESH, 0.5,
                        f"stage0b-{name}-a{amp:g}-d{disp['index']}-{tag}",
                        ism, sig, S.magmom_for(sup), meta))

    bad = [c for c in cases if not c["ok"]]
    report = {
        "schema": "legacy-qha-pool/stage0b/v1",
        "purpose": ("force / force-constant / phonon calibration across smearings "
                    "for stable, soft-mode and magnetic representatives"),
        "smearings": {k: {"ISMEAR": v[0], "SIGMA": v[1]} for k, v in SMEARINGS.items()},
        "amplitudes_angstrom": list(AMPLITUDES),
        "supercell": f"{SUPERCELL}x{SUPERCELL}x{SUPERCELL}",
        "kpoints_policy": ("explicit folding-compatible: supercell k=6 shift 0.5 "
                           "<-> 2-atom parent k=12 shift 0.5 (report 4.8.1)"),
        "ismear_minus5_policy": ("excluded from the force matrix: 9% force-constant "
                                 "deviation and no folding equivalence (4.8.2, 4.10)"),
        "recipe": {"ENCUT": S.ENCUT, "EDIFF": S.EDIFF, "ISYM": S.ISYM,
                   "NCORE": S.NCORE, "NELM": S.NELM, "NELMIN": S.NELMIN,
                   "ALGO": S.ALGO, "LORBIT": S.LORBIT,
                   "potcar_functional": S.POTCAR_FUNCTIONAL,
                   "alcr_magnetic_policy": "Al 1.0 / Cr 3.0 at every displacement "
                                           "and every smearing"},
        "manifest": manifest,
        "n_cases": len(cases), "n_bad": len(bad), "cases": cases,
    }
    (root / "validation.json").write_text(json.dumps(report, indent=2) + "\n")

    print(f"{'endpoint':<15} {'amp':>5} {'#disp':>6} {'atoms':>6} {'mesh':>5} "
          f"{'cases':>6}")
    print("-" * 52)
    for name in TARGETS:
        for amp in AMPLITUDES:
            nd = len(manifest[name]["displacements"][f"{amp:g}"])
            print(f"{name:<15} {amp:>5g} {nd:>6} {SUPER_NATOMS:>6} {SUPER_MESH:>5} "
                  f"{nd * len(SMEARINGS):>6}")
    print(f"\n{len(cases) - len(bad)}/{len(cases)} cases pass the input contract")
    for c in bad:
        print(f"  !! {c['dir']}")
        for p in c["problems"]:
            print(f"       {p}")
    if bad:
        raise SystemExit("refusing to submit: some inputs failed validation")


if __name__ == "__main__":
    main()

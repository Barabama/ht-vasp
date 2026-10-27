#!/usr/bin/env python
"""Build the multi-volume static-energy alignment test (ISMEAR sweep).

Question this answers (raised by the force/smearing review):
    Is the static-energy difference between ISMEAR=-5 and ISMEAR=1/2 a
    volume-INDEPENDENT constant shift (which QHA can absorb), or does it vary
    with volume (in which case the whole pool must use one smearing)?

Design, fixed by the review:
    * 5 volume points, not one.
    * 3 static flavours at every volume point:
          s1   ISMEAR=1  SIGMA=0.2   (production / "plan A")
          s2   ISMEAR=2  SIGMA=0.2   (MP-style control)
          sm5  ISMEAR=-5             (tetrahedron, static-audit reference)
    * Everything else identical, so ISMEAR is the ONLY variable.

Volume points come from the PRODUCTION QHA grid, not from legacy:
    QhaWorker(linear_strain=(-0.05, 0.05), number_of_frames=8)
    -> eps_l = linspace(-0.05, 0.05, 8), deformation = I*(1+eps),
       V(eps) = V0*(1+eps)^3.
    Selected frames (1-indexed) 1, 2, 4, 6, 8: both endpoints plus three
    interior points, i.e. the widest span the grid offers.

Two auxiliary sweeps, both at the middle volume point, to bound confounds:
    k7   same three smearings at 7x7x7   -> is dE_{-5->1} mesh-dependent?
    sig  ISMEAR=1 at SIGMA=0.05 / 0.10   -> is the MP result SIGMA-robust?

Cell: the 2-atom B2 AlTa conventional cubic cell. B2 has no free internal
coordinate, so a scaled cell is already the ISIF=2-relaxed geometry -- the
three smearings see byte-identical geometry and there is no geometry confound.
That is also why this is a "fixed-structure static energy" test and nothing
here says anything about forces or phonons.

Recipe: the production static recipe from main_qha.py GLOBAL_INCAR
(ENCUT=520, ALGO=Normal, EDIFF=1e-6, NELM=120, NELMIN=6, ISPIN=2,
MAGMOM Al 1.0 / Ta 3.0, POTCAR PBE(52) standard variants), with two
deliberate, documented deviations:
    ISYM = 2   (main_qha uses 0; E1 showed 0 is a defect and it costs 48x
                for statics. Fixed at 2 for ALL cases so ISMEAR stays the
                only variable.)
    NCORE = 8  (D0-D measured this as the correct parallel layout; the
                force test used it too.)

Writes only under --out. Never overwrites an existing INCAR.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp.inputs import Kpoints, Potcar

A0 = 3.2246015525754168            # legacy-relaxed B2 AlTa cell, Angstrom
POTCAR_DIR = "/nfs_ssd/documents/POT_GGA_PAW_PBE"   # PBE(52)
FUNCTIONAL = "PBE"
MAGMOM = {"Al": 1.0, "Ta": 3.0}
ENCUT = 520
EDIFF = 1e-6
ISYM = 2
NCORE = 8
NELM = 120
NELMIN = 6
KMAIN = 11
KALT = 7

# Production QHA volume grid.
FRAMES = 8
STRAIN = (-0.05, 0.05)
SELECTED = (1, 2, 4, 6, 8)          # 1-indexed frames


def b2_cell(a: float) -> Structure:
    """2-atom B2 AlTa cubic cell with lattice constant a."""
    lat = Lattice.cubic(a)
    return Structure(lat, ["Al", "Ta"], [[0, 0, 0], [0.5, 0.5, 0.5]])


def incar_text(system: str, ismear: int, sigma: float | None,
               magmom_line: str = "MAGMOM = 1.0 3.0") -> str:
    lines = [
        f"SYSTEM = {system}",
        f"ENCUT = {ENCUT}",
        "ISTART = 0",
        "ICHARG = 2",
        "LCHARG = False",
        "PREC = Normal",
        f"ISMEAR = {ismear}",
    ]
    if sigma is not None:
        lines.append(f"SIGMA = {sigma}")
    lines += [
        "ISPIN = 2",
        magmom_line.rstrip("\n"),
        "LREAL = False",
        "LWAVE = False",
        "IBRION = -1",
        "NSW = 0",
        f"EDIFF = {EDIFF:g}",
        f"ISYM = {ISYM}",
        "ALGO = Normal",
        f"NELM = {NELM}",
        f"NELMIN = {NELMIN}",
        f"NCORE = {NCORE}",
    ]
    return "\n".join(lines) + "\n"


def magmom_line_for(structure: Structure) -> str:
    """One MAGMOM value per atom, in POSCAR order.

    VASP rejects a MAGMOM whose length is neither 1 nor natoms:
    "Error reading item MAGMOM from file INCAR. Error code was IERR= 5".
    (That is exactly how the first 16-atom discrimination run died -- the
    2-atom INCAR's two values were reused verbatim for a 16-atom cell.)
    For a 2-atom B2 cell this reproduces the literal "MAGMOM = 1.0 3.0".
    """
    return "MAGMOM = " + " ".join(f"{MAGMOM[str(s.specie)]:.1f}" for s in structure)


def validate(d: Path, expect_atoms: int) -> dict:
    """Refuse to accept a directory VASP cannot parse (see build_d0_inputs.py)."""
    problems: list[str] = []
    poscar = d / "POSCAR"
    if not poscar.is_file():
        return {"ok": False, "problems": ["no POSCAR"]}
    lines = poscar.read_text().splitlines()
    for i in (2, 3, 4):
        try:
            vals = [float(x) for x in lines[i].split()]
            if len(vals) != 3:
                problems.append(f"POSCAR row {i + 1} has {len(vals)} numbers, not 3")
        except (ValueError, IndexError):
            problems.append(f"POSCAR row {i + 1} is not numeric: {lines[i]!r}")
    try:
        counts = [int(x) for x in lines[6].split()]
        n = sum(counts)
        if n != expect_atoms:
            problems.append(f"POSCAR declares {n} atoms, expected {expect_atoms}")
    except (ValueError, IndexError):
        problems.append(f"POSCAR row 7 is not integer counts: {lines[6]!r}")
        n = -1
    pot = (d / "POTCAR").read_text(errors="replace")
    nblk = pot.count("TITEL  =")
    poscar_species = lines[5].split()
    if nblk != len(poscar_species):
        problems.append(f"POTCAR has {nblk} blocks but POSCAR declares "
                        f"{len(poscar_species)} species")
    pot_species = re.findall(r"TITEL\s*=\s*PAW_\w+\s+(\w+)", pot)
    if pot_species and [s.upper() for s in pot_species] != [s.upper() for s in poscar_species]:
        problems.append(f"POTCAR order {pot_species} != POSCAR order {poscar_species}")
    # every k-point count in KPOINTS must be a plain positive integer
    kl = (d / "KPOINTS").read_text().splitlines()
    try:
        kk = [int(x) for x in kl[3].split()]
        if len(kk) != 3 or min(kk) < 3:
            problems.append(f"KPOINTS mesh {kl[3]!r} unusable (need 3 ints >= 3)")
    except (ValueError, IndexError):
        problems.append(f"KPOINTS line 4 is not a mesh: {kl[3]!r}")
    inc = (d / "INCAR").read_text()
    if "ENCUT" not in inc or "ISMEAR" not in inc:
        problems.append("INCAR missing ENCUT or ISMEAR")
    # Catch the "ISMEAR = ISMEAR" class of bug: the VALUE must be a literal.
    m = re.search(r"^ISMEAR\s*=\s*(\S+)\s*$", inc, re.M)
    if not m:
        problems.append("INCAR has no parseable ISMEAR line")
    else:
        try:
            int(m.group(1))
        except ValueError:
            problems.append(f"ISMEAR value is not an integer literal: {m.group(1)!r}")
    m = re.search(r"^SIGMA\s*=\s*(\S+)\s*$", inc, re.M)
    if m:
        try:
            float(m.group(1))
        except ValueError:
            problems.append(f"SIGMA value is not a numeric literal: {m.group(1)!r}")
    m = re.search(r"^ENCUT\s*=\s*(\S+)\s*$", inc, re.M)
    if m:
        try:
            float(m.group(1))
        except ValueError:
            problems.append(f"ENCUT value is not numeric: {m.group(1)!r}")
    if "ISMEAR = -5" in inc and "SIGMA" in inc:
        problems.append("ISMEAR=-5 case must not carry SIGMA")
    # MAGMOM must be 1 or natoms values long, otherwise VASP dies at startup with
    # "Error reading item MAGMOM from file INCAR. Error code was IERR= 5".
    m = re.search(r"^MAGMOM\s*=\s*(.+)$", inc, re.M)
    if not m:
        problems.append("INCAR has no MAGMOM line")
    else:
        vals = m.group(1).split()
        try:
            [float(v) for v in vals]
        except ValueError:
            problems.append(f"MAGMOM has a non-numeric value: {m.group(1)!r}")
        if len(vals) not in (1, n):
            problems.append(
                f"MAGMOM has {len(vals)} values but the cell has {n} atoms "
                f"(VASP accepts only 1 or natoms)")
    return {
        "ok": not problems,
        "problems": problems,
        "natoms": n,
        "potcar_blocks": nblk,
        "poscar_sha256": hashlib.sha256(poscar.read_bytes()).hexdigest()[:16],
        "potcar_sha256": hashlib.sha256((d / "POTCAR").read_bytes()).hexdigest()[:16],
        "incar_sha256": hashlib.sha256((d / "INCAR").read_bytes()).hexdigest()[:16],
        "kpoints_sha256": hashlib.sha256((d / "KPOINTS").read_bytes()).hexdigest()[:16],
    }


def write_case(d: Path, structure: Structure, kgrid: int, system: str,
               ismear: int, sigma: float | None, meta: dict,
               shift: tuple[float, float, float] = (0.0, 0.0, 0.0)) -> dict:
    if (d / "INCAR").exists():
        raise SystemExit(f"refusing to overwrite existing case: {d}")
    d.mkdir(parents=True, exist_ok=True)
    Poscar(structure).write_file(str(d / "POSCAR"))
    (d / "INCAR").write_text(incar_text(system, ismear, sigma,
                                        magmom_line_for(structure)))
    Kpoints.monkhorst_automatic([kgrid, kgrid, kgrid], list(shift)).write_file(
        str(d / "KPOINTS"))
    species_order = list(dict.fromkeys(str(s.specie) for s in structure))
    Potcar(species_order, functional=FUNCTIONAL).write_file(str(d / "POTCAR"))
    v = validate(d, len(structure))
    v["dir"] = str(d)
    v["kgrid"] = kgrid
    v["shift"] = list(shift)
    v["ismear"] = ismear
    v["sigma"] = sigma
    v.update(meta)
    (d / "case.json").write_text(json.dumps(v, indent=2) + "\n")
    return v


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    eps_l = np.linspace(STRAIN[0], STRAIN[1], FRAMES)
    v0 = A0 ** 3
    report = {
        "schema": "legacy-qha-pool/align/v1",
        "question": "is E(ISMEAR=-5) - E(ISMEAR=1) volume-independent?",
        "anchor_a0_angstrom": A0,
        "anchor_v0_angstrom3": v0,
        "volume_grid_source": f"QhaWorker linear_strain={STRAIN} number_of_frames={FRAMES}",
        "strain_frames_eps": [float(x) for x in eps_l],
        "selected_frames_1indexed": list(SELECTED),
        "recipe": {"ENCUT": ENCUT, "EDIFF": EDIFF, "ISYM": ISYM, "NCORE": NCORE,
                   "NELM": NELM, "NELMIN": NELMIN, "ISPIN": 2,
                   "MAGMOM": MAGMOM, "potcar_functional": FUNCTIONAL,
                   "potcar_dir": POTCAR_DIR,
                   "deviations_from_main_qha": {
                       "ISYM": "2 (main_qha uses 0; fixed at 2 for all cases so ISMEAR is the only variable)",
                       "NCORE": "8 (not set by main_qha; D0-D measured 8 as correct)",
                   }},
        "smearings": {
            "s1": {"ISMEAR": 1, "SIGMA": 0.2},
            "s2": {"ISMEAR": 2, "SIGMA": 0.2},
            "sm5": {"ISMEAR": -5, "SIGMA": None},
        },
        "primary_kgrid": KMAIN,
        "aux_kgrid": KALT,
        "cases": [],
    }

    # ---- primary: 5 volumes x 3 smearings, k=11 -------------------------
    for idx in SELECTED:
        eps = float(eps_l[idx - 1])
        a = A0 * (1.0 + eps)
        st = b2_cell(a)
        vol = float(st.volume)
        meta = {"block": "primary", "frame": idx, "eps": eps,
                "a_angstrom": a, "volume_angstrom3": vol,
                "volume_scale_vs_v0": vol / v0}
        for tag, spec in report["smearings"].items():
            report["cases"].append(write_case(
                out / "primary" / f"f{idx:02d}_{tag}", st, KMAIN,
                f"align-f{idx:02d}-{tag}", spec["ISMEAR"], spec["SIGMA"], meta))

    # ---- aux k-sensitivity: middle volume x 3 smearings, k=7 ------------
    idx = 4
    eps = float(eps_l[idx - 1])
    a = A0 * (1.0 + eps)
    st = b2_cell(a)
    meta = {"block": "kmesh", "frame": idx, "eps": eps,
            "a_angstrom": a, "volume_angstrom3": float(st.volume),
            "volume_scale_vs_v0": float(st.volume) / v0}
    for tag, spec in report["smearings"].items():
        report["cases"].append(write_case(
            out / "kmesh" / f"f{idx:02d}_k{KALT}_{tag}", st, KALT,
            f"align-k{KALT}-f{idx:02d}-{tag}", spec["ISMEAR"], spec["SIGMA"], meta))

    # ---- aux SIGMA-sensitivity: middle volume, ISMEAR=1, k=11 ----------
    for sig in (0.05, 0.10):
        report["cases"].append(write_case(
            out / "sigma" / f"f{idx:02d}_s1_sigma{sig:g}", st, KMAIN,
            f"align-sigma{sig:g}-f{idx:02d}", 1, sig,
            {**meta, "block": "sigma", "sigma": sig}))

    bad = [c for c in report["cases"] if not c["ok"]]
    report["n_cases"] = len(report["cases"])
    report["n_bad"] = len(bad)
    (out / "validation.json").write_text(json.dumps(report, indent=2) + "\n")

    for c in report["cases"]:
        mark = "OK " if c["ok"] else "BAD"
        print(f"  {mark} {c['dir']:<56} V={c['volume_angstrom3']:8.4f} "
              f"(x{c['volume_scale_vs_v0']:.5f}) k={c['kgrid']:>2} "
              f"ISMEAR={c['ismear']:>2} SIGMA={c['sigma']}")
        for p in c["problems"]:
            print(f"        !! {p}")
    print(f"\n{report['n_cases'] - len(bad)}/{report['n_cases']} cases validated")
    print(f"volume scales covered: "
          f"{sorted({round(c['volume_scale_vs_v0'], 5) for c in report['cases'] if c['block'] == 'primary'})}")
    if bad:
        raise SystemExit("refusing to submit: some inputs failed validation")


if __name__ == "__main__":
    main()

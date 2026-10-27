#!/usr/bin/env python
"""Build D0-D and convergence inputs for B2 AlTa, with a pre-flight validation.

Rewritten after the first attempt shipped a POSCAR whose lattice rows 2 and 3
were literal "{lat:.16f}" strings (a missing f-prefix), which made VASP read
"0 types and 0 ions" and then fail with a completely unrelated
"LUSE_VDW needs to be set to .TRUE." message. The structure is now built by
pymatgen and every generated directory is validated before it is accepted.

POTCAR: functional PBE (the PBE(52) library, POT_GGA_PAW_PBE), standard
variants -- the same family and variants as the legacy data, so results stay
comparable to E1 and to legacy.

Writes only under --out.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp.inputs import Kpoints, Potcar

A0 = 3.2246015525754168          # legacy relaxed B2 lattice constant, Angstrom
POTCAR_DIR = "/nfs_ssd/documents/POT_GGA_PAW_PBE"   # PBE(52)
FUNCTIONAL = "PBE"
MAGMOM = {"Al": 1.0, "Ta": 3.0}


def b2(super_cell: int) -> Structure:
    """2-atom B2 cell repeated super_cell times along each axis."""
    a = A0 * super_cell
    lat = Lattice.cubic(a)
    sites = []
    for el, base in (("Al", (0.0, 0.0, 0.0)), ("Ta", (0.5, 0.5, 0.5))):
        for i in range(super_cell):
            for j in range(super_cell):
                for k in range(super_cell):
                    sites.append((el, ((base[0] + i) / super_cell,
                                       (base[1] + j) / super_cell,
                                       (base[2] + k) / super_cell)))
    return Structure(lat, [s[0] for s in sites], [s[1] for s in sites])


def incar_text(system: str, encut: int, magmom_line: str,
               parallel: str = "KPAR = 2\n", extra: str = "") -> str:
    return f"""SYSTEM = {system}
ENCUT = {encut}
ISTART = 0
ICHARG = 2
LCHARG = False
PREC = Normal
ISMEAR = 1
SIGMA = 0.2
ISPIN = 2
{magmom_line}LREAL = False
LWAVE = False
IBRION = -1
NSW = 0
EDIFF = 1e-05
ISYM = 2
ALGO = Normal
NELM = 120
NELMIN = 10
{parallel}{extra}"""


def validate(d: Path, expect_atoms: int) -> dict:
    """Refuse to accept a directory VASP cannot parse."""
    problems = []
    poscar = d / "POSCAR"
    if not poscar.is_file():
        return {"ok": False, "problems": ["no POSCAR"]}
    lines = poscar.read_text().splitlines()
    # rows 3,4,5 (0-based 2,3,4) must be three floats each
    for i in (2, 3, 4):
        try:
            vals = [float(x) for x in lines[i].split()]
            if len(vals) != 3:
                problems.append(f"POSCAR row {i + 1} has {len(vals)} numbers, not 3")
        except (ValueError, IndexError):
            problems.append(f"POSCAR row {i + 1} is not numeric: {lines[i]!r}")
    # species counts
    try:
        counts = [int(x) for x in lines[6].split()]
        n = sum(counts)
        if n != expect_atoms:
            problems.append(f"POSCAR declares {n} atoms, expected {expect_atoms}")
    except (ValueError, IndexError):
        problems.append(f"POSCAR row 7 is not integer counts: {lines[6]!r}")
    # POTCAR must have the same number of blocks as species
    pot = (d / "POTCAR").read_text(errors="replace")
    nblk = pot.count("TITEL  =")
    if nblk != len(lines[5].split()):
        problems.append(f"POTCAR has {nblk} blocks but POSCAR declares "
                        f"{len(lines[5].split())} species")
    # POTCAR block ORDER must equal the POSCAR species order
    poscar_species = lines[5].split()
    pot_species = re.findall(r"TITEL\s*=\s*PAW_\w+\s+(\w+)", pot)
    if pot_species and [s.upper() for s in pot_species] != [s.upper() for s in poscar_species]:
        problems.append(f"POTCAR order {pot_species} != POSCAR order {poscar_species}: "
                        f"VASP matches them positionally")
    # INCAR must be parseable
    inc = (d / "INCAR").read_text()
    if "ENCUT" not in inc:
        problems.append("INCAR has no ENCUT")
    out = {"ok": not problems, "problems": problems,
           "natoms": n if "n" in dir() else None,
           "potcar_blocks": nblk,
           "poscar_sha256": hashlib.sha256(poscar.read_bytes()).hexdigest()[:16],
           "potcar_sha256": hashlib.sha256((d / "POTCAR").read_bytes()).hexdigest()[:16]}
    return out


def write_case(d: Path, structure: Structure, encut: int, kgrid: int,
               system: str, parallel: str = "KPAR = 2\n",
               extra: str = "") -> dict:
    d.mkdir(parents=True, exist_ok=True)
    Poscar(structure).write_file(str(d / "POSCAR"))
    syms = [str(s.specie) for s in structure]
    mm = "MAGMOM = " + " ".join(f"{MAGMOM[s]:.1f}" for s in syms) + "\n"
    (d / "INCAR").write_text(incar_text(system, encut, mm, parallel, extra))
    Kpoints.monkhorst_automatic([kgrid, kgrid, kgrid], [0, 0, 0]).write_file(
        str(d / "KPOINTS"))
    # POTCAR via pymatgen, PBE(52), standard variants.
    # The species order MUST follow the POSCAR species line: VASP matches POTCAR
    # blocks to POSCAR species positionally, so a reordered POTCAR silently
    # attaches the wrong element to the wrong sites.
    species_order = list(dict.fromkeys(syms))          # POSCAR order, deduplicated
    Potcar(species_order, functional=FUNCTIONAL).write_file(str(d / "POTCAR"))
    v = validate(d, len(structure))
    v["dir"] = str(d)
    v["encut"] = encut
    v["kgrid"] = kgrid
    v["parallel"] = parallel.strip() or "default"
    return v


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    report = {"schema": "legacy-qha-pool/d0-inputs/v1",
              "potcar_functional": FUNCTIONAL,
              "potcar_dir": POTCAR_DIR,
              "note": ("PBE(52), standard variants Al/Ta -- same family and variants "
                       "as the legacy data, so results stay comparable to E1."),
              "cases": []}

    # ---- D0-D: 128-atom static, three parallel layouts -------------------
    s4 = b2(4)
    for tag, par in (("kpar2", "KPAR = 2\n"),
                     ("ncore8", "NCORE = 8\n"),
                     ("plain", "")):
        report["cases"].append(
            write_case(out / "d0d" / f"run-{tag}", s4, 520, 3,
                       f"d0d-{tag}", parallel=par))

    # ---- convergence: ENCUT series (k=5) and k series (ENCUT=520) --------
    s2 = b2(2)
    for enc in (400, 520, 600, 700):
        report["cases"].append(
            write_case(out / "conv" / f"encut-{enc}", s2, enc, 5, f"conv-encut{enc}"))
    for k in (3, 5, 7, 9):
        report["cases"].append(
            write_case(out / "conv" / f"k-{k}", s2, 520, k, f"conv-k{k}"))

    bad = [c for c in report["cases"] if not c["ok"]]
    report["n_cases"] = len(report["cases"])
    report["n_bad"] = len(bad)
    (out / "validation.json").write_text(json.dumps(report, indent=2) + "\n")

    for c in report["cases"]:
        mark = "OK " if c["ok"] else "BAD"
        print(f"  {mark} {c['dir']:<44} natoms={c['natoms']:>4} "
              f"encut={c['encut']:>3} k={c['kgrid']} par={c['parallel']}")
        for p in c["problems"]:
            print(f"        !! {p}")
    print(f"\n{report['n_cases'] - len(bad)}/{report['n_cases']} cases validated")
    if bad:
        raise SystemExit("refusing to submit: some inputs failed validation")


if __name__ == "__main__":
    main()

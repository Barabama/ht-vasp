#!/usr/bin/env python
"""Regression tests for input_contract.py.

Every case here is a failure this project actually shipped, or a near-miss the
validator must keep catching.  The point is that a future edit to
``input_contract.py`` cannot silently drop a check: run this and the known-bad
inputs must still be rejected, the known-good ones must still pass.

The three real incidents (report §9.3):

  1. POSCAR lattice rows written as literal "{lat:.16f}" -> VASP read
     "POSCAR found : 0 types and 0 ions", then died with an unrelated
     LUSE_VDW message.
  2. INCAR got "ISMEAR = ISMEAR" because a dict was unpacked as a 2-tuple;
     20 directories passed a check that only looked for the substring.
  3. A 16-atom cell got the 2-atom MAGMOM -> VASP refused with
     "Error reading item MAGMOM ... IERR= 5", killing 14 cases.

Run:
    python test_input_contract.py            # static layer only, no VASP
    python test_input_contract.py --dirs D1 D2   # also check real trees

Exits non-zero on any failure.
"""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from input_contract import validate_dir  # noqa: E402

# ---- a minimal but VALID 2-atom B2 cell, used as the base for mutations ----
POSCAR_OK = """B2 AlTa
1.0
   3.2246015525754168    0.0000000000000000    0.0000000000000000
   0.0000000000000000    3.2246015525754168    0.0000000000000000
   0.0000000000000000    0.0000000000000000    3.2246015525754168
Al Ta
1 1
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000
   0.5000000000000000    0.5000000000000000    0.5000000000000000
"""

INCAR_OK = """SYSTEM = test
ENCUT = 520
ISTART = 0
ICHARG = 2
LCHARG = False
PREC = Normal
ISMEAR = 1
SIGMA = 0.2
ISPIN = 2
MAGMOM = 1.0 3.0
LREAL = False
LWAVE = False
IBRION = -1
NSW = 0
EDIFF = 1e-06
ISYM = 2
ALGO = Normal
NELM = 120
NELMIN = 6
NCORE = 8
"""

KPOINTS_OK = """Automatic kpoint scheme
0
Monkhorst-Pack
16 16 16
0.5 0.5 0.5
"""

# POTCAR blocks are assembled from a real file at run time; this is the order
# the POSCAR above requires.
POTCAR_AL = "   TITEL  = PAW_PBE Al 04Jan2001\n   ZVAL   =    3.00000000\n"
POTCAR_TA = "   TITEL  = PAW_PBE Ta 17Jan2003\n   ZVAL   =    5.00000000\n"
POTCAR_OK = POTCAR_AL + POTCAR_TA


def make(base: Path, **over) -> Path:
    d = base / "case"
    d.mkdir(parents=True, exist_ok=True)
    (d / "POSCAR").write_text(over.get("POSCAR", POSCAR_OK))
    (d / "INCAR").write_text(over.get("INCAR", INCAR_OK))
    (d / "KPOINTS").write_text(over.get("KPOINTS", KPOINTS_OK))
    (d / "POTCAR").write_text(over.get("POTCAR", POTCAR_OK))
    for name in over.get("remove", ()):
        f = d / name
        if f.exists():
            f.unlink()
    return d


# (name, overrides, must_fail, expected substring in a problem message)
CASES = [
    ("valid 2-atom B2 cell", {}, False, None),

    # ---- incident 1: POSCAR lattice rows are literal format strings --------
    ("POSCAR rows are literal '{lat:.16f}'", {
        "POSCAR": POSCAR_OK.replace(
            "   0.0000000000000000    3.2246015525754168    0.0000000000000000",
            "   0.0 {lat:.16f} 0.0")}, True, "not numeric"),

    # ---- incident 2: INCAR values are the tag names ------------------------
    ("INCAR ISMEAR = ISMEAR", {
        "INCAR": INCAR_OK.replace("ISMEAR = 1", "ISMEAR = ISMEAR")},
        True, "not an integer literal"),
    ("INCAR SIGMA = SIGMA", {
        "INCAR": INCAR_OK.replace("SIGMA = 0.2", "SIGMA = SIGMA")},
        True, "not numeric"),

    # ---- incident 3: MAGMOM length is neither 1 nor NIONS ------------------
    ("MAGMOM has 2 values for a 16-atom cell", {
        "POSCAR": POSCAR_OK.replace("Al Ta\n1 1\n", "Al Ta\n8 8\n").replace(
            "   0.0000000000000000    0.0000000000000000    0.0000000000000000\n"
            "   0.5000000000000000    0.5000000000000000    0.5000000000000000\n",
            "".join(f"   {i / 16:.6f}   0.0 0.0\n" for i in range(16)))},
        True, "MAGMOM"),

    # ---- POTCAR / POSCAR mismatches ---------------------------------------
    ("POTCAR block order reversed", {
        "POTCAR": POTCAR_TA + POTCAR_AL}, True, "element order"),
    ("POTCAR has one block, POSCAR two species", {
        "POTCAR": POTCAR_AL}, True, "blocks"),
    ("POSCAR declares 3 atoms, lists 2 coordinates", {
        "POSCAR": POSCAR_OK.replace("Al Ta\n1 1\n", "Al Ta\n1 2\n")},
        True, "coordinate lines"),

    # ---- KPOINTS ----------------------------------------------------------
    ("KPOINTS mesh is not integers", {
        "KPOINTS": KPOINTS_OK.replace("16 16 16", "16 16 x")},
        True, "not integers"),
    ("KPOINTS mesh has 2 entries", {
        "KPOINTS": KPOINTS_OK.replace("16 16 16", "16 16")},
        True, "not 3"),

    # ---- INCAR structural problems ---------------------------------------
    ("NCORE and KPAR both set", {
        "INCAR": INCAR_OK.replace("NCORE = 8", "NCORE = 8\nKPAR = 2")},
        True, "NCORE and KPAR"),
    ("IBRION is not an integer", {
        "INCAR": INCAR_OK.replace("IBRION = -1", "IBRION = eight")},
        True, "not an integer literal"),
    ("LCHARG is not a logical", {
        "INCAR": INCAR_OK.replace("LCHARG = False", "LCHARG = maybe")},
        True, "not a logical value"),
    ("ICHARG=1 without a CHGCAR", {
        "INCAR": INCAR_OK.replace("ICHARG = 2", "ICHARG = 1\nISTART = 0")},
        True, "needs a CHGCAR"),

    # ---- missing files ----------------------------------------------------
    ("no KPOINTS", {"remove": ("KPOINTS",)}, True, "no KPOINTS"),
    ("no POTCAR", {"remove": ("POTCAR",)}, True, "no POTCAR"),

    # ---- cases that must still PASS --------------------------------------
    ("ISMEAR=-5 with no SIGMA is fine", {
        "INCAR": INCAR_OK.replace("ISMEAR = 1\nSIGMA = 0.2\n", "ISMEAR = -5\n")},
        False, None),
    ("MAGMOM with a single value is fine", {
        "INCAR": INCAR_OK.replace("MAGMOM = 1.0 3.0", "MAGMOM = 2.0")},
        False, None),
    ("even mesh + shift 0.5 passes (contains Gamma)", {}, False, None),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="*", type=Path, default=[],
                    help="also validate these real trees; every case must PASS")
    args = ap.parse_args()

    npass = nfail = 0
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        for name, over, must_fail, expect in CASES:
            d = make(base / name.replace("/", "_"), **over)
            res = validate_dir(d)
            ok = (not res["ok"]) if must_fail else res["ok"]
            detail = ""
            if must_fail and res["ok"]:
                detail = " -> accepted but should have been rejected"
            elif must_fail and expect and not any(expect in p for p in res["problems"]):
                detail = (f" -> rejected, but for the wrong reason; expected "
                          f"{expect!r}, got {res['problems']}")
                ok = False
            elif not must_fail and not res["ok"]:
                detail = f" -> rejected: {res['problems']}"
            print(f"  {'PASS' if ok else 'FAIL'}  {name}{detail}")
            npass += ok
            nfail += not ok

    for tree in args.dirs:
        if not tree.is_dir():
            print(f"  SKIP  {tree} (not a directory)")
            continue
        n = bad = skipped = 0
        for incar in sorted(tree.rglob("INCAR")):
            parts = incar.parent.relative_to(tree).parts
            if any(p.startswith(".") for p in parts):
                continue
            # Archived trees are deliberately kept and are EXPECTED to fail:
            # this project never deletes a superseded input, it renames it with
            # a marker (report §9.3).  Counting them as failures would make this
            # check useless the moment anything is archived.
            if any(ARCHIVE_MARKERS in p for p in parts for ARCHIVE_MARKERS in
                   (".badmagmom", ".broken", ".oldmesh", ".partial", ".superseded")):
                skipped += 1
                continue
            n += 1
            res = validate_dir(incar.parent)
            if not res["ok"]:
                bad += 1
                print(f"  FAIL  {incar.parent} -> {res['problems']}")
        ok = bad == 0
        note = f" ({skipped} archived case(s) skipped)" if skipped else ""
        print(f"  {'PASS' if ok else 'FAIL'}  real tree {tree}: {n - bad}/{n} pass{note}")
        npass += ok
        nfail += not ok

    print(f"\n{npass} passed, {nfail} failed")
    if nfail:
        sys.exit(1)


if __name__ == "__main__":
    main()

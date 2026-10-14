"""Tests for the VASP input gate.

Every case here is a failure this project actually shipped, or a near-miss the
gate must keep catching.  The point is that a future edit cannot silently drop
a check: the known-bad inputs must still be rejected, the known-good ones must
still pass.

The three real incidents:

  1. POSCAR lattice rows written as literal ``"{lat:.16f}"`` -> VASP read
     "POSCAR found : 0 types and 0 ions", then died with an unrelated LUSE_VDW
     message.
  2. INCAR got ``ISMEAR = ISMEAR`` because a dict was unpacked as a 2-tuple, and
     20 directories passed a check that only looked for the substring.
  3. A 16-atom cell got the 2-atom MAGMOM -> VASP refused with "Error reading
     item MAGMOM ... IERR= 5", killing 14 cases.
"""

from __future__ import annotations

from pathlib import Path

import pydantic
import pytest

from htvasp.dfpt_qha.gate import iter_cases, read_kpoints, validate_dir

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

POTCAR_AL = "   TITEL  = PAW_PBE Al 04Jan2001\n   ZVAL   =    3.00000000\n"
POTCAR_TA = "   TITEL  = PAW_PBE Ta 17Jan2003\n   ZVAL   =    5.00000000\n"
POTCAR_OK = POTCAR_AL + POTCAR_TA

SIXTEEN_ATOM_POSCAR = (
    POSCAR_OK.replace("Al Ta\n1 1\n", "Al Ta\n8 8\n").replace(
        "   0.0000000000000000    0.0000000000000000    0.0000000000000000\n"
        "   0.5000000000000000    0.5000000000000000    0.5000000000000000\n",
        "".join(f"   {i / 16:.6f}   0.0 0.0\n" for i in range(16)),
    )
)


def make(base: Path, **overrides) -> Path:
    directory = base / "case"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "POSCAR").write_text(overrides.get("POSCAR", POSCAR_OK))
    (directory / "INCAR").write_text(overrides.get("INCAR", INCAR_OK))
    (directory / "KPOINTS").write_text(overrides.get("KPOINTS", KPOINTS_OK))
    (directory / "POTCAR").write_text(overrides.get("POTCAR", POTCAR_OK))
    for name in overrides.get("remove", ()):
        victim = directory / name
        if victim.exists():
            victim.unlink()
    return directory


# (label, overrides, must_fail, expected substring in a problem message)
CASES = [
    ("valid 2-atom B2 cell", {}, False, None),
    # ---- incident 1: POSCAR lattice rows are literal format strings --------
    (
        "POSCAR rows are literal '{lat:.16f}'",
        {
            "POSCAR": POSCAR_OK.replace(
                "   0.0000000000000000    3.2246015525754168    0.0000000000000000",
                "   0.0 {lat:.16f} 0.0",
            )
        },
        True,
        "not numeric",
    ),
    # ---- incident 2: INCAR values are the tag names ------------------------
    (
        "INCAR ISMEAR = ISMEAR",
        {"INCAR": INCAR_OK.replace("ISMEAR = 1", "ISMEAR = ISMEAR")},
        True,
        "not an integer literal",
    ),
    (
        "INCAR SIGMA = SIGMA",
        {"INCAR": INCAR_OK.replace("SIGMA = 0.2", "SIGMA = SIGMA")},
        True,
        "not numeric",
    ),
    # ---- incident 3: MAGMOM length is neither 1 nor NIONS ------------------
    (
        "MAGMOM has 2 values for a 16-atom cell",
        {"POSCAR": SIXTEEN_ATOM_POSCAR},
        True,
        "MAGMOM",
    ),
    # ---- POTCAR / POSCAR mismatches ---------------------------------------
    ("POTCAR block order reversed", {"POTCAR": POTCAR_TA + POTCAR_AL}, True, "element order"),
    ("POTCAR has one block, POSCAR two species", {"POTCAR": POTCAR_AL}, True, "blocks"),
    (
        "POSCAR declares 3 atoms, lists 2 coordinates",
        {"POSCAR": POSCAR_OK.replace("Al Ta\n1 1\n", "Al Ta\n1 2\n")},
        True,
        "coordinate lines",
    ),
    # ---- KPOINTS ----------------------------------------------------------
    ("KPOINTS mesh is not integers", {"KPOINTS": KPOINTS_OK.replace("16 16 16", "16 16 x")}, True, "not integers"),
    ("KPOINTS mesh has 2 entries", {"KPOINTS": KPOINTS_OK.replace("16 16 16", "16 16")}, True, "not 3"),
    # ---- INCAR structural problems ---------------------------------------
    ("NCORE and KPAR both set", {"INCAR": INCAR_OK.replace("NCORE = 8", "NCORE = 8\nKPAR = 2")}, True, "NCORE and KPAR"),
    ("IBRION is not an integer", {"INCAR": INCAR_OK.replace("IBRION = -1", "IBRION = eight")}, True, "not an integer literal"),
    ("LCHARG is not a logical", {"INCAR": INCAR_OK.replace("LCHARG = False", "LCHARG = maybe")}, True, "not a logical value"),
    ("ICHARG=1 without a CHGCAR", {"INCAR": INCAR_OK.replace("ICHARG = 2", "ICHARG = 1\nISTART = 0")}, True, "needs a CHGCAR"),
    # ---- missing files ----------------------------------------------------
    ("no KPOINTS", {"remove": ("KPOINTS",)}, True, "no KPOINTS"),
    ("no POTCAR", {"remove": ("POTCAR",)}, True, "no POTCAR"),
    # ---- cases that must still PASS --------------------------------------
    ("ISMEAR=-5 with no SIGMA is fine", {"INCAR": INCAR_OK.replace("ISMEAR = 1\nSIGMA = 0.2\n", "ISMEAR = -5\n")}, False, None),
    ("MAGMOM with a single value is fine", {"INCAR": INCAR_OK.replace("MAGMOM = 1.0 3.0", "MAGMOM = 2.0")}, False, None),
    ("even mesh + shift 0.5 passes (contains Gamma)", {}, False, None),
]


@pytest.mark.parametrize(
    ("label", "overrides", "must_fail", "needle"),
    CASES,
    ids=[case[0] for case in CASES],
)
def test_gate_case(tmp_path, label, overrides, must_fail, needle):
    result = validate_dir(make(tmp_path / label.replace("/", "_"), **overrides))
    if not must_fail:
        assert result.ok, f"{label} -> rejected: {result.problems}"
        return
    assert not result.ok, f"{label} -> accepted but should have been rejected"
    assert any(needle in problem for problem in result.problems), (
        f"{label} -> rejected, but for the wrong reason; expected {needle!r}, "
        f"got {list(result.problems)}"
    )


# --------------------------------------------------------------------------
# the Gamma rule the gate shares with the recipe contract
# --------------------------------------------------------------------------


def test_odd_mesh_with_half_shift_is_warned(tmp_path):
    """A gap the previous gate had: odd n + shift 0.5 excludes Gamma silently.

    Its condition was ``gamma_centred is False AND shift == 0``, so the
    uniform-half-shift case fell through every branch and produced no warning
    at all.
    """

    directory = make(
        tmp_path, KPOINTS=KPOINTS_OK.replace("16 16 16", "7 7 7")
    )
    result = validate_dir(directory)
    assert result.ok, "an unusual mesh is a warning, not a rejection"
    assert any("excludes Gamma" in warning for warning in result.warnings), result.warnings
    assert any("axis 0" in warning for warning in result.warnings), result.warnings


def test_unmodelled_shift_reports_unknown_not_excluded(tmp_path):
    """A mixed shift is not modelled, so membership is None, not False."""

    info = read_kpoints(make(tmp_path, KPOINTS=KPOINTS_OK.replace("0.5 0.5 0.5", "0.0 0.5 0.5")))
    assert info.gamma_centered is None
    assert "not modelled" in (info.gamma_rule or "")


def test_gate_result_is_frozen():
    """A gate whose verdict can be edited after the fact is not a gate."""

    result = validate_dir(Path("/nonexistent"))
    with pytest.raises(pydantic.ValidationError):
        result.ok = True


def test_missing_directory_reports_every_absent_file():
    result = validate_dir(Path("/nonexistent"))
    assert not result.ok
    for name in ("POSCAR", "POTCAR", "KPOINTS", "INCAR"):
        assert any(name in problem for problem in result.problems), result.problems


# --------------------------------------------------------------------------
# iter_cases and the dry run (the dry run needs a compute node)
# --------------------------------------------------------------------------


def test_iter_cases_skips_hidden_directories(tmp_path):
    """`--dryrun` writes its scratch copies under a dotted dir in the same tree."""

    make(tmp_path / "visible")
    make(tmp_path / ".dryscratch" / "leftover")
    found = [p.relative_to(tmp_path).as_posix() for p in iter_cases(tmp_path)]
    assert found == ["visible/case"], found


@pytest.mark.skipif(
    not Path("/nfs_ssd/tmp/legacy-qha-pool/stage0").is_dir(),
    reason="stage0 tree not present",
)
def test_real_stage0_tree_passes_except_archived():
    """The tree the gate was actually used on.

    Archived trees are deliberately kept and are EXPECTED to fail: this project
    never deletes a superseded input, it renames it with a marker.  Counting
    them as failures would make this check useless the moment anything is
    archived.
    """

    markers = (".badmagmom", ".broken", ".oldmesh", ".partial", ".superseded")
    for incar in sorted(Path("/nfs_ssd/tmp/legacy-qha-pool/stage0").rglob("INCAR")):
        parts = incar.parent.parts
        if any(marker in part for part in parts for marker in markers):
            continue
        result = validate_dir(incar.parent)
        assert result.ok, f"{incar.parent}: {list(result.problems)}"

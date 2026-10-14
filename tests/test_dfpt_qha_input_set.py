"""Tests for the DFPT-QHA input generation.

The generator and the gate are tested together on purpose: a generator whose
output the gate rejects is useless, and the gate's first real workload is
these inputs. Building them together is what surfaced that the gate did not
understand VASP's ``n*value`` repeat syntax, which pymatgen emits whenever
moments repeat.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pymatgen.core import Lattice, Structure

from htvasp.dfpt_qha.gate import parse_repeated_floats, validate_dir
from htvasp.dfpt_qha.input_set import (
    DFPT_ELECTRONIC,
    InputSetError,
    STAGES,
    STATIC_ELECTRONIC,
    incar_for,
    kpoints_for,
    magmom_for,
    mesh_for,
    potcar_symbols,
    round_even,
    write_stage_inputs,
)
from htvasp.dfpt_qha.recipe import Recipe

POOL_RECIPES = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")

B2 = Structure(Lattice.cubic(3.2246), ["Al", "Ta"], [[0, 0, 0], [0.5, 0.5, 0.5]])
ALNB = Structure(Lattice.cubic(3.3176), ["Al", "Nb"], [[0, 0, 0], [0.5, 0.5, 0.5]])

BASE_RECIPE = {
    "structure_class": "BCC_B2",
    "site_occupancy": {"1a": "Al", "1b": "Ta"},
    "representation": "conventional",
    "sym_reduce": True,
    "supercell": [2, 2, 2],
    "kpoints": {
        "grid": [16, 16, 16], "shift": [0.5, 0.5, 0.5], "cell_atoms": 2,
        "supercell_grid": [8, 8, 8], "supercell_shift": [0.5, 0.5, 0.5],
        "supercell_factor": 2, "folding_verified": True,
        "n3V": 16**3 * float(B2.volume),
    },
    "phonon_density": 16,
    "smearing": {"ISMEAR": 1, "SIGMA": 0.20},
    "displacement_angstrom": 0.01,
    "potcar": {"functional": "PBE", "variants": {"Nb": "Nb_pv", "default": "standard"}},
    "magnetic_policy": "nonmagnetic",
}


def recipe(**overrides) -> Recipe:
    return Recipe.from_dict({**BASE_RECIPE, **overrides})


# --------------------------------------------------------------------------
# the k-density invariant
# --------------------------------------------------------------------------


@pytest.mark.parametrize("density", [8, 12, 16, 20])
def test_mesh_density_invariant_holds_across_a_supercell(density):
    """n^3 * V is the invariant, and fold-matching is exactly what preserves it.

    A 2x2x2 supercell has 8x the volume, so it needs half the mesh in each
    direction to sample the same density -- which is also the fold-matching
    condition.  If these ever disagree, the supercell and the cell are
    sampling differently.
    """

    cell = mesh_for(float(B2.volume), density)
    supercell = mesh_for(float(B2.volume) * 8, density)
    assert cell == 2 * supercell
    assert cell % 2 == 0 and supercell % 2 == 0
    assert cell**3 * float(B2.volume) == pytest.approx(
        supercell**3 * float(B2.volume) * 8, rel=0.0, abs=cell**3 * float(B2.volume) * 0.02
    )


def test_mesh_uses_the_volume_not_the_lattice_constant():
    """A non-cubic cell must not be under-sampled by using `a`."""

    hexagonal = Structure(
        Lattice.hexagonal(3.2, 5.1), ["Al", "Ta"], [[0, 0, 0], [1 / 3, 2 / 3, 0.5]]
    )
    by_volume = mesh_for(float(hexagonal.volume), 16)
    by_lattice_constant = mesh_for(float(hexagonal.lattice.a) ** 3, 16)
    assert by_volume != by_lattice_constant


@pytest.mark.parametrize(
    ("value", "expected"),
    [(16.0, 16), (15.0, 16), (17.0, 16), (1.0, 2), (0.5, 2), (100.0, 100)],
)
def test_round_even(value, expected):
    assert round_even(value) == expected


def test_mesh_is_always_even_and_at_least_two():
    """Even is not cosmetic: only an even mesh takes shift 0.5 and still
    contains Gamma, which is what makes it fold-matchable."""

    for volume in (5.0, 16.0, 33.5, 200.0, 1000.0):
        assert mesh_for(volume, 16) % 2 == 0
        assert mesh_for(volume, 16) >= 2


# --------------------------------------------------------------------------
# POTCAR
# --------------------------------------------------------------------------


def test_nb_uses_the_recorded_variant():
    """potpaw_PBE has no plain Nb -- Potcar(["Nb"]) raises FileNotFoundError."""

    assert potcar_symbols(ALNB, recipe()) == ["Al", "Nb_pv"]


def test_other_elements_use_the_plain_symbol():
    assert potcar_symbols(B2, recipe()) == ["Al", "Ta"]


def test_potcar_symbols_follow_poscar_order():
    reversed_cell = Structure(
        Lattice.cubic(3.3176), ["Nb", "Al"], [[0, 0, 0], [0.5, 0.5, 0.5]]
    )
    assert potcar_symbols(reversed_cell, recipe()) == ["Nb_pv", "Al"]


# --------------------------------------------------------------------------
# magnetism is not guessed
# --------------------------------------------------------------------------


def test_unresolved_policy_refuses_to_produce_moments():
    with pytest.raises(InputSetError, match="unresolved"):
        magmom_for(B2, recipe(magnetic_policy="unresolved"))


def test_unresolved_policy_refuses_to_produce_inputs(tmp_path):
    with pytest.raises(InputSetError, match="unresolved"):
        write_stage_inputs(tmp_path / "x", B2, recipe(magnetic_policy="unresolved"), "static")


def test_ispin_is_written_explicitly_for_both_policies():
    """Left unset, ISPIN defaults to 1 and MAGMOM is silently ignored."""

    assert incar_for("static", recipe(magnetic_policy="nonmagnetic"))["ISPIN"] == 1
    assert incar_for("static", recipe(magnetic_policy="confirmed_branch"))["ISPIN"] == 2


def test_magmom_is_written_only_when_spin_polarised(tmp_path):
    nonmagnetic = write_stage_inputs(
        tmp_path / "nm", B2, recipe(magnetic_policy="nonmagnetic"), "static"
    )
    magnetic = write_stage_inputs(
        tmp_path / "m", B2, recipe(magnetic_policy="confirmed_branch"), "static"
    )
    assert nonmagnetic.ok and magnetic.ok
    assert "MAGMOM" not in (tmp_path / "nm" / "INCAR").read_text()
    assert "MAGMOM" in (tmp_path / "m" / "INCAR").read_text()


# --------------------------------------------------------------------------
# the two electronic families are not interchangeable
# --------------------------------------------------------------------------


def test_phonon_uses_the_dfpt_electronic_recipe():
    """DFPT ran with EDIFF=1e-5 and KPAR=2, the statics with 1e-6 and NCORE=8.

    These were transcribed from the INCARs that actually ran.  Using the
    static settings for the phonon stage would silently change the convergence
    of every frequency.
    """

    phonon = incar_for("phonon", recipe(magnetic_policy="confirmed_branch"))
    static = incar_for("static", recipe(magnetic_policy="confirmed_branch"))
    assert phonon["EDIFF"] == DFPT_ELECTRONIC["EDIFF"] == 1e-5
    assert static["EDIFF"] == STATIC_ELECTRONIC["EDIFF"] == 1e-6
    assert phonon["IBRION"] == 8
    # NCORE and KPAR are mutually exclusive in VASP; the gate rejects both.
    assert "KPAR" in phonon and "NCORE" not in phonon
    assert "NCORE" in static and "KPAR" not in static


@pytest.mark.parametrize("stage", STAGES)
def test_every_stage_passes_the_gate(stage, tmp_path):
    is_cubic = stage != "shape"
    result = write_stage_inputs(
        tmp_path / stage, B2, recipe(magnetic_policy="confirmed_branch"), stage,
        is_cubic=is_cubic,
    )
    assert result.ok, list(result.problems)
    assert result.warnings == (), list(result.warnings)


def test_shape_stage_is_refused_for_a_cubic_cell():
    with pytest.raises(InputSetError, match="non-cubic"):
        incar_for("shape", recipe(), is_cubic=True)


# --------------------------------------------------------------------------
# VASP's n*value syntax, which pymatgen emits
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("2*0.0", [0.0, 0.0]),
        ("1.0 3.0", [1.0, 3.0]),
        ("2*0.0 1.0 3*2.0", [0.0, 0.0, 1.0, 2.0, 2.0, 2.0]),
        ("2*abc", None),
        ("0*1.0", None),
        ("2*0.0 extra", None),
    ],
)
def test_parse_repeated_floats(text, expected):
    assert parse_repeated_floats(text) == expected


def test_repeat_syntax_is_counted_not_taken_as_one_value(tmp_path):
    """`2*0.0` is TWO moments, not one.

    Reading it as a single token both reports a bogus non-numeric value and
    makes a wrong-length MAGMOM look right -- the exact failure that killed 14
    cases when a 16-atom cell got the 2-atom MAGMOM.
    """

    directory = tmp_path / "case"
    directory.mkdir()
    (directory / "POSCAR").write_text(
        "B2 AlTa\n1.0\n3.2246 0 0\n0 3.2246 0\n0 0 3.2246\nAl Ta\n1 1\nDirect\n"
        "0 0 0\n0.5 0.5 0.5\n"
    )
    (directory / "KPOINTS").write_text(
        "Automatic\n0\nMonkhorst-Pack\n16 16 16\n0.5 0.5 0.5\n"
    )
    (directory / "POTCAR").write_text(
        "   TITEL  = PAW_PBE Al 04Jan2001\n   ZVAL   =    3.0\n"
        "   TITEL  = PAW_PBE Ta 17Jan2003\n   ZVAL   =    5.0\n"
    )
    incar = directory / "INCAR"
    incar.write_text("ISMEAR = 1\nSIGMA = 0.2\nISPIN = 2\nMAGMOM = 2*0.0\n")
    assert validate_dir(directory).ok, list(validate_dir(directory).problems)

    # and a run of the wrong length must still be caught
    incar.write_text("ISMEAR = 1\nSIGMA = 0.2\nISPIN = 2\nMAGMOM = 3*0.0\n")
    result = validate_dir(directory)
    assert not result.ok
    assert any("MAGMOM has 3 values" in problem for problem in result.problems), result.problems


# --------------------------------------------------------------------------
# against the real pool
# --------------------------------------------------------------------------


POOL_STRUCTURES = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_structures.json")


@pytest.mark.skipif(
    not (POOL_RECIPES.is_file() and POOL_STRUCTURES.is_file()),
    reason="pool not built",
)
def test_representative_pool_endpoints_yield_gate_clean_inputs(tmp_path):
    """Real recipes against the real cells, one per class plus the Nb case.

    The recipe's k-mesh and cell_atoms are specific to its own cell, so this
    has to use the real structure -- a stand-in cell would fail the gate for
    reasons that say nothing about the generator.
    """

    recipes = json.loads(POOL_RECIPES.read_text())["endpoints"]
    structures = json.loads(POOL_STRUCTURES.read_text())["endpoints"]

    wanted = ["BCC-Al-Ta-2", "FCC-Al-Cr-4", "SER-Zr-1", "SER-Al-1", "BCC-Al-Nb-2"]
    missing = [name for name in wanted if name not in recipes]
    assert not missing, missing

    for name in wanted:
        structure = Structure.from_dict(structures[name]["structure"])
        cell = structure.lattice
        is_cubic = pytest.approx(cell.angles, abs=1.0) == (90.0, 90.0, 90.0)
        parsed = Recipe.from_dict(
            {**recipes[name]["recipe"], "magnetic_policy": "confirmed_branch"}
        )
        stages = ("static", "phonon") if is_cubic else ("static", "phonon", "shape")
        for stage in stages:
            result = write_stage_inputs(
                tmp_path / name / stage, structure, parsed, stage, is_cubic=is_cubic
            )
            assert result.ok, (name, stage, list(result.problems))
            assert result.warnings == (), (name, stage, list(result.warnings))

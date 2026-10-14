"""Tests for the DFPT-QHA recipe contract.

Most of what follows is a negative case.  The contract's whole value is that it
FAILS when a field is absent or contradicts the structure, and a test that only
checks "a good recipe loads" passes just as well against a module that ignores
every field.

The real 84-endpoint pool is loaded at the end: a schema that only works on
hand-written fixtures has not been tested against the data it exists for.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from pymatgen.core import Lattice, Structure

from htvasp.dfpt_qha.recipe import (
    EndpointRecord,
    MissingRecipeField,
    Recipe,
    RecipeInconsistent,
    Validation,
    axis_contains_gamma,
    gamma_centred,
    load_pair,
)

POOL_RECIPES = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")
POOL_STRUCTURES = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_structures.json")

B2 = Structure(
    Lattice.cubic(3.2246), ["Al", "Ta"], [[0, 0, 0], [0.5, 0.5, 0.5]]
)
L12 = Structure(
    Lattice.cubic(4.0389),
    ["Al", "Cr", "Cr", "Cr"],
    [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
)

GOOD_RECIPE = {
    "structure_class": "BCC_B2",
    "site_occupancy": {"1a": "Al", "1b": "Ta"},
    "representation": "conventional",
    "sym_reduce": True,
    "supercell": [2, 2, 2],
    "kpoints": {
        "grid": [16, 16, 16],
        "shift": [0.5, 0.5, 0.5],
        "cell_atoms": 2,
        "supercell_grid": [8, 8, 8],
        "supercell_shift": [0.5, 0.5, 0.5],
        "supercell_factor": 2,
        "folding_verified": True,
        "n3V": 16**3 * float(B2.volume),
    },
    "phonon_density": 16,
    "smearing": {"ISMEAR": 1, "SIGMA": 0.20},
    "displacement_angstrom": 0.01,
    "potcar": {"functional": "PBE", "variants": {"Nb": "Nb_pv", "default": "standard"}},
    "magnetic_policy": "nonmagnetic",
}

GOOD_VALIDATION = {
    "stage0_coverage": "measured",
    "measured_on": ["BCC-Al-Ta-2"],
    "not_measured": [],
    "qha_applicability": "not_yet_assessed",
    "recipe_robustness": "phonon_recipe_robust",
    "recorded_at": "2026-10-13",
}

TOP_LEVEL_FIELDS = (
    "structure_class",
    "site_occupancy",
    "representation",
    "sym_reduce",
    "supercell",
    "kpoints",
    "phonon_density",
    "smearing",
    "displacement_angstrom",
    "potcar",
    "magnetic_policy",
)

KPOINT_FIELDS = (
    "grid",
    "shift",
    "cell_atoms",
    "supercell_grid",
    "supercell_shift",
    "supercell_factor",
    "folding_verified",
    "n3V",
)

VALIDATION_FIELDS = (
    "stage0_coverage",
    "measured_on",
    "not_measured",
    "qha_applicability",
    "recipe_robustness",
    "recorded_at",
)


def variant(recipe=None, validation=None) -> dict:
    """A deep copy every time.

    ``dict(GOOD_RECIPE, **overrides)`` is a SHALLOW copy, so nested dicts stay
    shared with the module-level fixture.  The negative tests below delete and
    overwrite keys inside those nested dicts, so with a shallow copy the first
    test to touch ``kpoints`` mutates it for every later test -- which is
    exactly what happened when this suite lived in scripts/, and it read as a
    module bug rather than a test bug.
    """

    merged = copy.deepcopy(GOOD_RECIPE)
    merged.update(copy.deepcopy(recipe or {}))
    valid = copy.deepcopy(GOOD_VALIDATION)
    valid.update(copy.deepcopy(validation or {}))
    return {"recipe": merged, "validation": valid}


def l12_recipe(**overrides) -> dict:
    base = {
        "structure_class": "FCC_L12",
        "site_occupancy": {"1a": "Al", "3c": ["Cr", "Cr", "Cr"]},
        "kpoints": dict(
            GOOD_RECIPE["kpoints"],
            cell_atoms=4,
            n3V=16**3 * float(L12.volume),
        ),
    }
    base.update(overrides)
    return base


# --------------------------------------------------------------------------
# absent fields must raise, and must name every one of them
# --------------------------------------------------------------------------


@pytest.mark.parametrize("field", TOP_LEVEL_FIELDS)
def test_missing_recipe_field_is_named(field):
    data = variant()["recipe"]
    del data[field]
    with pytest.raises(MissingRecipeField) as excinfo:
        Recipe.from_dict(data)
    assert field in excinfo.value.missing


@pytest.mark.parametrize("field", KPOINT_FIELDS)
def test_missing_kpoint_field_is_fully_qualified(field):
    data = variant()["recipe"]
    del data["kpoints"][field]
    with pytest.raises(MissingRecipeField) as excinfo:
        Recipe.from_dict(data)
    assert f"kpoints.{field}" in excinfo.value.missing


@pytest.mark.parametrize("field", VALIDATION_FIELDS)
def test_missing_validation_field_is_named(field):
    data = variant()["validation"]
    del data[field]
    with pytest.raises(MissingRecipeField) as excinfo:
        Validation.from_dict(data)
    assert field in excinfo.value.missing


def test_every_gap_is_reported_at_once():
    """A builder missing four fields should learn that in one round trip."""

    data = variant()["recipe"]
    for field in ("sym_reduce", "magnetic_policy", "phonon_density"):
        del data[field]
    del data["kpoints"]["grid"]
    with pytest.raises(MissingRecipeField) as excinfo:
        Recipe.from_dict(data)
    assert set(excinfo.value.missing) == {
        "sym_reduce",
        "magnetic_policy",
        "phonon_density",
        "kpoints.grid",
    }


def test_load_pair_needs_both_layers():
    with pytest.raises(MissingRecipeField) as excinfo:
        load_pair({"recipe": GOOD_RECIPE})
    assert "validation" in excinfo.value.missing


# --------------------------------------------------------------------------
# a recipe that parses must still be internally consistent
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("label", "override", "needle"),
    [
        (
            "supercell disagrees with supercell_factor",
            {"supercell": [2, 3, 2]},
            "supercell_factor",
        ),
        (
            "folding pairing broken",
            {"kpoints": dict(GOOD_RECIPE["kpoints"], supercell_grid=[7, 8, 8])},
            "fold-matched",
        ),
        (
            "cell mesh excludes Gamma",
            {"kpoints": dict(GOOD_RECIPE["kpoints"], shift=[0.0, 0.5, 0.5])},
            "kpoints[cell]",
        ),
        (
            "supercell mesh excludes Gamma",
            {"kpoints": dict(GOOD_RECIPE["kpoints"], supercell_shift=[0.0, 0.5, 0.5])},
            "kpoints[supercell]",
        ),
        (
            "smearing is not the frozen value",
            {"smearing": {"ISMEAR": -5, "SIGMA": 0.20}},
            "frozen",
        ),
        (
            "displacement is not the frozen value",
            {"displacement_angstrom": 0.05},
            "0.01",
        ),
        (
            "POTCAR variant contradicts the frozen choice",
            {"potcar": {"functional": "PBE", "variants": {"Nb": "Nb_sv"}}},
            "Nb_pv",
        ),
        (
            "POTCAR invents an exception",
            {"potcar": {"functional": "PBE", "variants": {"Ta": "Ta_pv"}}},
            "no recorded exception",
        ),
        (
            "POTCAR functional is not PBE",
            {"potcar": {"functional": "LDA", "variants": {"default": "standard"}}},
            "PBE",
        ),
        ("sym_reduce written as a string", {"sym_reduce": "true"}, "boolean"),
        ("structure_class outside the vocabulary", {"structure_class": "B2"}, "BCC_B2"),
        ("an unknown key", {"nonsense": 1}, "Extra inputs"),
    ],
)
def test_inconsistent_recipe_is_rejected(label, override, needle):
    with pytest.raises(RecipeInconsistent) as excinfo:
        Recipe.from_dict(variant(override)["recipe"])
    assert needle in str(excinfo.value), label


def test_gamma_truth_table():
    """Gamma is on an axis iff n is odd with shift 0, or even with shift 0.5."""

    assert axis_contains_gamma(7, 0.0) is True
    assert axis_contains_gamma(7, 0.5) is False
    assert axis_contains_gamma(8, 0.0) is False
    assert axis_contains_gamma(8, 0.5) is True
    # Gamma is in the MESH only when it is on every axis, so one bad axis is
    # enough -- the mixed-shift case that a whole-vector rule would miss.
    assert gamma_centred((16, 16, 16), (0.5, 0.5, 0.5)) is True
    assert gamma_centred((16, 16, 16), (0.0, 0.5, 0.5)) is False
    assert gamma_centred((8, 8, 8), (0.5, 0.5, 0.5)) is True


# --------------------------------------------------------------------------
# recipe vs the structure it claims to describe
# --------------------------------------------------------------------------


def test_good_recipe_describes_its_structure():
    recipe = Recipe.from_dict(GOOD_RECIPE)
    assert recipe.check_against(B2) == []


def test_recipe_from_dict_round_trips():
    recipe = Recipe.from_dict(GOOD_RECIPE)
    assert Recipe.from_dict(recipe.model_dump()) == recipe


@pytest.mark.parametrize(
    ("override", "structure", "needle"),
    [
        ({"kpoints": dict(GOOD_RECIPE["kpoints"], n3V=1.0)}, B2, "n3V"),
        ({"kpoints": dict(GOOD_RECIPE["kpoints"], cell_atoms=4)}, B2, "cell_atoms"),
        (
            {"site_occupancy": {"1a": "Cr", "1b": "Ta"}},
            B2,
            "site_occupancy",
        ),
    ],
)
def test_recipe_contradicting_its_structure(override, structure, needle):
    recipe = Recipe.from_dict(variant(override)["recipe"])
    problems = recipe.check_against(structure)
    assert any(needle in problem for problem in problems), problems


def test_class_and_cell_shape_must_agree():
    hexagonal = Structure(
        Lattice.hexagonal(3.2, 5.1), ["Al", "Ta"], [[0, 0, 0], [1 / 3, 2 / 3, 0.5]]
    )
    recipe = Recipe.from_dict(GOOD_RECIPE)
    assert any("expects a cubic" in problem for problem in recipe.check_against(hexagonal))


def test_l12_direction_is_checked_against_the_cell():
    recipe = Recipe.from_dict(variant(l12_recipe())["recipe"])
    assert recipe.check_against(L12) == []

    reversed_l12 = Structure(
        Lattice.cubic(4.0389),
        ["Cr", "Al", "Al", "Al"],
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
    )
    assert any(
        "1a" in problem for problem in recipe.check_against(reversed_l12)
    ), "the two L1_2 directions are different compounds"


def test_l12_requires_an_occupancy():
    """The class alone does not say which L1_2 direction this is."""

    data = variant(l12_recipe())["recipe"]
    data["site_occupancy"] = {}
    with pytest.raises(RecipeInconsistent) as excinfo:
        Recipe.from_dict(data)
    assert "site_occupancy" in str(excinfo.value)


# --------------------------------------------------------------------------
# the two validation axes must be free to disagree
# --------------------------------------------------------------------------


def test_axes_are_independent():
    record = Validation.from_dict(
        dict(
            GOOD_VALIDATION,
            qha_applicability="finite_q_soft_mode_diagnostic",
            recipe_robustness="phonon_recipe_robust",
        )
    )
    assert record.qha_applicability != "harmonic_QHA_applicable"
    assert record.recipe_robustness == "phonon_recipe_robust"


@pytest.mark.parametrize(
    ("override", "needle"),
    [
        ({"measured_on": []}, "empty measured_on"),
        ({"stage0_coverage": "untested", "measured_on": ["X"]}, "untested"),
        ({"measured_on": ["X", "X"]}, "duplicates"),
        ({"qha_applicability": "maybe"}, "qha_applicability"),
    ],
)
def test_validation_inconsistencies(override, needle):
    with pytest.raises(RecipeInconsistent) as excinfo:
        Validation.from_dict(dict(GOOD_VALIDATION, **override))
    assert needle in str(excinfo.value)


def test_untested_class_cannot_claim_robustness_and_applicability():
    with pytest.raises(RecipeInconsistent) as excinfo:
        Validation.from_dict(
            dict(
                GOOD_VALIDATION,
                stage0_coverage="untested",
                measured_on=[],
                qha_applicability="harmonic_QHA_applicable",
            )
        )
    assert "robustness claim needs a measurement" in str(excinfo.value)


# --------------------------------------------------------------------------
# the real pool
# --------------------------------------------------------------------------


@pytest.mark.skipif(not POOL_RECIPES.is_file(), reason="pool_recipes.json not built")
def test_real_pool_recipes_are_self_consistent():
    pool = json.loads(POOL_RECIPES.read_text())
    assert pool["n_endpoints"] == 84, pool["n_endpoints"]
    assert pool["n_problems"] == 0, pool["problems"][:3]

    by_class: dict[str, int] = {}
    for name, record in pool["endpoints"].items():
        parsed = EndpointRecord.from_dict(record)
        assert parsed.n_atoms > 0, name
        by_class[parsed.recipe.structure_class] = (
            by_class.get(parsed.recipe.structure_class, 0) + 1
        )
    assert sum(by_class.values()) == 84, by_class
    # the pool is 28 BCC/B2 + 49 FCC/L1_2 + 7 SER
    assert by_class.get("BCC_B2") == 28, by_class
    assert by_class.get("FCC_L12") == 49, by_class
    assert sum(v for k, v in by_class.items() if k.startswith("SER_")) == 7, by_class


@pytest.mark.skipif(
    not (POOL_RECIPES.is_file() and POOL_STRUCTURES.is_file()),
    reason="pool not built",
)
def test_real_pool_recipe_matches_the_real_structure():
    """The check that matters: 84 recipes against the 84 cells they describe."""

    recipes = json.loads(POOL_RECIPES.read_text())["endpoints"]
    structures = json.loads(POOL_STRUCTURES.read_text())["endpoints"]
    assert set(recipes) == set(structures), "the two pool files describe different sets"

    for name, record in recipes.items():
        structure = Structure.from_dict(structures[name]["structure"])
        parsed = EndpointRecord.from_dict(record)
        assert parsed.n_atoms == len(structure), name
        assert parsed.recipe.check_against(structure) == [], name


@pytest.mark.skipif(not POOL_RECIPES.is_file(), reason="pool_recipes.json not built")
def test_real_pool_round_trips():
    pool = json.loads(POOL_RECIPES.read_text())
    for name, record in pool["endpoints"].items():
        parsed = EndpointRecord.from_dict(record)
        assert EndpointRecord.from_dict(parsed.model_dump()) == parsed, name

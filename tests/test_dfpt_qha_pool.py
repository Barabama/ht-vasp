"""Every pool endpoint must generate a valid DFPT input set.

This is the batch's pre-flight, turned into a test: 84 endpoints, each loaded
from the pool's recipe and structure, each run through the same generators the
production flow uses.  A broken MAGMOM table, a POTCAR variant that does not
resolve, or a k-mesh that does not fold onto the supercell fails here rather
than after an allocation has been spent.

The pool is vendored under ``tests/fixtures/pool`` as a snapshot so the test is
hermetic; the live pool is ``/nfs_ssd/tmp/legacy-qha-pool/stage0``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from htvasp.dfpt_qha.input_set import DfptInputSetGenerator, DfptRelaxSetGenerator
from htvasp.dfpt_qha.phonons import supercell_structure
from htvasp.model.dfpt_qha_recipe import Recipe

POOL = Path(__file__).parent / "fixtures" / "pool"
SUPERCELL = ((2, 0, 0), (0, 2, 0), (0, 0, 2))


def _endpoints():
    from pymatgen.core import Structure

    recipes = json.loads((POOL / "pool_recipes.json").read_text())["endpoints"]
    structures = json.loads((POOL / "pool_structures.json").read_text())["endpoints"]
    return [
        (name, Recipe.from_dict(recipes[name]["recipe"]),
         Structure.from_dict(structures[name]["structure"]))
        for name in sorted(recipes)
    ]


def test_the_pool_has_84_endpoints_and_no_unresolved_magnetism():
    endpoints = _endpoints()
    assert len(endpoints) == 84
    unresolved = [name for name, recipe, _ in endpoints
                  if recipe.magnetic_policy == "unresolved"]
    assert unresolved == [], f"still unresolved: {unresolved}"


@pytest.mark.parametrize(
    ("name", "recipe", "structure"),
    _endpoints(),
    ids=[name for name, _, _ in _endpoints()],
)
def test_every_endpoint_generates_a_dfpt_input_set(name, recipe, structure):
    cell = supercell_structure(structure, SUPERCELL)
    generator = DfptInputSetGenerator(recipe=recipe, use_supercell_kpoints=True)
    generator.structure = cell

    incar = generator.incar
    assert incar["IBRION"] == 8
    assert incar["NSW"] == 1
    assert incar["ENCUT"] == pytest.approx(520.0)
    assert incar["GGA"].upper().startswith("PE")
    # MAGMOM must arrive expanded per site, not as a dict
    assert not isinstance(incar["MAGMOM"], dict)
    assert len(incar["MAGMOM"]) == len(cell)

    # the supercell mesh is the fold of the cell mesh
    kpts = generator.kpoints.kpts[0]
    assert tuple(kpts) == tuple(recipe.kpoints.supercell_grid)
    for axis in range(3):
        assert kpts[axis] * recipe.kpoints.supercell_factor == recipe.kpoints.grid[axis]

    assert generator.potcar_symbols


@pytest.mark.parametrize(
    ("name", "recipe", "structure"),
    _endpoints(),
    ids=[name for name, _, _ in _endpoints()],
)
def test_every_endpoint_generates_a_relax_input_set(name, recipe, structure):
    generator = DfptRelaxSetGenerator(recipe=recipe, isif=2)
    generator.structure = structure
    incar = generator.incar
    assert incar["IBRION"] == 2
    assert incar["ISIF"] == 2
    assert incar["NSW"] > 0
    assert np.all(np.asarray(incar["MAGMOM"], dtype=float) > -1.0)

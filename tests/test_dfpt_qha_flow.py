"""No-VASP graph construction tests for the DFPT volume flow."""

from __future__ import annotations

import json
from pathlib import Path

from jobflow import Flow
from pymatgen.core import Lattice, Structure

from htvasp.dfpt_qha.maker import DfptQhaMaker
from htvasp.model.dfpt_qha_recipe import Recipe

RECIPE_FILE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")
STRUCTURE = Structure(Lattice.cubic(3.2246), ["Al", "Ta"], [[0, 0, 0], [0.5, 0.5, 0.5]])


def recipe() -> Recipe:
    data = json.loads(RECIPE_FILE.read_text())["endpoints"]["BCC-Al-Ta-2"]["recipe"]
    return Recipe.from_dict(data)


def test_dfpt_qha_maker_does_not_use_finite_displacement_phonon_maker():
    maker = DfptQhaMaker(recipe=recipe())
    flow = maker.make(STRUCTURE, supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    assert isinstance(flow, Flow)
    assert [job.name for job in flow.jobs] == ["EOS Maker", "spawn_dfpt_volume_jobs"]
    assert maker.skip_analysis is True


def test_dfpt_qha_maker_requires_recipe():
    try:
        DfptQhaMaker()
    except ValueError as exc:
        assert "Recipe" in str(exc)
    else:
        raise AssertionError("DfptQhaMaker accepted no recipe")

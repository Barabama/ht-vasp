"""No-VASP graph construction tests for the DFPT volume flow."""

from __future__ import annotations

import json
from pathlib import Path

from jobflow import Flow, JobStore
from maggma.stores import MemoryStore
from pymatgen.core import Lattice, Structure

from htvasp.dfpt_qha.flow_jobs import spawn_dfpt_volume_jobs
from htvasp.dfpt_qha.input_set import DfptInputSetGenerator, DfptRelaxSetGenerator
from htvasp.dfpt_qha.maker import DfptQhaMaker
from htvasp.model.dfpt_qha_recipe import Recipe

RECIPE_FILE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")
STRUCTURE_FILE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_structures.json")
ENDPOINT = "BCC-Al-Ta-2"


def recipe() -> Recipe:
    data = json.loads(RECIPE_FILE.read_text())["endpoints"][ENDPOINT]["recipe"]
    return Recipe.from_dict(data)


def structure() -> Structure:
    data = json.loads(STRUCTURE_FILE.read_text())["endpoints"][ENDPOINT]["structure"]
    return Structure.from_dict(data)


def test_dfpt_qha_maker_builds_eos_then_the_dfpt_chain():
    maker = DfptQhaMaker(recipe=recipe())
    flow = maker.make(structure(), supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    assert isinstance(flow, Flow)
    assert [job.name for job in flow.jobs] == ["EOS Maker", "spawn_dfpt_volume_jobs"]
    assert maker.skip_analysis is False


def test_dfpt_qha_maker_requires_recipe():
    try:
        DfptQhaMaker()
    except ValueError as exc:
        assert "Recipe" in str(exc)
    else:
        raise AssertionError("DfptQhaMaker accepted no recipe")


def test_eos_stage_relaxes_instead_of_running_dfpt():
    """The EOS stage must move ions; IBRION=8 cannot.

    Wiring the DFPT generator into the EOS stage builds the E(V) curve from
    unrelaxed strained geometries, because VASP ignores ISIF/NSW for a
    dynamical-matrix run.
    """

    maker = DfptQhaMaker(recipe=recipe())
    initial = maker.initial_relax_maker.input_set_generator
    eos = maker.eos_relax_maker.input_set_generator
    assert isinstance(initial, DfptRelaxSetGenerator)
    assert isinstance(eos, DfptRelaxSetGenerator)
    assert initial.incar_updates["IBRION"] == 2
    assert initial.incar_updates["ISIF"] == 3
    assert eos.incar_updates["IBRION"] == 2
    assert eos.incar_updates["ISIF"] == 2
    assert eos.incar_updates["NSW"] > 0


def test_dfpt_stage_runs_on_the_recipe_supercell():
    """A DFPT run on the 2-atom cell is Gamma-only and cannot see the
    finite-q instability the recipe exists to characterise."""

    job = spawn_dfpt_volume_jobs(
        {"relax": {"structure": [structure()], "dir_name": ["/tmp/eos0"],
                   "energy": [-10.0], "volume": [float(structure().volume)]}},
        recipe=recipe().model_dump(mode="json"),
    )
    store = MemoryStore()
    store.connect()
    flow = job.run(JobStore(store)).replace
    vasp_jobs = [j for j in flow.jobs if j.name.startswith("DFPT volume")]
    assert len(vasp_jobs) == 1
    cell = vasp_jobs[0].function_args[0]
    assert len(cell) == 16, cell.formula
    generator = vasp_jobs[0].maker.input_set_generator
    generator.structure = cell
    assert generator.incar.get("IBRION") == 8
    # fold-matched supercell mesh, not the cell mesh
    assert tuple(generator.kpoints.kpts[0]) == tuple(recipe().kpoints.supercell_grid)


def test_dfpt_cell_is_phonopy_ordered_not_pymatgen_ordered():
    """The DFPT cell must be phonopy's supercell.

    pymatgen's make_supercell gives the same crystal in a different atom order.
    Force constants come back indexed by whatever order VASP read, so running
    on pymatgen's cell would attach them to the wrong atoms -- still nonzero,
    still symmetric, still plausible frequencies.
    """

    from htvasp.dfpt_qha.phonons import (
        ordering_fingerprint,
        supercell_structure,
    )

    cell = structure()
    matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    phonopy_cell = supercell_structure(cell, matrix)
    pymatgen_cell = cell.make_supercell(matrix, in_place=False)
    assert len(phonopy_cell) == len(pymatgen_cell) == 16
    assert ordering_fingerprint(phonopy_cell) != ordering_fingerprint(pymatgen_cell)


def test_flow_ends_in_the_qha_analysis_and_gates_on_validation():
    """The fail-closed gate must have something downstream to stop."""

    job = spawn_dfpt_volume_jobs(
        {"relax": {"structure": [structure()], "dir_name": ["/tmp/eos0"],
                   "energy": [-10.0], "volume": [float(structure().volume)]}},
        recipe=recipe().model_dump(mode="json"),
    )
    store = MemoryStore()
    store.connect()
    flow = job.run(JobStore(store)).replace
    names = [j.name for j in flow.jobs]
    assert "collect validated volumes" in names
    assert "DFPT thermal properties" in names
    assert "analyze DFPT-QHA" in names
    # the analysis is the flow's terminal job, so a stopped volume really does
    # halt the thermodynamics rather than merely reporting ok=False
    assert names.index("analyze DFPT-QHA") == len(names) - 2  # trailing store_inputs
    assert names.index("collect validated volumes") > names.index(
        "validate DFPT volume 1/1"
    )


def test_supercell_kpoints_are_the_fold_matched_ones():
    cell_generator = DfptInputSetGenerator(recipe=recipe(), use_supercell_kpoints=False)
    supercell_generator = DfptInputSetGenerator(recipe=recipe(), use_supercell_kpoints=True)
    cell_generator.structure = structure()
    supercell_generator.structure = structure()
    assert tuple(cell_generator.kpoints.kpts[0]) == tuple(recipe().kpoints.grid)
    assert tuple(supercell_generator.kpoints.kpts[0]) == tuple(
        recipe().kpoints.supercell_grid
    )
    assert (
        supercell_generator.kpoints.kpts[0][0] * recipe().kpoints.supercell_factor
        == cell_generator.kpoints.kpts[0][0]
    )


def test_potcar_symbols_are_resolved_by_atomate2_not_hardcoded():
    """No explicit variant table is set, so atomate2/pymatgen resolve them.

    Hard-coding was how a dead ``"default": "standard"`` key produced the wrong
    variant for every endpoint without any error.
    """

    generator = DfptInputSetGenerator(recipe=recipe())
    generator.structure = structure()
    assert generator.user_potcar_settings == {}
    assert generator.potcar_symbols == ["Al", "Ta_pv"]

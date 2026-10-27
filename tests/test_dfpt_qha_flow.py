"""No-VASP graph construction tests for the DFPT volume flow."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from jobflow import Flow, JobStore
from maggma.stores import MemoryStore
from pymatgen.core import Lattice, Structure

from htvasp.dfpt_qha.flow_jobs import (
    collect_validated_volumes,
    spawn_dfpt_volume_jobs,
)
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


def test_generators_build_a_real_incar():
    """``.incar``, not just ``incar_updates``.

    The first real Slurm run died here: pymatgen expands a per-element MAGMOM
    dict into a per-site list while building the INCAR from the settings dict,
    and it only consults ``user_incar_settings`` to do it.  A dict left in
    ``incar_updates`` -- which pymatgen applies afterwards -- reached
    ``Incar(...)`` unexpanded and raised KeyError: 0.  Asserting only
    ``incar_updates`` (as the earlier tests did) cannot see that.
    """

    cell = structure()
    for generator, ibrion in (
        (DfptRelaxSetGenerator(recipe=recipe(), isif=3), 2),
        (DfptRelaxSetGenerator(recipe=recipe(), isif=2), 2),
        (DfptInputSetGenerator(recipe=recipe()), 8),
    ):
        generator.structure = cell
        incar = generator.incar  # must not raise
        assert incar["IBRION"] == ibrion
        assert incar["ISPIN"] == 2
        # expanded per site, in POSCAR order: Al then Ta
        assert incar["MAGMOM"] == [1.0, 3.0], (type(generator).__name__, incar["MAGMOM"])
        assert incar["ISMEAR"] == 1 and abs(incar["SIGMA"] - 0.20) < 1e-12
        assert generator.potcar_symbols


def test_dfpt_generator_leaves_no_dict_valued_incar_key():
    """Any dict-valued setting other than a resolved one would crash pymatgen."""

    generator = DfptInputSetGenerator(recipe=recipe())
    generator.structure = structure()
    for key, value in generator.incar.items():
        assert not isinstance(value, dict), f"{key} is a dict: {value!r}"


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


def _spawned_flow():
    job = spawn_dfpt_volume_jobs(
        {"relax": {"structure": [structure()], "dir_name": ["/tmp/eos0"],
                   "energy": [-10.0], "volume": [float(structure().volume)]}},
        recipe=recipe().model_dump(mode="json"),
    )
    store = MemoryStore()
    store.connect()
    return job.run(JobStore(store)).replace


def test_stability_is_one_job_after_collection():
    """The verdict is per endpoint, so it runs after every volume is in.

    A soft volume fails the endpoint rather than being dropped from the fit: a
    hole in F(V,T) is invisible to the analysis, and the curve still looks
    smooth.  A per-volume job would also save nothing -- the VASP time is spent
    before either job runs.
    """

    names = [j.name for j in _spawned_flow().jobs]
    collected = names.index("collect validated volumes")
    classify = names.index("classify DFPT stability")
    thermal = names.index("DFPT thermal properties")
    analysis = names.index("analyze DFPT-QHA")
    assert collected < classify < thermal < analysis
    assert not any(name.startswith("classify DFPT volume") for name in names)


def test_validated_volumes_carry_where_the_matrix_came_from():
    """Provenance is part of the record, not a log line.

    The Hessian in vasprun.xml and the DYNMAT replay are both VASP output; a
    reader downstream has to be able to tell which one produced the frequencies.
    """

    document = {
        "force_constants": [],
        "force_constants_provenance": {
            "n_atoms": 16,
            "ordering_fingerprint": "x" * 64,
            "source": "vasprun.xml dynmat hessian",
        },
    }
    cells = [{"index": 0, "cell": {}, "supercell": {}, "volume_angstrom3": 1.0,
              "static_energy_eV": -1.0}]
    output = collect_validated_volumes(
        [{"ok": True, "dfpt_doc": document}], cells=cells
    ).run(JobStore(_memory())).output
    assert output["n_validated"] == 1
    assert output["volumes"][0]["force_constants_source"] == "vasprun.xml dynmat hessian"


def _memory():
    store = MemoryStore()
    store.connect()
    return store


def test_worker_falls_back_to_the_stability_verdict_when_the_flow_stopped():
    """A stopped flow still has to report why.

    Without the fallback the Worker returns ``None`` for every fail-closed
    endpoint -- i.e. exactly the ones whose reason matters.
    """

    from htvasp.workflows.dfpt_qha import DfptQhaWorker

    worker = DfptQhaWorker(recipe=recipe())
    store = JobStore(_memory())
    store.connect()
    payload = {
        "gate": "dfpt_stability",
        "ok": False,
        "qha_allowed": False,
        "n_volumes": 5,
        "n_blocked": 1,
        "blocked_volume_indices": [1],
        "reports": [
            {"index": 1, "status": "finite_q_soft_mode", "qha_allowed": False,
             "reason": "commensurate finite-q imaginary branch"},
        ],
    }
    store.update(
        {"uuid": "u0", "index": 7, "name": "classify DFPT stability",
         "output": payload},
        key="uuid",
        save=True,
    )
    worker.store = store

    result = worker.get_result()
    assert result is not None
    assert result["gate"] == "dfpt_stability"
    assert result["blocked_volume_indices"] == [1]
    assert result["reports"][0]["status"] == "finite_q_soft_mode"


def test_worker_reports_a_validation_failure_when_no_stability_job_exists():
    """A volume rejected before collection was previously reported as ``None``.

    That is the other fail-closed path -- wrong atom count, unusable shape --
    and it happens before the stability job exists.
    """

    from htvasp.workflows.dfpt_qha import DfptQhaWorker

    worker = DfptQhaWorker(recipe=recipe())
    store = JobStore(_memory())
    store.connect()
    store.update(
        {"uuid": "u0", "index": 5, "name": "validate DFPT volume 2/5",
         "output": {"ok": False,
                    "error": "DFPT ran on 2 atoms but the recipe's phonon "
                             "supercell is 16 atoms"}},
        key="uuid",
        save=True,
    )
    store.update(
        {"uuid": "u1", "index": 4, "name": "validate DFPT volume 1/5",
         "output": {"ok": True, "dfpt_doc": {}}},
        key="uuid",
        save=True,
    )
    worker.store = store

    result = worker.get_result()
    assert result is not None
    assert result["gate"] == "dfpt_validation"
    assert result["n_failed"] == 1
    assert "16 atoms" in result["failures"][0]["error"]


def test_the_thermal_job_consumes_the_stability_output():
    """Naming the jobs in the right order is not the same as wiring them.

    ``dfpt_thermal_properties(collected.output, ...)`` keeps the job names and
    their order, but bypasses the gate: every endpoint would then stop at the
    thermal job after all five VASP volumes had been paid for.
    """

    flow = _spawned_flow()
    stability = next(j for j in flow.jobs if j.name == "classify DFPT stability")
    thermal = next(j for j in flow.jobs if j.name == "DFPT thermal properties")
    source = thermal.function_args[0]
    assert getattr(source, "uuid", None) == stability.uuid


def test_worker_prefers_the_analysis_over_the_gate():
    """Both jobs exist on a successful run; the analysis is the result."""

    from htvasp.workflows.dfpt_qha import DfptQhaWorker

    worker = DfptQhaWorker(recipe=recipe())
    store = JobStore(_memory())
    store.connect()
    store.update(
        {"uuid": "a", "index": 7, "name": "classify DFPT stability",
         "output": {"gate": "dfpt_stability", "ok": True, "qha_allowed": True}},
        key="uuid", save=True,
    )
    store.update(
        {"uuid": "b", "index": 9, "name": "analyze DFPT-QHA",
         "output": {"ok": True, "eos": {"V0_angstrom3": 36.0}, "temperatures": [0.0]}},
        key="uuid", save=True,
    )
    worker.store = store
    result = worker.get_result()
    assert result["eos"]["V0_angstrom3"] == 36.0
    assert "gate" not in result


def test_worker_does_not_report_a_passed_gate_as_a_result():
    """A passed gate with no analysis means the flow died elsewhere.

    Returning the gate document there would be written to the endpoint's JSON,
    marked complete, and skipped forever -- with ok=true and no Gibbs curve.
    """

    from htvasp.workflows.dfpt_qha import DfptQhaWorker

    worker = DfptQhaWorker(recipe=recipe())
    store = JobStore(_memory())
    store.connect()
    store.update(
        {"uuid": "a", "index": 7, "name": "classify DFPT stability",
         "output": {"gate": "dfpt_stability", "ok": True, "qha_allowed": True}},
        key="uuid", save=True,
    )
    worker.store = store
    assert worker.get_result() is None


def test_the_magmom_table_reaches_the_child_generators():
    """The 16-element guess must be the one the DFPT INCAR is built from.

    A zero start on a Cr-containing cell converges to the non-magnetic
    solution, which is 0.28 eV/atom above the magnetic ground state; the table
    is only useful if it survives the Worker -> Maker -> generator chain.
    """

    from htvasp.dfpt_qha.maker import DfptQhaMaker
    from htvasp.workflows.dfpt_qha import DfptQhaWorker

    table = {"Al": 1.0, "Cr": 3.0}
    maker = DfptQhaMaker(recipe=recipe(), magmom=table)
    assert maker.initial_relax_maker.input_set_generator.magmom == table
    assert maker.eos_relax_maker.input_set_generator.magmom == table

    worker = DfptQhaWorker(recipe=recipe(), magmom=table)
    assert worker.flow_maker.magmom == table


def test_the_package_magmom_table_covers_all_sixteen_elements():
    from htvasp.dfpt_qha.input_set import MAGMOM_BY_ELEMENT

    assert len(MAGMOM_BY_ELEMENT) == 16
    for element in ("Al", "Cr", "Hf", "Nb", "Ta", "Ti", "V", "Zr"):
        assert element in MAGMOM_BY_ELEMENT


def test_incar_overrides_win_over_the_module_defaults():
    """A different project must not have to edit input_set.py to change ENCUT.

    The module's tables are a baseline for this pool; ``incar_overrides`` is
    applied last, after the frozen recipe and mode settings.
    """

    from htvasp.dfpt_qha.input_set import DfptInputSetGenerator

    generator = DfptInputSetGenerator(
        recipe=recipe(), incar_overrides={"ENCUT": 700, "NELM": 200}
    )
    generator.structure = structure()
    incar = generator.incar
    assert incar["ENCUT"] == pytest.approx(700.0)
    assert incar["NELM"] == 200
    # the mode defaults survive
    assert incar["IBRION"] == 8
    assert incar["NSW"] == 1


def test_ncore_and_kpar_reach_their_own_stage_and_not_the_other():
    """NCORE and KPAR are mutually exclusive, so they cannot share a dict.

    The relaxation family carries NCORE and the DFPT family carries KPAR.  A
    single ``incar_overrides`` applied to both would put KPAR next to NCORE in
    one INCAR, which VASP rejects.
    """

    from htvasp.dfpt_qha.maker import DfptQhaMaker

    maker = DfptQhaMaker(recipe=recipe(), relax_incar_overrides={"NCORE": 4},
                         dfpt_incar_overrides={"KPAR": 4})
    for child in (maker.initial_relax_maker, maker.eos_relax_maker):
        generator = child.input_set_generator
        assert generator.incar_overrides == {"NCORE": 4}
        generator.structure = structure()
        incar = generator.incar
        assert incar["NCORE"] == 4
        assert "KPAR" not in incar


def test_the_worker_takes_ncore_and_kpar_as_separate_knobs():
    from htvasp.workflows.dfpt_qha import DfptQhaWorker

    worker = DfptQhaWorker(recipe=recipe(), ncore=4, kpar=4)
    maker = worker.flow_maker
    for child in (maker.initial_relax_maker, maker.eos_relax_maker):
        assert child.input_set_generator.incar_overrides == {"NCORE": 4}
    assert maker.dfpt_incar_overrides == {"KPAR": 4}
    assert maker.relax_incar_overrides == {"NCORE": 4}

    with pytest.raises(ValueError, match="ncore must be positive"):
        DfptQhaWorker(recipe=recipe(), ncore=0)
    with pytest.raises(ValueError, match="kpar must be positive"):
        DfptQhaWorker(recipe=recipe(), kpar=-2)


def test_the_default_incar_still_carries_exactly_one_of_the_two_keys():
    """The defaults must not regress into an INCAR VASP would reject."""

    from htvasp.dfpt_qha.input_set import DfptInputSetGenerator, DfptRelaxSetGenerator

    relax = DfptRelaxSetGenerator(recipe=recipe(), isif=2)
    relax.structure = structure()
    assert relax.incar["NCORE"] == 8
    assert "KPAR" not in relax.incar

    dfpt = DfptInputSetGenerator(recipe=recipe())
    dfpt.structure = structure()
    assert dfpt.incar["KPAR"] == 2
    assert "NCORE" not in dfpt.incar


def test_the_maker_forwards_incar_overrides_to_every_child():
    from htvasp.dfpt_qha.maker import DfptQhaMaker

    maker = DfptQhaMaker(recipe=recipe(), incar_overrides={"ENCUT": 700})
    for child in (maker.initial_relax_maker, maker.eos_relax_maker):
        assert child.input_set_generator.incar_overrides == {"ENCUT": 700}


def test_job_names_come_from_one_place():
    """The builder and the Worker must agree on what a job is called."""

    from htvasp.dfpt_qha import names
    from htvasp.workflows.dfpt_qha import DfptQhaWorker

    flow = _spawned_flow()
    flow_names = [job.name for job in flow.jobs]
    assert names.STABILITY in flow_names
    assert names.ANALYSIS in flow_names
    assert any(name.startswith(names.VALIDATION) for name in flow_names)
    assert any(name.startswith(names.COLLECT) for name in flow_names)
    assert any(name.startswith(names.DFPT) for name in flow_names)

    worker = DfptQhaWorker(recipe=recipe())
    assert worker.get_result.__defaults__ == (names.ANALYSIS,)


def test_phonopy_construction_accepts_the_symmetry_parameters():
    """primitive_matrix decides the commensurate-q set; it must be adjustable."""

    from htvasp.dfpt_qha.phonons import phonopy_for

    cell = structure()
    default = phonopy_for(cell, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    explicit = phonopy_for(
        cell, [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
        primitive_matrix="auto", symprec=1e-5,
    )
    assert len(explicit.primitive) == len(default.primitive)
    # an explicit primitive matrix that keeps the conventional cell
    conventional = phonopy_for(
        cell, [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
        primitive_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    )
    assert len(conventional.primitive) == len(cell)


def test_a_missing_initial_moment_names_the_element():
    """An element outside the table must fail where the table is, not inside
    pymatgen's INCAR expansion (`KeyError: 0`)."""

    from htvasp.dfpt_qha.input_set import DfptInputSetGenerator

    generator = DfptInputSetGenerator(recipe=recipe(), magmom={"Ta": 3.0})
    generator.structure = structure()
    with pytest.raises(ValueError, match="no initial moment for \\['Al'\\]"):
        _ = generator.incar


# --------------------------------------------------------------------------
# antiferromagnetic endpoints
# --------------------------------------------------------------------------

def _afm_cell():
    from pymatgen.core import Lattice, Structure
    return Structure(Lattice.cubic(2.85), ["Cr", "Cr"],
                     [[0, 0, 0], [0.5, 0.5, 0.5]])


def test_expand_site_moments_puts_opposite_signs_on_the_two_sublattices():
    """The pattern is written on the 2-atom cell and used on the 16-atom supercell.

    A per-site list cannot be passed down unchanged: every atom has to be given
    the moment of the cell site it sits on.  If that expansion were wrong, the
    "antiferromagnetic" run would be ferromagnetic with half the moments
    flipped somewhere arbitrary -- and it would still converge.
    """
    import numpy as np

    from htvasp.dfpt_qha.input_set import expand_site_moments
    from htvasp.dfpt_qha.phonons import supercell_structure

    cell = _afm_cell()
    supercell = supercell_structure(cell, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    moments = expand_site_moments(supercell, cell, [3.0, -3.0])
    assert len(moments) == 16
    assert moments.count(3.0) == 8 and moments.count(-3.0) == 8

    # Every corner atom's eight nearest neighbours must be body atoms, so the
    # two signs must alternate across the bcc nearest-neighbour shell.
    frac = np.array([site.frac_coords for site in supercell]) * 2.0
    body = np.all(np.isclose(frac % 1.0, 0.5, atol=1e-3), axis=1)
    values = np.array(moments)
    assert set(values[~body]) == {3.0}
    assert set(values[body]) == {-3.0}
    distances = supercell.distance_matrix
    for index in np.where(~body)[0]:
        nearest = np.argsort(distances[index])[1:9]
        assert all(body[nearest])


def test_a_per_site_list_reaches_the_incar_at_every_stage():
    from htvasp.dfpt_qha.input_set import DfptInputSetGenerator, DfptRelaxSetGenerator
    from htvasp.dfpt_qha.phonons import supercell_structure

    cell = _afm_cell()
    generator = DfptRelaxSetGenerator(recipe=recipe(), isif=3, magmom=[3.0, -3.0],
                                      magmom_reference_cell=cell)
    generator.structure = cell
    assert generator.incar["MAGMOM"] == [3.0, -3.0]

    supercell = supercell_structure(cell, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    dfpt = DfptInputSetGenerator(recipe=recipe(), magmom=[3.0, -3.0],
                                 magmom_reference_cell=cell)
    dfpt.structure = supercell
    moments = dfpt.incar["MAGMOM"]
    assert len(moments) == 16
    assert moments.count(3.0) == 8 and moments.count(-3.0) == 8


def test_the_afm_maker_survives_jobflow_serialisation():
    """jobflow serialises makers through monty; a non-field attribute is lost."""
    from monty.json import MontyDecoder, jsanitize

    from htvasp.dfpt_qha.maker import DfptQhaMaker

    cell = _afm_cell()
    maker = DfptQhaMaker(recipe=recipe(), magmom=[3.0, -3.0], magmom_cell=cell)
    restored = MontyDecoder().process_decoded(
        jsanitize(maker, strict=True, allow_bson=True))
    assert restored.magmom == [3.0, -3.0]
    assert restored.magmom_cell is not None
    child = restored.initial_relax_maker.input_set_generator
    assert child.magmom_reference_cell is not None


def test_the_afm_endpoint_list_matches_the_spin_screen():
    """Only the endpoints whose screen chose afm may start staggered."""
    import main_dfpt_qha as portal

    assert portal.AFM_ENDPOINTS == {"BCC-Cr-Cr-2", "SER-Cr-1"}
    for name in portal.AFM_ENDPOINTS:
        structure, _ = portal.load_endpoint(name)
        moments, reference = portal.site_moments(name, structure)
        assert isinstance(moments, list) and reference is not None
        assert set(moments) == {3.0, -3.0}
    # a ferromagnetic Cr endpoint keeps the per-element table
    structure, _ = portal.load_endpoint("BCC-Cr-Hf-2")
    moments, reference = portal.site_moments("BCC-Cr-Hf-2", structure)
    assert isinstance(moments, dict) and reference is None
    assert moments["Cr"] == 3.0

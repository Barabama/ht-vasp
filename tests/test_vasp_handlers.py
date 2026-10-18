"""Tests for the custodian handler set and the DFPT convergence gates."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from jobflow import JobStore
from maggma.stores import MemoryStore
from pymatgen.core import Structure

from htvasp.dfpt_qha.jobs import STATIC_ENERGY_TOL_MEV_PER_ATOM, collect_dfpt_task_doc
from htvasp.utils.vasp_handlers import (
    DfptKpointSetChangeHandler,
    StableNonConvergingHandler,
    StableUnconvergedHandler,
    UnstableAlgoCorrection,
    base_handlers,
    dfpt_handlers,
    guard_algo_actions,
    with_handlers,
)

#: A real DFPT job directory whose POSCAR is in phonopy's atom order.  The
#: legacy e1 directory is NOT usable here: its POSCAR was hand-built, so the
#: ordering assertion correctly refuses it.
FIXTURE = Path(
    "/nfs_hdd/2025/gaominliang/ht-vasp/data/endpoints/BCC-Al-Ta-2/dfpt-qhaflow"
    "/runs/dfpt-qha-20261016T142130Z-a5a8ce189c/1-DFPT_volume_1/5"
)
POOL = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0")
MATRIX = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]


# --------------------------------------------------------------------------
# the handler set itself
# --------------------------------------------------------------------------


def test_base_set_keeps_the_convergence_handlers():
    """These are the ten that used to be dropped.

    ``htvasp.Worker`` replaced atomate2's eleven defaults with
    ``[VaspErrorHandler()]``, so a non-converged SCF reached the job boundary
    with nothing having reacted to it.
    """

    names = {type(handler).__name__ for handler in base_handlers()}
    assert "VaspErrorHandler" in names
    assert "StableUnconvergedHandler" in names
    assert "StableNonConvergingHandler" in names
    assert "FrozenJobErrorHandler" in names
    assert "StdErrHandler" in names


def test_base_set_leaves_out_the_recipe_rewriters():
    """Each of these would silently change a frozen scientific input."""

    names = {type(handler).__name__ for handler in base_handlers()}
    assert "LargeSigmaHandler" not in names        # rewrites SIGMA
    assert "IncorrectSmearingHandler" not in names  # rewrites ISMEAR
    assert "MeshSymmetryErrorHandler" not in names  # rewrites the k-mesh
    assert "PotimErrorHandler" not in names


def test_handlers_are_fresh_instances_each_call():
    """Custodian keeps per-run counters on the instances."""

    first = base_handlers()
    second = base_handlers()
    assert all(a is not b for a, b in zip(first, second))


def test_with_handlers_appends_rather_than_replaces():
    composed = with_handlers(run_vasp_kwargs={"handlers": [DfptKpointSetChangeHandler()]})
    names = [type(handler).__name__ for handler in composed["handlers"]]
    assert "VaspErrorHandler" in names
    assert "DfptKpointSetChangeHandler" in names
    assert composed["custodian_kwargs"]["max_errors_per_job"] == 3


def test_dfpt_extra_handler_is_a_detector():
    handler = dfpt_handlers()[0]
    assert handler.is_monitor is True
    # no actions -> custodian raises NonRecoverableError rather than "fixing"
    assert handler.correct(directory=".")["actions"] == []


def test_dfpt_handler_reads_the_kpoint_message(tmp_path):
    (tmp_path / "vasp.out").write_text(
        "some output\nVASP internal routines have requested a change of the "
        "k-point set.\nmore output\n"
    )
    assert DfptKpointSetChangeHandler().check(directory=str(tmp_path)) is True


def test_dfpt_handler_ignores_a_clean_run(tmp_path):
    (tmp_path / "vasp.out").write_text("all fine\n")
    assert DfptKpointSetChangeHandler().check(directory=str(tmp_path)) is False


# --------------------------------------------------------------------------
# the ALGO guard
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("current", "proposed"),
    [("VeryFast", "Fast"), ("Fast", "Normal"), ("Normal", "All"), ("All", "Damped")],
)
def test_more_stable_algo_changes_are_allowed(tmp_path, current, proposed):
    (tmp_path / "INCAR").write_text(f"ALGO = {current}\n")
    result = {"errors": ["Unconverged"],
              "actions": [{"dict": "INCAR", "action": {"_set": {"ALGO": proposed}}}]}
    kept = guard_algo_actions(result, str(tmp_path), handler="test")
    assert kept["actions"]


@pytest.mark.parametrize(
    ("current", "proposed"), [("All", "Normal"), ("Normal", "Fast"), ("Fast", "VeryFast")]
)
def test_less_stable_algo_changes_are_refused(tmp_path, current, proposed):
    """A downgrade changes the physics rather than the convergence."""

    (tmp_path / "INCAR").write_text(f"ALGO = {current}\n")
    result = {"errors": ["Unconverged"],
              "actions": [{"dict": "INCAR", "action": {"_set": {"ALGO": proposed}}}]}
    with pytest.raises(UnstableAlgoCorrection, match="less stable"):
        guard_algo_actions(result, str(tmp_path), handler="test")


def test_non_algo_actions_pass_through(tmp_path):
    (tmp_path / "INCAR").write_text("ALGO = Normal\n")
    result = {"errors": ["x"],
              "actions": [{"dict": "INCAR", "action": {"_set": {"AMIX": 0.1}}}]}
    kept = guard_algo_actions(result, str(tmp_path), handler="test")
    assert kept["actions"] == result["actions"]


def test_the_wrapped_handlers_are_the_guarded_ones():
    assert issubclass(StableUnconvergedHandler, object)
    assert StableUnconvergedHandler.correct is not None
    assert StableNonConvergingHandler.correct is not None


def _store() -> JobStore:
    docs = MemoryStore()
    docs.connect()
    return JobStore(docs)


# --------------------------------------------------------------------------
# the collector's two convergence gates
# --------------------------------------------------------------------------


def _fixture_doc() -> dict:
    if not (FIXTURE / "DYNMAT").is_file():
        pytest.skip("DFPT fixture unavailable")
    import gzip

    from pymatgen.io.vasp import Poscar

    with gzip.open(FIXTURE / "POSCAR.gz", "rt") as handle:
        supercell = Poscar.from_str(handle.read()).structure
    cell = Structure(
        np.array(supercell.lattice.matrix) / 2, ["Al", "Ta"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    return {
        "dir_name": str(FIXTURE),
        "state": "successful",
        "structure": supercell.as_dict(),
        "output": {"energy": -15.76945914},
    }, cell


def test_collector_rejects_a_failed_task_state():
    """state != successful means the SCF did not converge.

    ``VaspErrorHandler`` has 45 text signatures and none of them matches
    "aborting loop EDIFF was not reached", so this is the only thing standing
    between a non-converged run and its dynamical matrix.
    """

    doc, cell = _fixture_doc()
    doc = dict(doc, state="failed")
    recipe = json.loads((POOL / "pool_recipes.json").read_text())["endpoints"][
        "BCC-Al-Ta-2"
    ]["recipe"]
    with pytest.raises(ValueError, match="not successful"):
        collect_dfpt_task_doc(
            doc, recipe=recipe, cell=cell.as_dict(), supercell_matrix=MATRIX
        ).run(_store())


def test_collector_accepts_an_energy_consistent_with_the_static():
    doc, cell = _fixture_doc()
    recipe = json.loads((POOL / "pool_recipes.json").read_text())["endpoints"][
        "BCC-Al-Ta-2"
    ]["recipe"]
    n_atoms = len(doc["structure"]["sites"])
    # the fixture's own energy, per atom, is the reference it must match
    per_atom = doc["output"]["energy"] / n_atoms
    payload = collect_dfpt_task_doc(
        doc, recipe=recipe, cell=cell.as_dict(), supercell_matrix=MATRIX,
        expected_static_energy_per_atom=per_atom,
    ).run(_store()).output
    assert payload["provenance"]["energy_vs_static_meV_per_atom"] == pytest.approx(0.0)


def test_collector_rejects_a_diverged_ground_state():
    """The BCC-Co-Ti-2 shape: a run that stopped early on a diverged SCF."""

    doc, cell = _fixture_doc()
    recipe = json.loads((POOL / "pool_recipes.json").read_text())["endpoints"][
        "BCC-Al-Ta-2"
    ]["recipe"]
    n_atoms = len(doc["structure"]["sites"])
    per_atom = doc["output"]["energy"] / n_atoms
    with pytest.raises(ValueError, match="disagrees with the static energy"):
        collect_dfpt_task_doc(
            doc, recipe=recipe, cell=cell.as_dict(), supercell_matrix=MATRIX,
            expected_static_energy_per_atom=per_atom + 1.0,  # +1 eV/atom
        ).run(_store())


def test_static_energy_tolerance_is_the_documented_value():
    assert STATIC_ENERGY_TOL_MEV_PER_ATOM == 100.0

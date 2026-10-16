"""No-VASP tests for the DFPT TaskDoc collector."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from jobflow import JobStore
from maggma.stores import MemoryStore
from pymatgen.core import Structure

from htvasp.dfpt_qha.jobs import collect_dfpt_task_doc, validate_dfpt_task_doc
from htvasp.dfpt_qha.task_doc import DfptTaskDoc

#: A real DFPT run: 2x2x2 supercell of B2 AlTa, IBRION=8, VASP wrote
#: FORCE_CONSTANTS in phonopy format.  Used as the fixture precisely because it
#: is real -- a hand-built one hid two fatal path bugs the first time round.
DFPT_DIR = Path("/nfs_ssd/tmp/legacy-qha-pool/e1-dfpt-al-ta")
RECIPE_FILE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")
ENDPOINT = "BCC-Al-Ta-2"
HAVE_FIXTURE = (DFPT_DIR / "FORCE_CONSTANTS").is_file()
needs_fixture = pytest.mark.skipif(not HAVE_FIXTURE, reason="DFPT fixture unavailable")


def _store() -> JobStore:
    docs = MemoryStore()
    docs.connect()
    return JobStore(docs)


def recipe_dict() -> dict:
    return json.loads(RECIPE_FILE.read_text())["endpoints"][ENDPOINT]["recipe"]


def _task_doc() -> dict:
    """A minimal TaskDoc stand-in carrying the run's directory and structure."""

    return {
        "dir_name": str(DFPT_DIR),
        "state": "successful",
        "structure": Structure.from_file(str(DFPT_DIR / "CONTCAR")).as_dict(),
    }


@needs_fixture
def test_collector_reads_force_constants_into_the_contract():
    payload = collect_dfpt_task_doc(_task_doc(), recipe=recipe_dict()).run(_store()).output
    parsed = DfptTaskDoc.model_validate(payload)
    provenance = parsed.force_constants_provenance
    assert provenance.source == "FORCE_CONSTANTS"
    assert provenance.parser == "phonopy.file_IO.parse_FORCE_CONSTANTS"
    assert provenance.n_atoms == 16
    assert provenance.shape == (16, 16, 3, 3)
    assert provenance.ordering_fingerprint
    assert len(parsed.force_constants) == 16


@needs_fixture
def test_collected_matrix_satisfies_the_acoustic_sum_rule():
    """The one cheap check that a matrix is a force-constant matrix at all.

    The vasprun <dynmat> block also has shape (16, 16, 3, 3) and is symmetric,
    but it violates this rule by 65% of its own scale, so using it would give
    plausible-looking frequencies that are simply wrong.
    """

    payload = collect_dfpt_task_doc(_task_doc(), recipe=recipe_dict()).run(_store()).output
    parsed = DfptTaskDoc.model_validate(payload)
    residual = parsed.force_constants_provenance.max_acoustic_sum_residual
    assert residual is not None
    scale = max(abs(x) for row in parsed.force_constants for m in row for v in m for x in v)
    assert residual / scale < 1e-3, (residual, scale)


def test_collector_requires_the_force_constants_file(tmp_path):
    with pytest.raises(ValueError, match="FORCE_CONSTANTS"):
        collect_dfpt_task_doc(
            {"dir_name": str(tmp_path), "state": "successful", "structure": {}},
            recipe=recipe_dict(),
        ).run(_store())


def test_collector_strips_the_hostname_prefix(tmp_path):
    """TaskDoc.dir_name is an emmet URI: "host:/abs/path".

    Path("host:/abs/path") is relative and never exists, so every volume's
    collector used to fail after VASP had already consumed the allocation.
    """

    if not HAVE_FIXTURE:
        pytest.skip("DFPT fixture unavailable")
    prefixed = f"{__import__('socket').gethostname()}:{DFPT_DIR}"
    payload = collect_dfpt_task_doc(
        {"dir_name": prefixed, "state": "successful", "structure": _task_doc()["structure"]},
        recipe=recipe_dict(),
    ).run(_store()).output
    assert DfptTaskDoc.model_validate(payload).force_constants_provenance.n_atoms == 16


def test_validation_stops_on_a_malformed_document():
    result = validate_dfpt_task_doc({"bad": True}, recipe=recipe_dict()).run(_store())
    assert result.output["ok"] is False
    assert result.stop_children is True


@needs_fixture
def test_validation_rejects_a_matrix_from_the_wrong_cell():
    """A Gamma-only result from the 2-atom cell is square and non-empty.

    Without the atom-count check against the recipe it would be recorded as a
    valid volume, and the finite-q instability would be lost silently.
    """

    payload = collect_dfpt_task_doc(_task_doc(), recipe=recipe_dict()).run(_store()).output
    # A valid recipe describing a DIFFERENT cell: 4-atom cell, factor 2, so the
    # phonon supercell is 32 atoms rather than the 16 the matrix holds.
    wrong_recipe = dict(recipe_dict())
    wrong_recipe["kpoints"] = dict(
        wrong_recipe["kpoints"],
        grid=[8, 8, 8], supercell_grid=[4, 4, 4], cell_atoms=4,
        n3V=8**3 * float(Structure.from_dict(_task_doc()["structure"]).volume) / 8,
    )
    result = validate_dfpt_task_doc(payload, recipe=wrong_recipe).run(_store())
    assert result.output["ok"] is False
    assert result.stop_children is True
    assert "phonon supercell" in result.output["error"]


@needs_fixture
def test_validation_accepts_the_real_collected_document():
    collected = collect_dfpt_task_doc(_task_doc(), recipe=recipe_dict()).run(_store()).output
    result = validate_dfpt_task_doc(collected, recipe=recipe_dict()).run(_store())
    assert result.output["ok"] is True
    assert result.stop_children is False
    # the validated document travels with the verdict so the thermal job can
    # consume it; a gate that returns only {ok: bool} strands the data
    assert result.output["dfpt_doc"]["force_constants"]

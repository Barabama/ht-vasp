"""No-VASP tests for the DFPT TaskDoc collector."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from jobflow import Flow, JobStore, job
from maggma.stores import MemoryStore
from htvasp.dfpt_qha.jobs import collect_dfpt_task_doc, validate_dfpt_task_doc
from htvasp.dfpt_qha.task_doc import DfptTaskDoc

DFPT_DIR = Path("/nfs_ssd/tmp/legacy-qha-pool/e1-dfpt-al-ta")
RECIPE = json.loads(
    Path("/nfs_ssd/tmp/legacy-qha-pool/e1-dfpt-al-ta/recipe.json").read_text()
)


def _store() -> JobStore:
    docs = MemoryStore()
    docs.connect()
    return JobStore(docs)


@pytest.mark.skipif(not (DFPT_DIR / "vasprun.xml").is_file(), reason="DFPT fixture unavailable")
def test_collector_reads_vasprun_dynmat_into_contract():
    task_doc = {"dir_name": str(DFPT_DIR), "state": "successful"}
    job = collect_dfpt_task_doc(task_doc, recipe=RECIPE, vasprun_path=DFPT_DIR)
    result = job.run(_store()) if hasattr(job, "run") else None
    assert result is not None
    payload = result.output
    parsed = DfptTaskDoc.model_validate(payload)
    assert parsed.force_constants_provenance.source == "vasprun.xml"
    assert parsed.force_constants_provenance.parser == "pymatgen.Vasprun.force_constants"
    assert parsed.force_constants_provenance.n_atoms == 16
    assert parsed.force_constants_provenance.shape == (16, 16, 3, 3)
    assert len(parsed.force_constants) == 16


def test_collector_requires_vasprun(tmp_path):
    with pytest.raises(ValueError, match="no vasprun.xml"):
        collect_dfpt_task_doc(
            {"dir_name": str(tmp_path), "state": "successful"}, recipe=RECIPE
        ).run(_store())


def test_validation_job_stops_on_invalid_collector_output():
    result = validate_dfpt_task_doc({"bad": True}).run(_store())
    assert result.output["ok"] is False
    assert result.stop_children is True


def test_validation_job_accepts_collector_output():
    if not (DFPT_DIR / "vasprun.xml").is_file():
        pytest.skip("DFPT fixture unavailable")
    collected = collect_dfpt_task_doc(
        {"dir_name": str(DFPT_DIR), "state": "successful"},
        recipe=RECIPE,
        vasprun_path=DFPT_DIR,
    ).run(_store()).output
    result = validate_dfpt_task_doc(collected).run(_store())
    assert result.output["ok"] is True
    assert result.stop_children is False

"""Offline tests for the SQS calculation Flow and observation job."""

from __future__ import annotations

import hashlib
import json

import pytest
from jobflow.core.reference import OutputReference
from pymatgen.core import Lattice, Structure

from htvasp.sqs.jobs import make_static_observation
from htvasp.sqs.maker import SqsCalculationMaker
from htvasp.sqs.task_doc import SqsStaticObservation
from htvasp.workflows import SqsMixingWorker


def _structure() -> Structure:
    return Structure(Lattice.cubic(2.86), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])


def _fake_static(state: str = "successful") -> dict:
    return {
        "state": state,
        "energy": -20.0,
        "structure": _structure().as_dict(),
        "dir_name": "/tmp/static",
    }


def _contains_reference(value, uuid: str, attribute: str) -> bool:
    if isinstance(value, OutputReference):
        expected = () if attribute == "output" else (("a", attribute),)
        return value.uuid == uuid and value.attributes == expected
    if isinstance(value, dict):
        return any(_contains_reference(item, uuid, attribute) for item in value.values())
    if isinstance(value, (list, tuple)):
        return any(_contains_reference(item, uuid, attribute) for item in value)
    return False


def test_worker_builds_relax_static_observation_flow():
    worker = SqsMixingWorker(vasp_args={"vasp_cmd": "echo test"}, calculation_id="BCC-Al-Nb-x0.25")
    flow = worker._make_flow(_structure())
    assert [job.name for job in flow.jobs] == ["sqs relax", "sqs static", "sqs observation"]
    assert flow.output.uuid == flow.jobs[2].uuid
    assert _contains_reference(flow.jobs[1].function_args, flow.jobs[0].uuid, "structure")
    assert _contains_reference(flow.jobs[1].function_kwargs, flow.jobs[0].uuid, "dir_name")
    assert _contains_reference(flow.jobs[2].function_args, flow.jobs[1].uuid, "output")


def test_relax_maker_uses_full_cell_relaxation():
    worker = SqsMixingWorker(
        vasp_args={"vasp_cmd": "echo test"},
        relax_incar={"ISIF": 2},
    )
    assert worker.flow_maker.relax_maker.input_set_generator.user_incar_settings["ISIF"] == 3
    assert type(worker.flow_maker.relax_maker).__name__ == "RelaxMaker"


def test_static_maker_uses_clean_scf_after_relax():
    static_maker = SqsMixingWorker(vasp_args={"vasp_cmd": "echo test"}).flow_maker.static_maker
    settings = static_maker.input_set_generator.user_incar_settings
    assert settings["ISTART"] == 0
    assert settings["ICHARG"] == 2
    assert settings["ALGO"] == "Normal"
    assert settings["NELM"] == 200
    assert static_maker.copy_vasp_kwargs == {}


def test_observation_standardizes_static_task_document():
    task = _fake_static()
    response = make_static_observation(task, calculation_id="BCC-Al-Nb-x0.25").function(
        task, "BCC-Al-Nb-x0.25"
    )
    observation = SqsStaticObservation.model_validate(response)
    assert observation.energy_total_eV == -20.0
    assert observation.energy_per_atom_eV == -10.0
    assert observation.n_sites == 2
    assert len(observation.task_document_sha256) == 64
    canonical = json.dumps(
        {
            "state": "successful",
            "energy": -20.0,
            "structure": _structure().as_dict(),
            "dir_name": "/tmp/static",
        },
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    assert observation.task_document_sha256 == hashlib.sha256(canonical).hexdigest()


def test_observation_accepts_nested_model_dump_shape():
    task = {"output": _fake_static()}
    response = make_static_observation(task, calculation_id="ref-Fe").function(task, "ref-Fe")
    assert response["calculation_id"] == "ref-Fe"


@pytest.mark.parametrize(
    "task, message",
    [
        ({**_fake_static("failed")}, "did not succeed"),
        ({"state": "successful", "energy": -1.0}, "no final structure"),
        ({"state": "successful", "structure": _structure().as_dict()}, "no total"),
    ],
)
def test_observation_rejects_invalid_static_task(task, message):
    job = make_static_observation(task, calculation_id="bad")
    with pytest.raises(ValueError, match=message):
        job.function(task, "bad")


def test_worker_uses_observation_as_default_result_name():
    worker = SqsMixingWorker(vasp_args={"vasp_cmd": "echo test"})
    assert worker.get_result.__defaults__ == ("sqs observation",)


def test_maker_accepts_calculation_spec_metadata():
    maker = SqsCalculationMaker(run_vasp_kwargs={"vasp_cmd": "echo test"})
    flow = maker.make(_structure(), {"calculation_id": "endpoint-Fe"})
    assert flow.jobs[-1].function_kwargs["calculation_id"] == "endpoint-Fe"

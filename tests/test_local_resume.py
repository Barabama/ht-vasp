"""Dynamic sub-flows must be resumable.

A DFPT endpoint expands into a ``Response(replace=Flow(...))`` sub-flow of one
VASP job per volume.  On resume the parent job is re-executed (its stored
output is None when the sub-flow did not finish), which rebuilds the sub-flow
with brand-new job uuids.  Forcing ``resume=False`` on that sub-flow re-ran
every volume whose TaskDoc was already in the store -- hours of VASP per
interrupted endpoint.
"""

from __future__ import annotations

from jobflow import Flow, JobStore, Response, job
from maggma.stores import MemoryStore

from htvasp.utils.local import run_locally_custom

EXECUTIONS = {"first": 0, "second": 0, "fail": True}


@job
def first():
    EXECUTIONS["first"] += 1
    return {"value": 1}


@job
def second(previous):
    EXECUTIONS["second"] += 1
    if EXECUTIONS["fail"]:
        raise RuntimeError("simulated interruption")
    return {"value": previous["value"] + 1}


@job
def parent():
    a = first()
    b = second(a.output)
    return Response(replace=Flow([a, b], output=b.output, name="sub"))


def _store() -> JobStore:
    store = MemoryStore()
    store.connect()
    return JobStore(store)


def test_dynamic_sub_flow_reuses_completed_jobs_on_resume(tmp_path):
    EXECUTIONS.update({"first": 0, "second": 0, "fail": True})
    store = _store()

    run_locally_custom(parent(), store=store, root_dir=tmp_path, resume=True)
    assert EXECUTIONS["first"] == 1
    assert EXECUTIONS["second"] == 1

    # the interruption is over; resuming must not re-run the completed job
    EXECUTIONS["fail"] = False
    run_locally_custom(parent(), store=store, root_dir=tmp_path, resume=True)
    assert EXECUTIONS["first"] == 1, "completed sub-flow job was executed again"
    assert EXECUTIONS["second"] == 2

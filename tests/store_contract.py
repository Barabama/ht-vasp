"""The store contract: a finished store must be readable by a fresh process.

Every jobflow flow this project runs is *supposed* to be resumable -- ``Worker``
defaults to ``resume=True`` and ``run_locally_custom`` matches stored jobs by
``(name, index)``.  That promise has one precondition nothing was checking: the
store directory must be self-contained.  Anything a document points at has to
live on disk beside it.

It did not.  ``Worker.run_flow`` created the "data" additional store as a
``MemoryStore``, and atomate2's ``@vasp_job`` passes ``data=_DATA_OBJECTS``, so
every VASP TaskDoc (force constants, DOS, band structure ...) was replaced by a
``blob_uuid`` pointer whose target existed only in the writer's RAM.  A resume
is by definition a new process, so the pointer dangled, ``get_output(load=True)``
raised ``KeyError``, ``_validate_job_output`` read that as "unfinished" and
re-ran the job.  FCC-Cr-Nb-4 recomputed six DFPT volumes that way on 2026-10-29.

Two things are deliberately *not* treated as violations, because an interrupted
flow leaves them behind by design:

* a row with ``output: null`` -- a job that never finished;
* a row whose output holds an ``OutputReference`` to such a job -- jobflow
  cannot resolve a reference to something that was never written, and that is
  exactly what "the run stopped here" looks like.

The defect is narrower and unambiguous: a ``blob_uuid`` pointer that has no
target on disk.  That row *did* finish; its payload is simply gone.

``assert_store_resumable`` is the assertion that would have caught the bug.  It
is deliberately built on the documented layout rather than on the production
store factory: if that layout changes, this must fail and force a deliberate
update instead of silently following along.
"""

from __future__ import annotations

import json
from pathlib import Path

DOCS_NAME = "store.json"
#: Where offloaded blobs are expected to live.  Hard-coded on purpose -- see the
#: module docstring.  A store that keeps them anywhere else is not resumable by
#: the contract this project relies on.
BLOBS_NAME = "store-data.json"


def dangling_blobs(store_dir: Path) -> list[tuple[str, str]]:
    """``(job name, blob_uuid)`` for every pointer whose target is not on disk.

    Walks the raw JSON rather than going through ``get_output``: the point is to
    separate "this pointer has no target" from "this reference points at a job
    that never ran", and only the first one is a defect.
    """
    store_dir = Path(store_dir)
    documents = json.loads((store_dir / DOCS_NAME).read_text())
    available: set[str] = set()
    blobs_path = store_dir / BLOBS_NAME
    if blobs_path.is_file():
        payload = json.loads(blobs_path.read_text() or "[]")
        if isinstance(payload, list):
            available = {
                entry.get("blob_uuid") for entry in payload if entry.get("blob_uuid")
            }

    dangling: list[tuple[str, str]] = []
    for document in documents:
        for blob_uuid in _blob_uuids(document.get("output")):
            if blob_uuid not in available:
                dangling.append((document.get("name") or "?", blob_uuid))
    return dangling


def _blob_uuids(node) -> list[str]:
    """Every ``blob_uuid`` recorded anywhere under ``node``."""
    found: list[str] = []
    if isinstance(node, dict):
        if isinstance(node.get("blob_uuid"), str):
            found.append(node["blob_uuid"])
        for value in node.values():
            found.extend(_blob_uuids(value))
    elif isinstance(node, list):
        for value in node:
            found.extend(_blob_uuids(value))
    return found


def assert_store_resumable(store_dir: Path, *, allow: tuple[str, ...] = ()) -> None:
    """Fail unless every blob pointer in the store resolves from disk alone.

    Args:
        store_dir: The directory a finished (or interrupted) flow left behind.
        allow: Job names whose dangling pointers are expected.  Empty by
            default: a partial *run* is normal, but a finished job whose payload
            evaporated is not.
    """
    store_dir = Path(store_dir)
    assert (store_dir / DOCS_NAME).is_file(), f"no {DOCS_NAME} in {store_dir}"
    failures = [row for row in dangling_blobs(store_dir) if row[0] not in allow]
    if not failures:
        return
    detail = "\n".join(f"  {name}: blob {uuid}" for name, uuid in failures)
    hint = ""
    if not (store_dir / BLOBS_NAME).is_file():
        hint = (
            f"\n{BLOBS_NAME} was never written.  That is what a process-local "
            "additional store produces: the pointers survive in store.json, "
            "their targets die with the process that wrote them."
        )
    raise AssertionError(
        f"{len(failures)} blob pointer(s) do not resolve, so a resume would "
        f"silently re-run those jobs:\n{detail}{hint}"
    )

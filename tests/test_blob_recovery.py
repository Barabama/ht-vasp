"""Rebuilding a lost blob store from what ``store.json`` still holds.

``Worker.run_flow`` used to create the "data" additional store as a
``MemoryStore``, so atomate2's offloaded VASP payloads (``force_constants``,
``normalmode_eigenvecs``) died with the process that wrote them.  The pointers
survived; on the next process they dangled, ``_validate_job_output`` read that
as "unfinished", and the volumes were re-run.  FCC-Cr-Nb-4 lost six that way.

The payloads are not actually gone.  The collector embeds the whole resolved
TaskDoc in its own output, and the VASP directory is still archived, so a
dangling pointer can be repaired from either.  ``scripts/recover_blob_store.py``
does that and demands the two sources agree before writing anything --
otherwise it is a guess, not a reconstruction.

These tests use the vendored DFPT job directory, so the archive check runs on
a real ``vasprun.xml`` rather than on a mock.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from htvasp.utils.local import _validate_job_output

FIXTURE = Path(__file__).parent / "fixtures" / "BCC-Al-Ta-2-dfpt-volume-1"

#: Where the collector keeps the resolved TaskDoc inside its output, and the
#: two fields atomate2 offloads.  Both are chosen by atomate2, not here.
TASK_DOC = "task_doc"
CALC_PATH = ("calcs_reversed", 0, "output")
OFFLOADED = ("force_constants", "normalmode_eigenvecs")


def _payloads() -> dict[str, list]:
    from pymatgen.io.vasp.outputs import Vasprun

    parsed = Vasprun(
        FIXTURE / "vasprun.xml.gz",
        parse_potcar_file=False,
        parse_dos=False,
        parse_eigen=False,
    )
    return {field: np.asarray(getattr(parsed, field)).tolist() for field in OFFLOADED}


def _store_row(name: str, uuid: str, output) -> dict:
    return {"uuid": uuid, "index": 1, "name": name, "output": output, "completed_at": None}


def _broken_run(tmp_path: Path, *, with_archive: bool = True) -> tuple[Path, dict]:
    """A run directory shaped like the one FCC-Cr-Nb-4 was left with.

    The DFPT volume job's output holds ``blob_uuid`` pointers where its payload
    used to be (that is literally what jobflow writes); the collector's output
    holds the same values inline, because it received the resolved TaskDoc.
    """
    payloads = _payloads()
    run = tmp_path / "run"

    pointers = {}
    for field in OFFLOADED:
        pointers[field] = {
            "@class": "",
            "@module": "",
            "blob_uuid": f"blob-{field}",
            "store": "data",
        }

    calc = {"output": dict(pointers)}
    volume_output = {"calcs_reversed": [calc]}
    collect_doc = {"calcs_reversed": [{"output": dict(payloads)}]}

    rows = [
        _store_row("DFPT volume 1/1", "uuid-volume", volume_output),
        _store_row("collect DFPT volume 1/1", "uuid-collect",
                   {TASK_DOC: collect_doc}),
    ]
    run.mkdir()
    (run / "store.json").write_text(json.dumps(rows))

    if with_archive:
        volume_dir = run / "1-DFPT_volume_1" / "7"
        volume_dir.mkdir(parents=True)
        shutil.copy(FIXTURE / "vasprun.xml.gz", volume_dir / "vasprun.xml.gz")
    return run, payloads


def _run_recovery(run: Path, *extra: str) -> tuple[int, str]:
    """Call the script's ``main`` in-process, returning its exit code and output."""
    import importlib.util

    script = Path(__file__).resolve().parent.parent / "scripts" / "recover_blob_store.py"
    spec = importlib.util.spec_from_file_location("recover_blob_store", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    import contextlib
    import io

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        code = module.main([str(run), *extra])
    return code, buffer.getvalue()


def test_the_recovered_store_resolves_every_pointer(tmp_path):
    """The point of the exercise: a resume can now load the offloaded output."""
    run, payloads = _broken_run(tmp_path)
    code, output = _run_recovery(run, "--write")
    assert code == 0, output

    from jobflow import JobStore
    from maggma.stores import JSONStore

    store = JobStore(
        JSONStore(str(run / "store.json"), read_only=True),
        additional_stores={
            "data": JSONStore(str(run / "store-data.json"), read_only=True)
        },
    )
    store.connect()
    loaded = store.get_output(uuid="uuid-volume", which=1, load=True)

    for field in OFFLOADED:
        recovered = loaded["calcs_reversed"][0]["output"][field]
        assert np.array_equal(np.asarray(recovered), np.asarray(payloads[field])), (
            f"{field} was not restored to the value the pointer referred to"
        )


def test_a_missing_recovery_is_a_dangling_pointer(tmp_path):
    """Without the blob store, the volume reads as unfinished -- the original bug."""
    run, _payloads_ = _broken_run(tmp_path, with_archive=False)
    from jobflow import JobStore
    from maggma.stores import JSONStore

    (run / "store-data.json").write_text("[]")
    store = JobStore(
        JSONStore(str(run / "store.json"), read_only=False),
        additional_stores={
            "data": JSONStore(str(run / "store-data.json"), read_only=False)
        },
    )
    store.connect()

    class _Job:
        uuid = "uuid-volume"
        index = 1
        name = "DFPT volume 1/1"

    assert _validate_job_output(_Job, store) is False


def test_the_archive_check_refuses_a_payload_that_does_not_match(tmp_path):
    """A reconstruction that disagrees with vasprun.xml is a guess, not a recovery."""
    run, _ = _broken_run(tmp_path)
    rows = json.loads((run / "store.json").read_text())
    inline = rows[1]["output"][TASK_DOC]["calcs_reversed"][0]["output"]["force_constants"]
    inline[0][0][0][0] += 1.0
    (run / "store.json").write_text(json.dumps(rows))

    code, output = _run_recovery(run, "--write")
    assert code == 1
    assert "differs from" in output
    assert not (run / "store-data.json").exists(), (
        "an unverified reconstruction must not be written"
    )


def test_a_missing_archive_is_reported_rather_than_skipped(tmp_path):
    """No vasprun means the payload cannot be independently checked."""
    run, _ = _broken_run(tmp_path, with_archive=False)
    code, output = _run_recovery(run, "--write")
    assert code == 1
    assert "no archived vasprun.xml" in output
    assert not (run / "store-data.json").exists()


def test_an_existing_blob_store_is_never_overwritten(tmp_path):
    """The script repairs a store; it does not replace one that already works."""
    run, _ = _broken_run(tmp_path)
    (run / "store-data.json").write_text("[]")
    code, output = _run_recovery(run, "--write")
    assert code == 0
    assert (run / "store-data.json").read_text() == "[]"
    assert "already exists" in output


def test_a_pointer_without_a_donor_row_is_refused(tmp_path):
    """A job whose payload nothing else recorded must fail loudly, not silently."""
    run, _ = _broken_run(tmp_path)
    rows = json.loads((run / "store.json").read_text())
    rows[1]["name"] = "collect DFPT volume 2/1"      # no longer the donor
    (run / "store.json").write_text(json.dumps(rows))

    code, output = _run_recovery(run, "--write")
    assert code == 1
    assert "no donor row" in output
    assert not (run / "store-data.json").exists()


def test_the_sidecar_records_where_each_payload_came_from(tmp_path):
    """A repaired store has to say so: the provenance is part of the evidence."""
    run, _ = _broken_run(tmp_path)
    code, output = _run_recovery(run, "--write")
    assert code == 0, output

    sidecar = json.loads((run / "store-data.recovery.json").read_text())
    assert sidecar["schema"] == "htvasp-blob-recovery/v1"
    assert sidecar["verified_against"] == "archived vasprun.xml"
    assert {item["job"] for item in sidecar["blobs"]} == {"collect DFPT volume 1/1"}

    blobs = json.loads((run / "store-data.json").read_text())
    assert {blob["blob_uuid"] for blob in blobs} == {f"blob-{f}" for f in OFFLOADED}

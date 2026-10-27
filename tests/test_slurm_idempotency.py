"""Offline regression tests for durable Slurm intent submission.

These tests monkeypatch the subprocess boundary.  They never invoke sbatch,
squeue, sacct, VASP, or a cluster scheduler.
"""

from __future__ import annotations

import json
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from htvasp.slurm.manager import (
    SlurmJobManager,
    SlurmConfig,
    atomic_json_write,
    artifact_marker_matches,
    write_completion_marker,
)


class FakeCompleted:
    def __init__(self, stdout="", returncode=0):
        self.stdout = stdout
        self.stderr = ""
        self.returncode = returncode


class FakeSlurm:
    def __init__(self, *, state="PENDING", fail_after_submit=False):
        self.state = state
        self.fail_after_submit = fail_after_submit
        self.sbatch_calls = []
        self.jobs = {}
        self._next_id = 100
        self.lock = threading.Lock()

    def __call__(self, args, **kwargs):
        command = " ".join(args) if isinstance(args, list) else str(args)
        if command.startswith("sbatch"):
            with self.lock:
                self.sbatch_calls.append(command)
                job_id = str(self._next_id)
                self._next_id += 1
                marker = next(
                    (item.split("=", 1)[1] for item in args if item.startswith("--comment=")),
                    "",
                )
                self.jobs[job_id] = marker
            return FakeCompleted(f"Submitted batch job {job_id}\n")
        if isinstance(args, list) and args and args[0] == "squeue":
            if "--jobs" in args:
                return FakeCompleted(f"{self.state}\n" if self.state != "UNKNOWN" else "")
            lines = [
                f"{job_id}|{self.state}|{comment}"
                for job_id, comment in self.jobs.items()
                if self.state != "UNKNOWN"
            ]
            return FakeCompleted("\n".join(lines) + ("\n" if lines else ""))
        if isinstance(args, list) and args and args[0] == "sacct":
            if self.state == "UNKNOWN":
                return FakeCompleted("", returncode=1)
            if "--format=JobIDRaw,State,Comment" in args:
                lines = [
                    f"{job_id}|{self.state}|{comment}"
                    for job_id, comment in self.jobs.items()
                ]
                return FakeCompleted("\n".join(lines) + ("\n" if lines else ""))
            return FakeCompleted(f"{self.state}\n")
        raise AssertionError(f"unexpected scheduler command: {args}")


def _request(tmp_path: Path, *, stage="endpoints", run_id=None, dependency=()):
    config = SlurmConfig(job_name="test", output_log="test.log", dependency="")
    return dict(
        system_name="BCC-Al-Nb",
        stage=stage,
        recipe_version="test-v1",
        settings_fingerprint="settings-a",
        manifest_fingerprint="manifest-a",
        command="python worker.py --stage " + stage,
        config=config,
        ledger_path=tmp_path / "slurm_submissions.json",
        work_dir=tmp_path,
        run_id=run_id,
        dependency_job_ids=list(dependency),
    )


def test_identical_request_submits_exactly_once(monkeypatch, tmp_path):
    fake = FakeSlurm(state="PENDING")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    manager = SlurmJobManager()
    first = manager.submit_intent(**_request(tmp_path))
    second = manager.submit_intent(**_request(tmp_path))
    assert first.submitted is True
    assert second.existing is True
    assert second.job_id == first.job_id
    assert len(fake.sbatch_calls) == 1
    ledger = json.loads((tmp_path / "slurm_submissions.json").read_text())
    assert ledger["schema_version"] == 1
    assert len(ledger["entries"]) == 1
    assert ledger["entries"][0]["state"] == "active"


def test_concurrent_requests_submit_exactly_once(monkeypatch, tmp_path):
    fake = FakeSlurm(state="PENDING")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    request = _request(tmp_path)
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(lambda _: SlurmJobManager().submit_intent(**request), range(2)))
    assert len(fake.sbatch_calls) == 1
    assert {result.job_id for result in results} == {"100"}


def test_active_completed_failed_unknown_reconciliation(monkeypatch, tmp_path):
    fake = FakeSlurm(state="PENDING")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    manager = SlurmJobManager()
    request = _request(tmp_path)
    submitted = manager.submit_intent(**request)
    assert submitted.job_id == "100"

    fake.state = "COMPLETED"
    completed = manager.submit_intent(**request)
    assert completed.state == "completed"
    assert completed.job_id == "100"
    assert len(fake.sbatch_calls) == 1

    fake.state = "FAILED"
    failed = manager.submit_intent(**request)
    assert failed.state == "failed"
    assert len(fake.sbatch_calls) == 1

    unknown_request = _request(tmp_path, stage="points", run_id="unknown")
    fake.state = "UNKNOWN"
    unknown = manager.submit_intent(**unknown_request)
    assert unknown.state == "submitted"
    fake.state = "UNKNOWN"
    unknown_again = manager.submit_intent(**unknown_request)
    assert unknown_again.state == "reconcile_required"
    assert unknown_again.reconcile_required is True
    assert len(fake.sbatch_calls) == 2


def test_failed_job_requires_explicit_retry_and_keeps_retry_of(monkeypatch, tmp_path):
    fake = FakeSlurm(state="FAILED")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    manager = SlurmJobManager()
    request = _request(tmp_path)
    first = manager.submit_intent(**request)
    assert first.submitted is True
    no_retry = manager.submit_intent(**request)
    assert no_retry.state == "failed"
    retry = manager.submit_intent(**{**request, "retry": True})
    assert retry.submitted is True
    assert retry.job_id == "101"
    entries = SlurmJobManager.ledger_entries(request["ledger_path"])
    assert len(entries) == 2
    assert entries[1]["retry_of"] == entries[0]["intent_key"]
    assert entries[1]["attempt"] == 2


def test_successful_sbatch_with_ledger_failure_is_reconcile_required(monkeypatch, tmp_path):
    fake = FakeSlurm(state="PENDING")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    original = SlurmJobManager._write_ledger
    calls = {"count": 0}

    def fail_final_write(path, payload):
        calls["count"] += 1
        if calls["count"] == 3:
            raise OSError("simulated crash after sbatch")
        return original(path, payload)

    monkeypatch.setattr(SlurmJobManager, "_write_ledger", staticmethod(fail_final_write))
    result = SlurmJobManager().submit_intent(**_request(tmp_path))
    assert result.state == "reconcile_required"
    assert result.job_id == "100"
    assert len(fake.sbatch_calls) == 1


def test_crash_after_sbatch_is_recovered_by_intent_comment(monkeypatch, tmp_path):
    fake = FakeSlurm(state="PENDING")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    original = SlurmJobManager._write_ledger
    calls = {"count": 0}

    def fail_once(path, payload):
        calls["count"] += 1
        if calls["count"] == 3:
            raise OSError("simulated crash after sbatch")
        return original(path, payload)

    monkeypatch.setattr(SlurmJobManager, "_write_ledger", staticmethod(fail_once))
    request = _request(tmp_path)
    first = SlurmJobManager().submit_intent(**request)
    assert first.state == "reconcile_required"
    assert first.job_id == "100"

    monkeypatch.setattr(SlurmJobManager, "_write_ledger", staticmethod(original))
    recovered = SlurmJobManager().submit_intent(**request)
    assert recovered.state == "active"
    assert recovered.job_id == "100"
    assert recovered.existing is True
    assert len(fake.sbatch_calls) == 1
    entry = SlurmJobManager.ledger_entries(request["ledger_path"])[0]
    assert entry["job_id"] == "100"
    assert entry["state"] == "active"


def test_unknown_is_not_silently_resubmitted(monkeypatch, tmp_path):
    fake = FakeSlurm(state="UNKNOWN")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    manager = SlurmJobManager()
    request = _request(tmp_path)
    first = manager.submit_intent(**request)
    second = manager.submit_intent(**request)
    assert first.submitted is True
    assert second.state == "reconcile_required"
    assert len(fake.sbatch_calls) == 1


def test_dependency_is_added_to_effective_command_config(monkeypatch, tmp_path):
    fake = FakeSlurm(state="PENDING")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    result = SlurmJobManager().submit_intent(**_request(tmp_path, stage="points", dependency=("100", "101")))
    assert result.submitted is True
    assert "--dependency=afterok:100:101" in fake.sbatch_calls[0]
    entry = SlurmJobManager.ledger_entries(tmp_path / "slurm_submissions.json")[0]
    assert entry["dependency_job_ids"] == ["100", "101"]


def test_intent_key_changes_for_command_settings_manifest_and_rerun(tmp_path):
    manager = SlurmJobManager()
    base = _request(tmp_path)
    fields = (
        "system_name", "stage", "recipe_version", "settings_fingerprint",
        "manifest_fingerprint", "command", "config",
    )
    key = manager.intent_key(**{field: base[field] for field in fields})
    command_key = manager.intent_key(**{**{field: base[field] for field in (
        "system_name", "stage", "recipe_version", "settings_fingerprint",
        "manifest_fingerprint", "config",
    )}, "command": "python other.py"})
    settings_key = manager.intent_key(**{**{field: base[field] for field in (
        "system_name", "stage", "recipe_version", "manifest_fingerprint",
        "command", "config",
    )}, "settings_fingerprint": "settings-b"})
    rerun_key = manager.intent_key(**{field: base[field] for field in (
        "system_name", "stage", "recipe_version", "settings_fingerprint",
        "manifest_fingerprint", "command", "config",
    )}, run_id="run-2")
    assert len({key, command_key, settings_key, rerun_key}) == 4


def test_atomic_json_write_replaces_valid_json(tmp_path):
    path = tmp_path / "ledger.json"
    atomic_json_write(path, {"ok": True})
    assert json.loads(path.read_text()) == {"ok": True}
    assert not list(tmp_path.glob("*.tmp"))


def test_artifact_marker_rejects_stale_output_and_accepts_matching_intent(tmp_path):
    artifact = tmp_path / "result.json"
    artifact.write_text("{}")
    entry = {
        "intent_key": "intent-a",
        "system_name": "system",
        "stage": "static",
        "recipe_version": "recipe-v1",
        "settings_fingerprint": "settings-a",
        "manifest_fingerprint": "manifest-a",
    }
    assert not artifact_marker_matches(artifact, entry)
    write_completion_marker(artifact, **entry)
    assert artifact_marker_matches(artifact, entry)
    assert not artifact_marker_matches(artifact, {**entry, "settings_fingerprint": "settings-b"})


def test_legacy_ledger_is_a_reconcile_stop_not_a_submit_hint(monkeypatch, tmp_path):
    fake = FakeSlurm(state="UNKNOWN")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    path = tmp_path / "slurm_submissions.json"
    path.write_text(
        json.dumps([{"job_id": "77", "stage": "endpoints", "name": "BCC-Al-Nb"}])
    )
    result = SlurmJobManager().submit_intent(**_request(tmp_path))
    assert result.state == "reconcile_required"
    assert result.reconcile_required is True
    assert len(fake.sbatch_calls) == 0


def test_rerun_run_id_creates_a_second_explicit_intent(monkeypatch, tmp_path):
    fake = FakeSlurm(state="PENDING")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    request = _request(tmp_path, run_id="run-1")
    first = SlurmJobManager().submit_intent(**request)
    second = SlurmJobManager().submit_intent(**{**request, "run_id": "run-2"})
    assert first.submitted and second.submitted
    assert first.intent_key != second.intent_key
    assert len(fake.sbatch_calls) == 2
    entries = SlurmJobManager.ledger_entries(request["ledger_path"])
    assert [entry["run_id"] for entry in entries] == ["run-1", "run-2"]


def test_changed_command_settings_or_manifest_submit_new_intent(monkeypatch, tmp_path):
    fake = FakeSlurm(state="PENDING")
    monkeypatch.setattr("htvasp.slurm.manager.subprocess.run", fake)
    manager = SlurmJobManager()
    base = _request(tmp_path)
    manager.submit_intent(**base)
    manager.submit_intent(**{**base, "command": "python worker.py --different"})
    manager.submit_intent(**{**base, "settings_fingerprint": "settings-b"})
    manager.submit_intent(**{**base, "manifest_fingerprint": "manifest-b"})
    assert len(fake.sbatch_calls) == 4
    assert len(SlurmJobManager.ledger_entries(base["ledger_path"])) == 4

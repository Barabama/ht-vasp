"""Offline submission-contract tests for the three production portals."""

from __future__ import annotations

from pathlib import Path


class FakeResult:
    state = "submitted"
    job_id = "9001"
    intent_key = "intent-entrypoint"
    submitted = True
    existing = False
    reconcile_required = False

    @property
    def ok(self):
        return True


class FakeManager:
    all_calls = []

    def __init__(self):
        self.calls = []
        type(self).all_calls.append(self.calls)

    def get_cpu_config(self, **kwargs):
        from htvasp.slurm import SlurmConfig

        return SlurmConfig(**kwargs)

    def submit_intent(self, **kwargs):
        self.calls.append(kwargs)
        return FakeResult()


def test_qha_submission_uses_intent_api(monkeypatch):
    import main_qha

    monkeypatch.setattr(main_qha, "SlurmJobManager", FakeManager)
    main_qha.STRUCTURE_NAMES = ["SER-Al"]
    manager = FakeManager()
    result = main_qha.submit_job("SER-Al", "static", manager)
    assert result == "9001"
    call = manager.calls[-1]
    assert call["stage"] == "static"
    assert call["recipe_version"] == main_qha.QHA_RECIPE_VERSION
    assert call["manifest_fingerprint"]
    assert "--run-id" in call["command"]


def test_w90_submission_uses_intent_api_and_run_id(monkeypatch):
    import main_w90

    monkeypatch.setattr(main_w90, "SlurmJobManager", FakeManager)
    FakeManager.all_calls.clear()
    main_w90.STRUCTURE_NAMES = ["SER-Co"]
    main_w90.submit_jobs(force=False)
    call = FakeManager.all_calls[-1][-1]
    assert call["stage"] == "tb2j"
    assert call["recipe_version"] == "main_w90-v2"
    assert "--run-id" in call["command"]
    assert "--rerun" not in call["command"]


def test_oj_submission_targets_oj_not_relax(monkeypatch):
    import main_oj

    monkeypatch.setattr(main_oj, "SlurmJobManager", FakeManager)
    FakeManager.all_calls.clear()
    main_oj.STRUCTURE_NAMES = ["SER-Co"]
    main_oj.submit_jobs(rerun=False)
    call = FakeManager.all_calls[-1][-1]
    assert call["stage"] == "oj"
    assert "--oj" in call["command"]
    assert "--relax" not in call["command"]
    assert "--run-id" in call["command"]

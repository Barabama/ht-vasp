"""CLI smoke tests for the non-VASP SQS stages."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def _run(repo_root, *args):
    return subprocess.run(
        [sys.executable, "main_sqs.py", *args],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )


def test_positional_command_interface_is_accepted():
    repo_root = Path(__file__).resolve().parents[1]
    result = _run(repo_root, "generate", "--name", "BCC-Al-Nb", "--dry-run")
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["name"] == "BCC-Al-Nb"


def test_positional_command_defaults_to_batch():
    repo_root = Path(__file__).resolve().parents[1]
    result = _run(repo_root, "generate", "--dry-run")
    assert result.returncode == 0, result.stderr
    payloads = [json.loads(line) for line in result.stdout.splitlines() if line.strip()]
    assert [item["name"] for item in payloads] == [
        "BCC-Al-Nb",
        "BCC-Al-Ti",
        "BCC-Al-V",
        "BCC-Nb-Ti",
        "BCC-Nb-V",
        "BCC-Ti-V",
    ]


def test_positional_status_and_validate_accept_name():
    repo_root = Path(__file__).resolve().parents[1]
    status = _run(repo_root, "status", "--name", "BCC-Al-Nb")
    assert status.returncode == 0, status.stderr
    status_payload = json.loads(status.stdout)
    assert status_payload["points_failed"] in {0, 3}
    assert status_payload["points_successful"] + status_payload["points_pending"] + status_payload["points_failed"] == status_payload["points"]

    validate = _run(repo_root, "validate", "--name", "BCC-Al-Nb")
    assert validate.returncode == 0, validate.stderr
    assert validate.stdout.strip() == "SQS manifest valid"


def test_read_only_stage_commands_accept_name_for_default_root():
    repo_root = Path(__file__).resolve().parents[1]
    expected = repo_root / "data" / "endmembers" / "BCC-Al-Nb" / "sqsflow"
    status = _run(repo_root, "--status", "--name", "BCC-Al-Nb")
    assert status.returncode == 0, status.stderr
    assert json.loads(status.stdout)["endpoints"] == 2
    assert expected.is_dir()


def test_compact_batch_dry_run_uses_configured_structure_names():
    repo_root = Path(__file__).resolve().parents[1]
    result = _run(repo_root, "--batch", "--generate", "--dry-run")
    assert result.returncode == 0, result.stderr
    payloads = [json.loads(line) for line in result.stdout.splitlines() if line.strip()]
    assert [item["name"] for item in payloads] == [
        "BCC-Al-Nb",
        "BCC-Al-Ti",
        "BCC-Al-V",
        "BCC-Nb-Ti",
        "BCC-Nb-V",
        "BCC-Ti-V",
    ]


def test_compact_endpoint_dry_run_is_slurm_safe():
    repo_root = Path(__file__).resolve().parents[1]
    result = _run(repo_root, "--name", "BCC-Al-Nb", "--slurm", "--endpoints", "--dry-run")
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["stage"] == "endpoints"
    assert payload["submit"] is False


def test_positional_endpoint_dry_run_is_slurm_safe():
    repo_root = Path(__file__).resolve().parents[1]
    result = _run(
        repo_root,
        "endpoints",
        "--slurm",
        "--dry-run",
        "--name",
        "BCC-Al-Nb",
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    assert payload["stage"] == "endpoints"
    assert payload["submit"] is False


def test_positional_points_dry_run_does_not_construct_slurm_manager(monkeypatch):
    import main_sqs

    class UnexpectedSlurmConstruction:
        def __init__(self, *args, **kwargs):
            raise AssertionError("dry-run must not construct a Slurm submitter")

    monkeypatch.setattr("htvasp.slurm.SlurmJobManager", UnexpectedSlurmConstruction)
    assert main_sqs.main(
        ["points", "--name", "BCC-Al-Nb", "--slurm", "--dry-run"]
    ) == 0


def test_positional_collect_dry_run_is_read_only():
    repo_root = Path(__file__).resolve().parents[1]
    result = _run(repo_root, "collect", "--name", "BCC-Al-V", "--dry-run")
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == {
        "name": "BCC-Al-V",
        "action": "collect",
        "write": False,
    }


def test_compact_help_hides_internal_recipe_parameters():
    repo_root = Path(__file__).resolve().parents[1]
    result = _run(repo_root, "--help")
    assert result.returncode == 0
    assert "--global-incar" not in result.stdout
    assert "--static-incar" not in result.stdout
    assert "--batch" in result.stdout


def test_compact_recipe_is_internal():
    from main_sqs import SQS_RECIPE_VERSION, _calculation_settings, _make_plan

    system, plan, calculation = _make_plan("BCC-Al-Nb")
    _, _, _, fingerprint = _calculation_settings()
    assert system.canonical_name == "BCC-Al-Nb"
    assert plan.supercell == (2, 2, 2)
    assert calculation["recipe_version"] == SQS_RECIPE_VERSION
    assert calculation["settings_fingerprint"] == fingerprint


def test_slurm_uses_one_structure_level_log_without_error_file():
    from main_sqs import _slurm_log_path
    from htvasp.model import SqsSystem

    system = SqsSystem.from_name("BCC-Al-Nb")
    assert _slurm_log_path(system) == (
        Path(__file__).resolve().parents[1]
        / "data"
        / "endmembers"
        / "BCC-Al-Nb"
        / "slurm-sqsflow.log"
    )


def test_points_dependency_gate_rejects_unreconciled_endpoint():
    import main_sqs

    class FakeManager:
        def reconcile_ledger(self, path):
            return []

        def ledger_entries(self, path):
            return [{
                "system_name": "BCC-Al-Nb",
                "stage": "endpoints",
                "manifest_fingerprint": "manifest-a",
                "job_id": "123",
            }]

        def reconcile_entry(self, entry):
            return "active", entry["job_id"]

    with __import__("pytest").raises(Exception, match="not successfully reconciled"):
        main_sqs._endpoint_dependency_job_ids(
            FakeManager(),
            Path("/tmp/unused-ledger.json"),
            name="BCC-Al-Nb",
            manifest_fingerprint="manifest-a",
        )


def test_points_dependency_gate_returns_completed_endpoint_job_ids():
    import main_sqs

    class FakeManager:
        def reconcile_ledger(self, path):
            return []

        def ledger_entries(self, path):
            return [
                {
                    "system_name": "BCC-Al-Nb",
                    "stage": "endpoints",
                    "manifest_fingerprint": "manifest-a",
                    "job_id": "123",
                    "attempt": 1,
                }
            ]

        def reconcile_entry(self, entry):
            return "completed", entry["job_id"]

    assert main_sqs._endpoint_dependency_job_ids(
        FakeManager(),
        Path("/tmp/unused-ledger.json"),
        name="BCC-Al-Nb",
        manifest_fingerprint="manifest-a",
    ) == ["123"]


def test_reconcile_dry_run_does_not_call_scheduler(monkeypatch):
    import main_sqs

    class UnexpectedManager:
        def __init__(self, *args, **kwargs):
            raise AssertionError("dry-run reconcile must not construct scheduler manager")

    monkeypatch.setattr(main_sqs, "SlurmJobManager", UnexpectedManager)
    assert main_sqs.main(["reconcile", "--name", "BCC-Al-Nb", "--dry-run"]) == 0

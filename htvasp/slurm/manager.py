"""Slurm submission and durable intent management for HT-VASP."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shlex
import subprocess
import tempfile
import time
import uuid
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, Literal, Mapping

log = logging.getLogger(__name__)

LEDGER_SCHEMA_VERSION = 1
LEDGER_KIND = "htvasp.slurm_intent_ledger"
ACTIVE_STATES = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", "RESIZING"}
COMPLETED_STATES = {"COMPLETED"}
FAILED_STATES = {"CANCELLED", "DEADLINE", "FAILED", "NODE_FAIL", "OUT_OF_MEMORY", "OOM", "PREEMPTED", "REVOKED", "TIMEOUT"}
INTENT_STATES = {"planned", "submitting", "submitted", "active", "completed", "failed", "reconcile_required"}


def payload_fingerprint(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def new_run_id(prefix: str = "run") -> str:
    return f"{prefix}-{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}-{uuid.uuid4().hex[:10]}"


def stable_run_id(prefix: str, payload: Any) -> str:
    """Derive a repeatable run directory name for an ordinary request."""

    return f"{prefix}-{payload_fingerprint(payload)[:20]}"


def atomic_json_write(path: Path | str, payload: Any) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=True, default=str)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        try:
            directory_fd = os.open(destination.parent, os.O_DIRECTORY)
        except (AttributeError, OSError):
            directory_fd = None
        if directory_fd is not None:
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        return destination
    finally:
        if temporary.exists():
            temporary.unlink()


@contextmanager
def file_lock(path: Path | str, *, timeout: float = 60.0):
    import fcntl

    lock_path = Path(path)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as handle:
        deadline = time.monotonic() + timeout
        while True:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise TimeoutError(f"Timed out acquiring lock: {lock_path}")
                time.sleep(0.05)
        try:
            yield handle
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


@dataclass(frozen=True)
class SubmissionResult:
    intent_key: str
    state: str
    job_id: str | None = None
    submitted: bool = False
    existing: bool = False
    reconcile_required: bool = False

    @property
    def ok(self) -> bool:
        return self.state in {"submitted", "active", "completed"} and not self.reconcile_required


def completion_marker_path(artifact_path: Path | str) -> Path:
    artifact = Path(artifact_path)
    return artifact.with_name(f".{artifact.name}.intent.json")


def write_completion_marker(artifact_path: Path | str, *, intent_key: str, system_name: str, stage: str, recipe_version: str, settings_fingerprint: str, manifest_fingerprint: str, run_id: str | None = None) -> Path:
    artifact = Path(artifact_path)
    if not artifact.is_file():
        raise FileNotFoundError(f"Cannot mark missing artifact: {artifact}")
    return atomic_json_write(completion_marker_path(artifact), {
        "schema_version": LEDGER_SCHEMA_VERSION,
        "intent_key": intent_key,
        "system_name": system_name,
        "stage": stage,
        "recipe_version": recipe_version,
        "settings_fingerprint": settings_fingerprint,
        "manifest_fingerprint": manifest_fingerprint,
        "run_id": run_id,
        "artifact": str(artifact.resolve()),
        "artifact_mtime_ns": artifact.stat().st_mtime_ns,
    })


def artifact_marker_matches(artifact_path: Path | str, entry: Mapping[str, Any]) -> bool:
    artifact = Path(artifact_path)
    marker = completion_marker_path(artifact)
    if not artifact.is_file() or not marker.is_file():
        return False
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    fields = ("intent_key", "system_name", "stage", "recipe_version", "settings_fingerprint", "manifest_fingerprint")
    return all(payload.get(field) == entry.get(field) for field in fields) and payload.get("artifact") == str(artifact.resolve())


def mark_artifact_from_ledger(ledger_path: Path | str, artifact_path: Path | str, *, intent_key: str | None = None) -> Path | None:
    key = intent_key or os.environ.get("HTVASP_INTENT_KEY")
    if not key:
        return None
    for entry in SlurmJobManager.ledger_entries(ledger_path):
        if entry.get("intent_key") == key:
            return write_completion_marker(
                artifact_path,
                intent_key=key,
                system_name=str(entry.get("system_name", "")),
                stage=str(entry.get("stage", "")),
                recipe_version=str(entry.get("recipe_version", "")),
                settings_fingerprint=str(entry.get("settings_fingerprint", "")),
                manifest_fingerprint=str(entry.get("manifest_fingerprint", "")),
                run_id=entry.get("run_id"),
            )
    log.warning("No ledger entry found for artifact intent %s", key)
    return None


def current_artifact_matches(
    ledger_path: Path | str,
    artifact_path: Path | str,
    *,
    intent_key: str | None = None,
) -> bool:
    """Check an artifact against the intent propagated to a child process."""

    key = intent_key or os.environ.get("HTVASP_INTENT_KEY")
    if not key:
        return False
    for entry in SlurmJobManager.ledger_entries(ledger_path):
        if entry.get("intent_key") == key:
            return artifact_marker_matches(artifact_path, entry)
    return False


@dataclass
class SlurmConfig:
    job_name: str = "vasp-job"
    output_log: str = "job.log"
    error_log: str = ""
    memory: str = "20G"
    nodes: int = 1
    ntasks: int = 1
    ntasks_per_node: int = 0
    nodelist: str = ""
    cpus_per_task: int = 1
    gpus_per_task: int = 0
    partition: Literal["partCPU", "partGPU"] = "partCPU"
    time_limit: str = "100:00:00"
    dependency: str = ""
    vasp_cmd: str = "srun vasp_std"
    vasp_gam_cmd: str = "srun vasp_gam"
    conda_prefix: str = "/opt/miniconda3"
    conda_env: str = ""
    module_name: str | list[str] = ""
    extra_commands: list[str] = field(default_factory=list)

    def init_module(self, module_name: str | list[str] = "") -> list[str]:
        names = module_name if isinstance(module_name, list) else [module_name]
        return [". /etc/profile.d/modules.sh", "ulimit -s unlimited", "module purge"] + [f"module load {name}" for name in names if name]

    def init_conda(self, conda_env: str = "") -> list[str]:
        if conda_env.startswith("./"):
            conda_env = str(Path.cwd() / conda_env[2:])
        elif conda_env.startswith("."):
            conda_env = str(Path.cwd() / conda_env[1:])
        return [f". {self.conda_prefix}/etc/profile.d/conda.sh", f"conda activate {shlex.quote(conda_env)}"]


class SlurmJobManager:
    def __init__(self):
        self.submitted_jobs: list[dict[str, Any]] = []

    @staticmethod
    def config_fingerprint(config: SlurmConfig) -> str:
        return payload_fingerprint(asdict(config))

    @staticmethod
    def intent_key(*, system_name: str, stage: str, recipe_version: str, settings_fingerprint: str, manifest_fingerprint: str, command: str, config: SlurmConfig, conda_env: str = "", module_name: str | list[str] = "", run_id: str | None = None, retry_attempt: int | None = None) -> str:
        payload: dict[str, Any] = {
            "system_name": str(system_name), "stage": str(stage), "recipe_version": str(recipe_version),
            "settings_fingerprint": str(settings_fingerprint), "manifest_fingerprint": str(manifest_fingerprint),
            "command_fingerprint": payload_fingerprint(str(command)), "config_fingerprint": SlurmJobManager.config_fingerprint(config),
            "conda_env": str(conda_env), "module_name": module_name,
        }
        if run_id is not None:
            payload["run_id"] = str(run_id)
        if retry_attempt is not None:
            payload["retry_attempt"] = int(retry_attempt)
        return payload_fingerprint(payload)

    @staticmethod
    def _now() -> str:
        return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    @classmethod
    def _ledger_payload(cls, raw: Any) -> dict[str, Any]:
        if raw is None:
            return {"schema_version": LEDGER_SCHEMA_VERSION, "kind": LEDGER_KIND, "entries": []}
        if isinstance(raw, list):
            now = cls._now()
            entries = []
            for item in raw:
                item = dict(item)
                entries.append({**item, "legacy": True, "intent_key": item.get("intent_key"), "state": item.get("state") if item.get("state") in INTENT_STATES else "reconcile_required", "job_id": str(item["job_id"]) if item.get("job_id") else None, "created_at": item.get("created_at", now), "updated_at": item.get("updated_at", now), "attempt": int(item.get("attempt", 1)), "retry_of": item.get("retry_of")})
            return {"schema_version": LEDGER_SCHEMA_VERSION, "kind": LEDGER_KIND, "entries": entries}
        if not isinstance(raw, dict) or raw.get("kind") != LEDGER_KIND or int(raw.get("schema_version", -1)) != LEDGER_SCHEMA_VERSION:
            raise ValueError("Unsupported Slurm submission ledger")
        if not isinstance(raw.get("entries", []), list):
            raise ValueError("Slurm submission ledger entries must be a list")
        return {"schema_version": LEDGER_SCHEMA_VERSION, "kind": LEDGER_KIND, "entries": raw.get("entries", [])}

    @classmethod
    def read_ledger(cls, ledger_path: Path | str) -> dict[str, Any]:
        source = Path(ledger_path)
        if not source.exists():
            return cls._ledger_payload(None)
        try:
            return cls._ledger_payload(json.loads(source.read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            raise RuntimeError(f"Cannot safely read Slurm ledger {source}: {exc}") from exc

    @classmethod
    def ledger_entries(cls, ledger_path: Path | str) -> list[dict[str, Any]]:
        return [dict(entry) for entry in cls.read_ledger(ledger_path)["entries"]]

    @staticmethod
    def _write_ledger(ledger_path: Path | str, payload: Mapping[str, Any]) -> None:
        atomic_json_write(ledger_path, payload)

    @staticmethod
    def _upsert_entry(ledger: Mapping[str, Any], entry: Mapping[str, Any]) -> dict[str, Any]:
        entries = [dict(item) for item in ledger.get("entries", [])]
        for index, old in enumerate(entries):
            if entry.get("intent_key") and old.get("intent_key") == entry.get("intent_key"):
                entries[index] = dict(entry)
                break
        else:
            entries.append(dict(entry))
        return {"schema_version": LEDGER_SCHEMA_VERSION, "kind": LEDGER_KIND, "entries": entries}

    @staticmethod
    def _intent_family_matches(entry: Mapping[str, Any], candidate: Mapping[str, Any]) -> bool:
        fields = ("system_name", "stage", "recipe_version", "settings_fingerprint", "manifest_fingerprint", "command_fingerprint", "config_fingerprint", "conda_env", "module_name", "run_id")
        return not entry.get("legacy") and all(entry.get(field) == candidate.get(field) for field in fields)

    def _query_squeue_state(self, job_id: str) -> str | None:
        try:
            result = subprocess.run(["squeue", "--noheader", "--jobs", str(job_id), "--format=%T"], capture_output=True, text=True, check=False)
        except OSError:
            return None
        if result.returncode != 0:
            return None
        lines = result.stdout.strip().splitlines()
        return lines[0].strip().upper() if lines and lines[0].strip() else None

    def _query_sacct_state(self, job_id: str) -> str | None:
        try:
            result = subprocess.run(["sacct", "--noheader", "--parsable2", "--allocations", "--jobs", str(job_id), "--format=State"], capture_output=True, text=True, check=False)
        except OSError:
            return None
        if result.returncode != 0:
            return None
        for line in result.stdout.strip().splitlines():
            state = line.strip().split("|", 1)[0].split()[0].upper()
            if state:
                return state.split("+", 1)[0]
        return None

    def get_job_state(self, job_id: str) -> str | None:
        return self._query_squeue_state(job_id) or self._query_sacct_state(job_id)

    def find_jobs_by_intent(self, intent_key: str) -> list[dict[str, str]]:
        marker = f"htvasp-intent:{intent_key}"
        found: list[dict[str, str]] = []
        try:
            result = subprocess.run(["squeue", "--noheader", "--format=%i|%T|%k"], capture_output=True, text=True, check=False)
        except OSError:
            result = None
        if result is not None and result.returncode == 0:
            for line in result.stdout.strip().splitlines():
                fields = line.strip().split("|", 2)
                if len(fields) == 3 and fields[2].strip() == marker:
                    found.append({"job_id": fields[0].strip(), "state": fields[1].strip().upper()})
        try:
            result = subprocess.run(["sacct", "--noheader", "--parsable2", "--allocations", "--starttime=1970-01-01", "--format=JobIDRaw,State,Comment"], capture_output=True, text=True, check=False)
        except OSError:
            result = None
        if result is not None and result.returncode == 0:
            seen = {item["job_id"] for item in found}
            for line in result.stdout.strip().splitlines():
                fields = line.strip().split("|", 2)
                if len(fields) == 3 and fields[2].strip() == marker:
                    job_id = fields[0].strip().split(".", 1)[0]
                    if job_id and job_id not in seen:
                        found.append({"job_id": job_id, "state": fields[1].strip().upper()})
                        seen.add(job_id)
        return found

    def reconcile_entry(self, entry: Mapping[str, Any]) -> tuple[str, str | None]:
        raw_job_id = entry.get("job_id")
        job_id = str(raw_job_id).strip() if raw_job_id is not None else None
        if not job_id:
            job_id = None
        if not job_id:
            if entry.get("state") == "submitting" and entry.get("intent_key"):
                matches = self.find_jobs_by_intent(str(entry["intent_key"]))
                if len(matches) != 1:
                    return "reconcile_required", None
                job_id = matches[0]["job_id"]
                state = matches[0].get("state") or self.get_job_state(job_id)
            elif entry.get("state") in {"completed", "failed"}:
                return str(entry["state"]), None
            else:
                return "reconcile_required", None
        else:
            state = self.get_job_state(job_id)
        if state is None:
            return "reconcile_required", job_id
        normalized = str(state).upper().split("+", 1)[0]
        if normalized in ACTIVE_STATES:
            return "active", job_id
        if normalized in COMPLETED_STATES:
            return "completed", job_id
        if normalized in FAILED_STATES or normalized.startswith(("CANCELLED", "FAILED")):
            return "failed", job_id
        return "reconcile_required", job_id

    def reconcile_ledger(self, ledger_path: Path | str, *, intent_key: str | None = None, artifact_checker: Any = None) -> list[dict[str, Any]]:
        path = Path(ledger_path).resolve()
        with file_lock(path.with_name(f".{path.name}.lock")):
            ledger = self.read_ledger(path)
            entries = []
            for entry in ledger["entries"]:
                if intent_key and entry.get("intent_key") != intent_key:
                    entries.append(entry)
                    continue
                state, job_id = self.reconcile_entry(entry)
                updated = dict(entry, state=state, updated_at=self._now())
                if job_id:
                    updated["job_id"] = job_id
                if state == "completed" and artifact_checker is not None and not artifact_checker(updated):
                    updated["state"] = "reconcile_required"
                entries.append(updated)
            self._write_ledger(path, {**ledger, "entries": entries})
            return entries

    def submit_intent(self, *, system_name: str, stage: str, recipe_version: str, settings_fingerprint: str, manifest_fingerprint: str, command: str, config: SlurmConfig, ledger_path: Path | str, conda_env: str = "", module_name: str | list[str] = "", work_dir: Path | str = "", run_id: str | None = None, retry: bool = False, dependency_job_ids: list[str] | tuple[str, ...] = (), artifact_checker: Any = None, artifact_shortcut: bool = False, artifact_path: Path | str | None = None) -> SubmissionResult:
        if not all(str(value).strip() for value in (system_name, stage, recipe_version, settings_fingerprint, manifest_fingerprint)):
            raise ValueError("intent identity fields must be non-empty")
        path = Path(ledger_path).resolve()
        dependencies = tuple(str(item) for item in dependency_job_ids if str(item))
        effective_config = replace(config, dependency=f"afterok:{':'.join(dependencies)}") if dependencies and not config.dependency else config

        def make_key(attempt: int | None = None) -> str:
            return self.intent_key(system_name=system_name, stage=stage, recipe_version=recipe_version, settings_fingerprint=settings_fingerprint, manifest_fingerprint=manifest_fingerprint, command=command, config=effective_config, conda_env=conda_env, module_name=module_name, run_id=run_id, retry_attempt=attempt)

        requested_key = make_key()
        with file_lock(path.with_name(f".{path.name}.lock")):
            ledger = self.read_ledger(path)
            if run_id is None:
                legacy = next((entry for entry in ledger["entries"] if entry.get("legacy") and (entry.get("system_name") or entry.get("name")) == system_name and entry.get("stage") == stage), None)
                if legacy is not None:
                    state, job_id = self.reconcile_entry(legacy)
                    return SubmissionResult(requested_key, state if state == "active" else "reconcile_required", job_id, existing=True, reconcile_required=state != "active")
            identity = {"system_name": system_name, "stage": stage, "recipe_version": recipe_version, "settings_fingerprint": settings_fingerprint, "manifest_fingerprint": manifest_fingerprint, "command_fingerprint": payload_fingerprint(command), "config_fingerprint": self.config_fingerprint(effective_config), "conda_env": conda_env, "module_name": module_name, "run_id": run_id}
            family = [entry for entry in ledger["entries"] if self._intent_family_matches(entry, identity)]
            matching = max(family, key=lambda entry: int(entry.get("attempt", 1)), default=None)
            retry_of = None
            if matching is not None:
                state, job_id = self.reconcile_entry(matching)
                updated = dict(matching, state=state, updated_at=self._now())
                if job_id:
                    updated["job_id"] = job_id
                if state == "completed" and artifact_checker is not None and not artifact_checker(updated):
                    state = "reconcile_required"
                    updated["state"] = state
                self._write_ledger(path, self._upsert_entry(ledger, updated))
                if state in {"active", "submitted", "completed"}:
                    return SubmissionResult(requested_key, state, updated.get("job_id"), existing=True)
                if state == "reconcile_required":
                    return SubmissionResult(requested_key, state, updated.get("job_id"), existing=True, reconcile_required=True)
                if state == "failed" and not retry:
                    return SubmissionResult(requested_key, state, updated.get("job_id"), existing=True)
                attempt = int(updated.get("attempt", 1)) + 1
                requested_key = make_key(attempt)
                retry_of = updated.get("intent_key")
            else:
                attempt = 1
                candidate = {"intent_key": requested_key, **identity, "command": command, "config": asdict(effective_config), "attempt": attempt, "artifact_path": str(artifact_path) if artifact_path is not None else None}
                if artifact_shortcut and artifact_checker is not None and artifact_checker(candidate):
                    now = self._now()
                    self._write_ledger(path, self._upsert_entry(ledger, {**candidate, "state": "completed", "job_id": None, "created_at": now, "updated_at": now, "retry_of": None, "dependency_job_ids": list(dependencies)}))
                    return SubmissionResult(requested_key, "completed", existing=True)
            now = self._now()
            entry = {"schema_version": LEDGER_SCHEMA_VERSION, "intent_key": requested_key, **identity, "command": command, "command_fingerprint": payload_fingerprint(command), "config_fingerprint": self.config_fingerprint(effective_config), "config": asdict(effective_config), "conda_env": conda_env, "module_name": module_name, "run_id": run_id, "state": "submitting", "job_id": None, "created_at": now, "updated_at": now, "attempt": attempt, "retry_of": retry_of, "dependency_job_ids": list(dependencies), "artifact_path": str(artifact_path) if artifact_path is not None else None}
            entry["state"] = "planned"
            self._write_ledger(path, self._upsert_entry(ledger, entry))
            entry["state"] = "submitting"
            entry["updated_at"] = self._now()
            self._write_ledger(path, self._upsert_entry(self.read_ledger(path), entry))
            try:
                job_id = self.submit_command(
                    command=command,
                    config=effective_config,
                    conda_env=conda_env,
                    module_name=module_name,
                    intent_key=requested_key,
                    work_dir=work_dir,
                )
            except Exception as exc:
                # sbatch may have accepted the job before the client observed
                # an exception.  Keep the durable ``submitting`` marker so a
                # later call can find the job by its Slurm comment.
                log.exception(
                    "Submission call failed for intent %s; reconciliation required: %s",
                    requested_key,
                    exc,
                )
                return SubmissionResult(
                    requested_key,
                    "reconcile_required",
                    reconcile_required=True,
                )
            if not job_id:
                entry["state"] = "failed"
                entry["updated_at"] = self._now()
                self._write_ledger(path, self._upsert_entry(self.read_ledger(path), entry))
                return SubmissionResult(requested_key, "failed")
            entry.update({"job_id": str(job_id), "state": "submitted", "updated_at": self._now()})
            try:
                self._write_ledger(path, self._upsert_entry(self.read_ledger(path), entry))
            except Exception as exc:
                log.critical("sbatch succeeded for intent %s but ledger persistence failed: %s", requested_key, exc)
                return SubmissionResult(requested_key, "reconcile_required", str(job_id), submitted=True, reconcile_required=True)
            return SubmissionResult(requested_key, "submitted", str(job_id), submitted=True)

    def submit_command(self, command: str, config: SlurmConfig | None = None, conda_env: str = "", module_name: str | list[str] = "", work_dir: Path | str = "", intent_key: str | None = None) -> str | None:
        config = config or self._default_cpu_config()
        working_directory = Path(work_dir or Path.cwd()).resolve()
        working_directory.mkdir(parents=True, exist_ok=True)
        effective_conda = conda_env or config.conda_env
        effective_module = module_name or config.module_name
        init = config.init_module(effective_module) if effective_module else []
        if effective_conda:
            init += config.init_conda(effective_conda)
        init += config.extra_commands
        if intent_key:
            init.append(f"export HTVASP_INTENT_KEY={shlex.quote(intent_key)}")
        full_command = "; ".join(init + [command])
        cmd = ["sbatch", f"--job-name={config.job_name}", f"--output={config.output_log}", f"--mem={config.memory}", f"--nodes={config.nodes}", f"--ntasks={config.ntasks}", f"--cpus-per-task={config.cpus_per_task}", f"--partition={config.partition}", f"--time={config.time_limit}", f"--wrap={full_command}"]
        if config.error_log:
            cmd.append(f"--error={config.error_log}")
        if config.ntasks_per_node:
            cmd.append(f"--ntasks-per-node={config.ntasks_per_node}")
        if config.nodelist:
            cmd.append(f"--nodelist={config.nodelist}")
        if config.dependency:
            cmd.append(f"--dependency={config.dependency}")
        if config.gpus_per_task > 0:
            cmd.append(f"--gpus-per-task={config.gpus_per_task}")
        if intent_key:
            cmd.append(f"--comment=htvasp-intent:{intent_key}")
        try:
            result = subprocess.run(cmd, cwd=working_directory, capture_output=True, text=True, check=True)
        except (OSError, subprocess.CalledProcessError) as exc:
            log.error("sbatch failed: %s", exc)
            return None
        for line in result.stdout.strip().splitlines():
            if "Submitted batch job" in line:
                job_id = line.split()[-1]
                self.submitted_jobs.append({"job_id": job_id, "command": command, "config": config, "work_dir": str(working_directory), "intent_key": intent_key})
                return job_id
        if intent_key:
            raise RuntimeError(
                "sbatch returned success but no job ID was parsed; reconciliation required"
            )
        return None

    def _default_cpu_config(self) -> SlurmConfig:
        return SlurmConfig(job_name="vasp-cpu", output_log="job.log", ntasks=48, memory="96G", partition="partCPU")

    def get_cpu_config(self, **kwargs: Any) -> SlurmConfig:
        config = self._default_cpu_config()
        for key, value in kwargs.items():
            if hasattr(config, key):
                setattr(config, key, value)
        return config

    def get_gpu_config(self, **kwargs: Any) -> SlurmConfig:
        return SlurmConfig(job_name=kwargs.pop("job_name", "vasp-gpu"), output_log=kwargs.pop("output_log", "job.log"), ntasks=kwargs.pop("ntasks", 1), gpus_per_task=kwargs.pop("gpus_per_task", 1), memory=kwargs.pop("memory", "20G"), partition="partGPU", module_name=kwargs.pop("module_name", "vasp-gpu"), **kwargs)

    def get_job_status(self, job_id: str) -> dict[str, Any] | None:
        try:
            result = subprocess.run(["scontrol", "show", "job", str(job_id)], capture_output=True, text=True, check=True)
        except (OSError, subprocess.CalledProcessError):
            return None
        status: dict[str, Any] = {}
        for line in result.stdout.strip().splitlines():
            for item in line.split():
                if "=" in item:
                    key, value = item.split("=", 1)
                    status[key] = value
        return status

    def is_job_completed(self, job_id: str) -> bool:
        state = self.get_job_state(job_id)
        return state in COMPLETED_STATES | FAILED_STATES if state is not None else False

    def cancel_job(self, job_id: str) -> bool:
        try:
            subprocess.run(["scancel", str(job_id)], capture_output=True, text=True, check=True)
            return True
        except (OSError, subprocess.CalledProcessError):
            return False

    def list_jobs(self, user: str = "") -> list[dict[str, Any]]:
        import getpass
        try:
            result = subprocess.run(["squeue", "-u", user or getpass.getuser()], capture_output=True, text=True, check=True)
        except (OSError, subprocess.CalledProcessError):
            return []
        lines = result.stdout.strip().splitlines()
        if not lines:
            return []
        header = lines[0].split()
        return [dict(zip(header, line.split())) for line in lines[1:]]


__all__ = ["ACTIVE_STATES", "COMPLETED_STATES", "FAILED_STATES", "INTENT_STATES", "LEDGER_KIND", "LEDGER_SCHEMA_VERSION", "SlurmConfig", "SlurmJobManager", "SubmissionResult", "artifact_marker_matches", "atomic_json_write", "completion_marker_path", "current_artifact_matches", "file_lock", "mark_artifact_from_ledger", "new_run_id", "payload_fingerprint", "stable_run_id", "write_completion_marker"]

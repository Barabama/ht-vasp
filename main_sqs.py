#!/usr/bin/env python3
"""SQS-mixing entry point for the HT-VASP endmember portal.

Recommended usage:
    python main_sqs.py generate
    python main_sqs.py endpoints --slurm
    python main_sqs.py points --slurm
    python main_sqs.py collect
    python main_sqs.py status
    python main_sqs.py validate
    python main_sqs.py --name BCC-Al-Nb --generate
    python main_sqs.py --batch --generate
    python main_sqs.py --batch --slurm --endpoints
    python main_sqs.py --batch --slurm --points
    python main_sqs.py --batch --collect

The calculation recipe is intentionally configured in this file, following
the ``main_qha.py``/``main_w90.py`` pattern.  This is the only supported CLI;
stale manifests must be regenerated with ``--rerun``.
"""

from __future__ import annotations

import argparse
import json
import logging
import shlex
import sys
from pathlib import Path
from typing import Any, Mapping

from pymatgen.core import Structure

from htvasp.model import (
    SqsPlan,
    SqsSystem,
    TemplateRecord,
    TemplateSnapshot,
    payload_checksum,
    structure_checksum,
)
from htvasp.sqs.orchestrator import SqsOrchestrationError, SqsOrchestrator
from htvasp.slurm import SlurmJobManager, atomic_json_write, file_lock, new_run_id
from htvasp.sqs.task_doc import SqsStaticObservation

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger("main_sqs")

PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_ENDMEMBER_ROOT = PROJECT_ROOT / "data" / "endmembers"
DEFAULT_TEMPLATE_ROOT = PROJECT_ROOT / "data" / "poscars"

# First production scope: binary BCC_A2 systems in the Al-Nb-Ti-V database.
# FCC/HCP and magnetic Co-Fe-Mn-Ni systems require a separate CALPHAD mapping
# decision and are deliberately not included in the default SQS batch.
STRUCTURE_NAMES = [
    "BCC-Al-Nb",
    "BCC-Al-Ti",
    "BCC-Al-V",
    "BCC-Nb-Ti",
    "BCC-Nb-V",
    "BCC-Ti-V",
]

PORTAL_ACTIONS = (
    "generate",
    "endpoints",
    "points",
    "collect",
    "status",
    "validate",
    "reconcile",
)

SQS_ENGINE = "sqsgenerator"
SQS_ITERATIONS = 10000
SQS_SEED = 0
SQS_POTCAR_FUNCTIONAL = "PBE"
SQS_RELAX_RECIPROCAL_DENSITY = 64
SQS_STATIC_RECIPROCAL_DENSITY = 64
SQS_FLOW_ROOT = "/nfs_ssd/tmp"
SQS_VASP_COMMAND = "srun vasp_std"
SQS_CONDA_ENV = str(PROJECT_ROOT / ".conda")
SQS_MODULE = "vasp-cpu"
SQS_SLURM_PARTITION = "partCPU"
SQS_SLURM_NTASKS = 32
SQS_SLURM_MEMORY = "64G"
SQS_SLURM_TIME = "48:00:00"
SQS_RECIPE_VERSION = "main_sqs-v3"
SQS_SLURM_LOG_NAME = "slurm-sqsflow.log"

DEFAULT_GLOBAL_INCAR: dict[str, Any] = {
    "ENCUT": 400,
    "ISTART": 0,
    "ICHARG": 2,
    "PREC": "Normal",
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ALGO": "Fast",
    "NELM": 100,
    "NELMIN": 6,
    "NELMDL": -6,
    "IBRION": 2,
    "ISIF": 2,
    "NSW": 10,
    "POTIM": 0.2,
    "EDIFF": 1e-5,
    "EDIFFG": -0.05,
    "ISPIN": 2,
    "ISYM": 0,
    "LREAL": "Auto",
    "SYMPREC": 1e-5,
    "LWAVE": False,
    "LCHARG": False,
    "LORBIT": None,
    "LOPTICS": False,
    "LVTOT": False,
    "LORBIT": 10,
    "GGA": "PE",
}
DEFAULT_RELAX_INCAR: dict[str, Any] = {
    "NSW": 80,
    "EDIFFG": -0.05,
}
DEFAULT_STATIC_INCAR: dict[str, Any] = {
    "ENCUT": 520,
    "PREC": "Accurate",
    "EDIFF": 1e-6,
    "LORBIT": 10,
}


def _batch_base_root(system: SqsSystem) -> Path:
    return (DEFAULT_ENDMEMBER_ROOT / system.canonical_name / "sqsflow").resolve()


def _batch_root(system: SqsSystem, run_id: str | None = None) -> Path:
    """Return an immutable run root; older runs are never removed."""

    base = _batch_base_root(system)
    if run_id:
        run_component = str(run_id).strip()
        if not run_component or Path(run_component).name != run_component:
            raise SqsOrchestrationError(f"Invalid SQS run_id: {run_id!r}")
        return base / "runs" / run_component
    pointer = base / "current_run.json"
    if pointer.is_file():
        try:
            payload = json.loads(pointer.read_text(encoding="utf-8"))
            candidate = (base / Path(str(payload["relative_root"]))).resolve()
            candidate.relative_to(base)
            if (candidate / "batch.json").is_file():
                return candidate
        except (OSError, KeyError, TypeError, ValueError):
            log.warning("Ignoring invalid SQS current-run pointer: %s", pointer)
    return base


def _write_current_run(system: SqsSystem, run_id: str, root: Path) -> None:
    base = _batch_base_root(system)
    atomic_json_write(
        base / "current_run.json",
        {
            "schema_version": 1,
            "run_id": run_id,
            "relative_root": str(root.relative_to(base)),
        },
    )


def _slurm_log_path(system: SqsSystem, run_id: str | None = None) -> Path:
    """Return the single combined stdout/stderr log for one SQS structure."""

    if run_id:
        return _batch_root(system, run_id) / SQS_SLURM_LOG_NAME
    return _batch_base_root(system).parent / SQS_SLURM_LOG_NAME


def _template_snapshot(system: SqsSystem) -> TemplateSnapshot:
    path = DEFAULT_TEMPLATE_ROOT / f"{system.canonical_name}.vasp"
    if not path.is_file():
        raise FileNotFoundError(
            f"Local template does not exist: {path}"
        )
    structure = Structure.from_file(path)
    record = TemplateRecord(
        phase=system.phase,
        path=str(path),
        source_species=tuple(str(site.specie) for site in structure),
        structure_checksum=structure_checksum(structure),
    )
    return TemplateSnapshot(record=record, structure=structure)


def _calculation_settings() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], str]:
    global_incar = dict(DEFAULT_GLOBAL_INCAR)
    relax_incar = dict(DEFAULT_RELAX_INCAR)
    static_incar = dict(DEFAULT_STATIC_INCAR)
    effective_relax = {**global_incar, "LWAVE": True, **relax_incar, "ISIF": 3}
    effective_static = {
        **global_incar,
        "ISTART": 0,
        "ICHARG": 2,
        "IBRION": -1,
        "NSW": 0,
        "ALGO": "Normal",
        "NELM": 200,
        "LWAVE": False,
        "LCHARG": False,
        "LAECHG": False,
        **static_incar,
    }
    payload = {
        "workflow": "relax_then_static",
        "potcar_functional": SQS_POTCAR_FUNCTIONAL,
        "effective_relax_incar": effective_relax,
        "effective_static_incar": effective_static,
        "relax_reciprocal_density": SQS_RELAX_RECIPROCAL_DENSITY,
        "static_reciprocal_density": SQS_STATIC_RECIPROCAL_DENSITY,
        "recipe_version": SQS_RECIPE_VERSION,
    }
    return global_incar, relax_incar, static_incar, payload_checksum(payload)


def _calculation_provenance(settings_fingerprint: str) -> dict[str, Any]:
    global_incar, relax_incar, static_incar, _ = _calculation_settings()
    return {
        "settings_fingerprint": settings_fingerprint,
        "recipe_version": SQS_RECIPE_VERSION,
        "workflow": {
            "worker": "SqsMixingWorker",
            "flow": "relax_then_static",
            "relax_maker": "RelaxMaker",
            "static_maker": "StaticMaker",
        },
        "global_incar": global_incar,
        "relax_incar": {**global_incar, "LWAVE": True, **relax_incar, "ISIF": 3},
        "static_incar": {
            **global_incar,
            "ISTART": 0,
            "ICHARG": 2,
            "IBRION": -1,
            "NSW": 0,
            "ALGO": "Normal",
            "NELM": 200,
            "LWAVE": False,
            "LCHARG": False,
            "LAECHG": False,
            **static_incar,
        },
        "relax_incar_overrides": relax_incar,
        "static_incar_overrides": static_incar,
        "kpoints": {
            "relax_reciprocal_density": SQS_RELAX_RECIPROCAL_DENSITY,
            "static_reciprocal_density": SQS_STATIC_RECIPROCAL_DENSITY,
        },
        "vasp": {"command": SQS_VASP_COMMAND},
        "potcar": {"functional": SQS_POTCAR_FUNCTIONAL},
    }


def _make_plan(name: str) -> tuple[SqsSystem, SqsPlan, dict[str, Any]]:
    system = SqsSystem.from_name(name)
    template = _template_snapshot(system)
    _, _, _, settings_fingerprint = _calculation_settings()
    plan = SqsPlan.default(
        system,
        template,
        supercell=None,
        iterations=SQS_ITERATIONS,
        seed=SQS_SEED,
        engine=SQS_ENGINE,
        engine_config={},
    )
    return system, plan, _calculation_provenance(settings_fingerprint)


def _structure_for_spec(root: Path, spec) -> Structure:
    path = Path(spec.structure.path)
    if not path.is_absolute():
        path = root / path
    return Structure.from_file(path)


def _run_spec(
    root: Path,
    spec,
    calculation: Mapping[str, Any] | None = None,
    *,
    rerun: bool = False,
) -> Mapping[str, Any] | SqsStaticObservation:
    from htvasp.workflows import SqsMixingWorker

    calculation = dict(calculation or {})
    global_incar = dict(calculation.get("global_incar", {}))
    relax_incar = dict(calculation.get("relax_incar_overrides", {}))
    static_incar = dict(calculation.get("static_incar_overrides", {}))
    if not global_incar or not relax_incar or not static_incar:
        fallback_global, fallback_relax, fallback_static, fallback_fp = _calculation_settings()
        global_incar = global_incar or fallback_global
        relax_incar = relax_incar or fallback_relax
        static_incar = static_incar or fallback_static
        if calculation.get("settings_fingerprint") and calculation["settings_fingerprint"] != fallback_fp:
            raise SqsOrchestrationError(
                "Execution settings do not match the batch fingerprint; "
                "regenerate with the intended settings or provide manifest provenance"
            )
    kpoints = dict(calculation.get("kpoints", {}))
    vasp_command = str(calculation.get("vasp", {}).get("command", SQS_VASP_COMMAND))
    potcar_functional = str(
        calculation.get("potcar", {}).get("functional", SQS_POTCAR_FUNCTIONAL)
    )
    kind = spec.kind
    task_root = root / "runs" / kind / spec.calculation_id
    flow_root = Path(SQS_FLOW_ROOT).resolve()
    worker = SqsMixingWorker(
        worker_name=spec.calculation_id,
        vasp_args={"vasp_cmd": vasp_command},
        potcar_functional=potcar_functional,
        global_incar=global_incar,
        relax_incar=relax_incar,
        static_incar=static_incar,
        relax_reciprocal_density=kpoints.get(
            "relax_reciprocal_density", SQS_RELAX_RECIPROCAL_DENSITY
        ),
        static_reciprocal_density=kpoints.get(
            "static_reciprocal_density", SQS_STATIC_RECIPROCAL_DENSITY
        ),
        calculation_spec=spec,
    )
    worker.run_flow(
        name=spec.calculation_id,
        structure=_structure_for_spec(root, spec),
        flow_dir=flow_root,
        store_dir=task_root,
        resume=not rerun,
        ensure_success=True,
        raise_immediately=False,
    )
    output = worker.get_result()
    if not output:
        raise RuntimeError(f"No observation output for {spec.calculation_id}")
    return output


def _selected_names(args: argparse.Namespace) -> list[str]:
    """Resolve one requested system or the configured production batch."""

    if args.name:
        return [SqsSystem.from_name(args.name).canonical_name]
    # The compact command interface is batch-first.  Keep --batch as an
    # accepted compatibility flag, but do not require it for new commands.
    return list(STRUCTURE_NAMES)


def _current_calculation() -> dict[str, Any]:
    _, _, _, settings_fingerprint = _calculation_settings()
    return _calculation_provenance(settings_fingerprint)


def _assert_current_batch(manifest: Any, name: str) -> None:
    expected = _current_calculation()
    actual = dict(manifest.calculation)
    if actual.get("recipe_version") != SQS_RECIPE_VERSION:
        raise SqsOrchestrationError(
            f"SQS batch for {name} uses an old recipe; regenerate with --rerun"
        )
    if actual.get("settings_fingerprint") != expected["settings_fingerprint"]:
        raise SqsOrchestrationError(
            f"SQS batch for {name} has a conflicting calculation fingerprint; regenerate with --rerun"
        )


def _portal_generate(
    name: str,
    *,
    dry_run: bool = False,
    rerun: bool = False,
    run_id: str | None = None,
) -> int:
    system = SqsSystem.from_name(name)
    if rerun:
        run_id = run_id or new_run_id("sqs")
    root = _batch_root(system, run_id)
    if rerun and run_id and (root / "batch.json").exists():
        raise SqsOrchestrationError(
            f"SQS run already exists: {root}; choose a new run_id for --rerun"
        )
    if dry_run:
        print(json.dumps({"name": name, "root": str(root), "action": "generate"}, separators=(",", ":")))
        return 0
    with file_lock(_batch_base_root(system) / ".generate.lock"):
        orchestrator = SqsOrchestrator(root)
        if orchestrator.manifest_path.exists() and not rerun:
            _assert_current_batch(orchestrator.load(), name)
            errors = orchestrator.validate()
            if errors:
                raise SqsOrchestrationError(
                    f"Existing SQS batch is invalid for {name}: " + "; ".join(errors)
                )
            log.info("SQS batch for %s already generated; use --rerun to rebuild", name)
            return 0
        _, plan, calculation = _make_plan(name)
        if run_id:
            calculation = {**calculation, "run_id": run_id}
        manifest = orchestrator.generate(
            [plan],
            settings_fingerprint=calculation["settings_fingerprint"],
            calculation=calculation,
            overwrite=rerun,
        )
        if run_id:
            _write_current_run(system, run_id, root)
    print(json.dumps({"name": name, "root": str(root), **manifest.status()}, indent=2))
    return 0


def _portal_execute(
    name: str,
    stage: str,
    *,
    dry_run: bool = False,
    rerun: bool = False,
    run_id: str | None = None,
) -> int:
    system = SqsSystem.from_name(name)
    if rerun and run_id and (_batch_root(system, run_id) / "batch.json").exists():
        raise SqsOrchestrationError(
            f"SQS run already exists: {_batch_root(system, run_id)}; choose a new run_id"
        )
    if rerun and not run_id:
        # ``--rerun`` is an explicit new run request.  Generate a fresh
        # immutable manifest first; this never removes the previous run.
        _portal_generate(name, rerun=True)
        run_id = None
        rerun = False
    root = _batch_root(system, run_id)
    orchestrator = SqsOrchestrator(root)
    manifest = orchestrator.load()
    _assert_current_batch(manifest, name)
    if dry_run:
        specs = orchestrator.pending_endpoints() if stage == "endpoints" else orchestrator.pending_points()
        print(json.dumps({"name": name, "stage": stage, "pending": [spec.calculation_id for spec in specs]}, indent=2))
        return 0
    if stage == "points":
        specs = orchestrator.pending_points()
    else:
        specs = orchestrator.pending_endpoints()
    for spec in specs:
        try:
            output = _run_spec(root, spec, manifest.calculation, rerun=rerun)
            if isinstance(output, SqsStaticObservation):
                output = output.model_copy(
                    update={"batch_fingerprint": manifest.batch_fingerprint}
                )
            else:
                output = {
                    **dict(output),
                    "batch_fingerprint": manifest.batch_fingerprint,
                }
            if stage == "endpoints":
                orchestrator.mark_endpoint_success(spec.metadata["cache_key"], output)
            else:
                orchestrator.mark_point_success(spec.calculation_id, output)
        except Exception as exc:
            log.exception("%s failed for %s", stage, spec.calculation_id)
            if stage == "endpoints":
                orchestrator.mark_endpoint_failure(spec.metadata["cache_key"], str(exc))
            else:
                orchestrator.mark_point_failure(spec.calculation_id, str(exc))
    status = orchestrator.status()
    print(json.dumps({"name": name, **status}, indent=2))
    return 1 if status["endpoints_failed"] or status["points_failed"] else 0


def _portal_collect(name: str, *, dry_run: bool = False, run_id: str | None = None) -> int:
    system = SqsSystem.from_name(name)
    orchestrator = SqsOrchestrator(_batch_root(system, run_id))
    manifest = orchestrator.load()
    _assert_current_batch(manifest, name)
    if dry_run:
        print(
            json.dumps(
                {"name": name, "action": "collect", "write": False},
                separators=(",", ":"),
            )
        )
        return 0
    data = orchestrator.collect(system.canonical_name)
    print(json.dumps({"name": name, "dataset_id": data.dataset_id, "points": len(data.points)}, indent=2))
    return 0


def _portal_status(name: str, *, validate: bool = False, run_id: str | None = None) -> int:
    orchestrator = SqsOrchestrator(_batch_root(SqsSystem.from_name(name), run_id))
    manifest = orchestrator.load()
    _assert_current_batch(manifest, name)
    if validate:
        errors = orchestrator.validate()
        if errors:
            for error in errors:
                print(f"ERROR: {error}", file=sys.stderr)
            return 1
        print("SQS manifest valid")
        return 0
    print(json.dumps({"name": name, **orchestrator.status()}, indent=2))
    return 0


def _portal_reconcile(name: str, *, run_id: str | None = None) -> int:
    """Reconcile Slurm state for this run without submitting a job."""

    system = SqsSystem.from_name(name)
    root = _batch_root(system, run_id)
    ledger_path = root / "slurm_submissions.json"
    manager = SlurmJobManager()
    entries = manager.reconcile_ledger(ledger_path)
    print(
        json.dumps(
            {
                "name": name,
                "root": str(root),
                "entries": [
                    {
                        "intent_key": entry.get("intent_key"),
                        "stage": entry.get("stage"),
                        "job_id": entry.get("job_id"),
                        "state": entry.get("state"),
                    }
                    for entry in entries
                ],
            },
            indent=2,
        )
    )
    return 0 if all(entry.get("state") != "reconcile_required" for entry in entries) else 1


def _portal_submit(
    name: str,
    stage: str,
    *,
    rerun: bool = False,
    dry_run: bool = False,
    run_id: str | None = None,
) -> int:
    """Submit one stage job for one structure using the internal recipe."""

    system = SqsSystem.from_name(name)
    if rerun and run_id and (_batch_root(system, run_id) / "batch.json").exists():
        raise SqsOrchestrationError(
            f"SQS run already exists: {_batch_root(system, run_id)}; choose a new run_id"
        )
    if rerun:
        # A rerun is a new immutable batch, not an overwrite of the current
        # run.  Generate the manifest first; the actual stage submit then
        # uses the new pointer and a normal (non-rerun) child invocation.
        _portal_generate(name, rerun=True, run_id=run_id)
        run_id = None
        rerun = False
    root = _batch_root(system, run_id)
    orchestrator = SqsOrchestrator(root)
    manifest = orchestrator.load()
    _assert_current_batch(manifest, name)
    if stage == "points":
        manifest.assert_points_ready(root)
    else:
        manifest.assert_valid(root)

    if dry_run:
        print(json.dumps({"name": name, "stage": stage, "root": str(root), "submit": False}, separators=(",", ":")))
        return 0

    script = Path(__file__).resolve()
    command_args = [
        sys.executable,
        str(script),
        "--name",
        system.canonical_name,
        f"--{stage}",
    ]
    manifest_run_id = manifest.calculation.get("run_id")
    if run_id and manifest_run_id and str(run_id) != str(manifest_run_id):
        raise SqsOrchestrationError(
            f"Requested run_id {run_id!r} does not match manifest run_id {manifest_run_id!r}"
        )
    effective_run_id = str(manifest_run_id or run_id) if (manifest_run_id or run_id) else None
    if effective_run_id:
        command_args.extend(["--run-id", effective_run_id])
    command = " ".join(shlex.quote(item) for item in command_args)
    manager = SlurmJobManager()
    config = manager.get_cpu_config(
        job_name=f"{system.canonical_name}-sqs-{stage}",
        output_log=str(_slurm_log_path(system, effective_run_id)),
        ntasks=SQS_SLURM_NTASKS,
        memory=SQS_SLURM_MEMORY,
        time_limit=SQS_SLURM_TIME,
        partition=SQS_SLURM_PARTITION,
    )
    record_path = root / "slurm_submissions.json"

    def artifact_complete(entry: Mapping[str, Any]) -> bool:
        """Only skip a completed Slurm record when this batch is complete."""

        if entry.get("manifest_fingerprint") != manifest.batch_fingerprint:
            return False
        current = SqsOrchestrator(root).load()
        if current.batch_fingerprint != manifest.batch_fingerprint:
            return False
        records = current.endpoints if stage == "endpoints" else tuple(
            point for system_record in current.systems for point in system_record.points
        )
        if not records or any(record.state != "successful" for record in records):
            return False
        for record in records:
            observation_path = getattr(record, "observation_path", None)
            if not observation_path:
                return False
            path = root / observation_path
            if not path.is_file():
                return False
            try:
                observation = SqsStaticObservation.model_validate(
                    json.loads(path.read_text(encoding="utf-8"))
                )
            except Exception:
                return False
            if observation.batch_fingerprint != manifest.batch_fingerprint:
                return False
        return True

    dependency_job_ids: list[str] = []
    if stage == "points":
        # The manifest gate above is authoritative.  Carry successful
        # endpoint job IDs as afterok dependencies when available, so a
        # points job cannot start before the endpoint job's Slurm completion.
        dependency_job_ids = _endpoint_dependency_job_ids(
            manager,
            record_path,
            name=name,
            manifest_fingerprint=manifest.batch_fingerprint,
        )


    result = manager.submit_intent(
        system_name=system.canonical_name,
        stage=stage,
        recipe_version=SQS_RECIPE_VERSION,
        settings_fingerprint=str(manifest.calculation.get("settings_fingerprint", "")),
        manifest_fingerprint=manifest.batch_fingerprint,
        command=command,
        config=config,
        ledger_path=record_path,
        conda_env=SQS_CONDA_ENV,
        module_name=SQS_MODULE,
        work_dir=PROJECT_ROOT,
        run_id=effective_run_id,
        retry=rerun,
        dependency_job_ids=dependency_job_ids,
        artifact_checker=artifact_complete,
    )
    if not result.ok:
        log.error(
            "Refused to submit %s for %s: state=%s intent=%s job_id=%s",
            stage,
            name,
            result.state,
            result.intent_key,
            result.job_id,
        )
        return 1
    print(
        json.dumps(
            {
                "name": name,
                "stage": stage,
                "job_id": result.job_id,
                "intent_key": result.intent_key,
                "state": result.state,
                "submitted": result.submitted,
                "root": str(root),
            },
            indent=2,
        )
    )
    return 0


def _endpoint_dependency_job_ids(
    manager: Any,
    record_path: Path,
    *,
    name: str,
    manifest_fingerprint: str,
) -> list[str]:
    """Return completed endpoint job IDs or refuse an unsafe points submit."""

    # Reconcile first so the subsequent dependency set is derived from the
    # persisted state, not a stale in-memory snapshot.
    manager.reconcile_ledger(record_path)
    entries = manager.ledger_entries(record_path)
    endpoint_entries = [
        entry
        for entry in entries
        if (
            entry.get("system_name") == name
            and entry.get("stage") == "endpoints"
            and entry.get("manifest_fingerprint") == manifest_fingerprint
        )
    ]
    if any(entry.get("legacy") for entry in entries if entry.get("stage") == "endpoints"):
        raise SqsOrchestrationError(
            "Legacy endpoint ledger cannot be used for points; reconcile it explicitly"
        )
    # A failed job followed by an explicit retry leaves multiple ledger
    # entries for the same endpoint intent.  Only the latest attempt is the
    # current dependency; historical attempts must not block a valid retry.
    latest_by_intent: dict[tuple[str, str], dict[str, Any]] = {}
    for entry in endpoint_entries:
        key = (str(entry.get("system_name")), str(entry.get("stage")))
        previous = latest_by_intent.get(key)
        if previous is None or int(entry.get("attempt", 1)) > int(previous.get("attempt", 1)):
            latest_by_intent[key] = entry

    dependencies: list[str] = []
    for entry in latest_by_intent.values():
        if not entry.get("job_id"):
            raise SqsOrchestrationError(
                "Endpoint ledger has no job ID; refusing points submission"
            )
        state, _ = manager.reconcile_entry(entry)
        if state != "completed":
            raise SqsOrchestrationError(
                "Endpoint ledger is not successfully reconciled; refusing points submission"
            )
        dependencies.append(str(entry["job_id"]))
    return dependencies


def _action_from_args(args: argparse.Namespace) -> str | None:
    """Resolve the positional command and legacy boolean action flags."""

    flag_actions = [action for action in PORTAL_ACTIONS if getattr(args, action, False)]
    command = getattr(args, "command", None)
    if command:
        conflicting = [action for action in flag_actions if action != command]
        if conflicting:
            joined = ", ".join(conflicting)
            raise ValueError(f"command {command!r} conflicts with --{joined}")
        return command
    return flag_actions[0] if flag_actions else None


def run_portal(args: argparse.Namespace) -> int:
    """Dispatch the compact main_*-style interface."""

    names = _selected_names(args)
    action = _action_from_args(args)
    if action is None:
        raise ValueError("provide a SQS command or legacy action flag")
    if args.slurm and action not in {"endpoints", "points"}:
        raise ValueError("--slurm is only valid with --endpoints or --points")
    exit_code = 0
    for name in names:
        try:
            if action == "generate":
                result = _portal_generate(
                    name,
                    dry_run=args.dry_run,
                    rerun=args.rerun,
                    run_id=args.run_id,
                )
            elif action in {"endpoints", "points"}:
                if args.slurm:
                    result = _portal_submit(
                        name,
                        action,
                        rerun=args.rerun,
                        dry_run=args.dry_run,
                        run_id=args.run_id,
                    )
                else:
                    result = _portal_execute(
                        name,
                        action,
                        dry_run=args.dry_run,
                        rerun=args.rerun,
                        run_id=args.run_id,
                    )
            elif action == "collect":
                result = _portal_collect(name, dry_run=args.dry_run, run_id=args.run_id)
            elif action == "reconcile":
                if args.dry_run:
                    print(
                        json.dumps(
                            {"name": name, "action": "reconcile", "write": False},
                            separators=(",", ":"),
                        )
                    )
                    result = 0
                else:
                    result = _portal_reconcile(name, run_id=args.run_id)
            else:
                result = _portal_status(
                    name, validate=action == "validate", run_id=args.run_id
                )
            exit_code = max(exit_code, int(result or 0))
        except (ValueError, FileNotFoundError, SqsOrchestrationError) as exc:
            log.error("%s: %s", name, exc)
            exit_code = 1
    return exit_code


def _add_portal_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "command",
        nargs="?",
        choices=PORTAL_ACTIONS,
        metavar="COMMAND",
        help="SQS stage (defaults to the legacy --* action flags)",
    )
    parser.add_argument("--name", "-n", help="Run one configured SQS structure")
    parser.add_argument("--batch", "-b", action="store_true", help="Run all STRUCTURE_NAMES")
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument("--generate", action="store_true", help="Generate SQS structures")
    actions.add_argument("--endpoints", action="store_true", help="Run pure-element endpoints")
    actions.add_argument("--points", action="store_true", help="Run the three SQS composition points")
    actions.add_argument("--collect", action="store_true", help="Collect completed points into the SQS JSON")
    actions.add_argument("--status", "--check", dest="status", action="store_true", help="Show batch status")
    actions.add_argument("--validate", action="store_true", help="Validate batch manifests")
    parser.add_argument("--slurm", "-s", action="store_true", help="Submit endpoint/point jobs to Slurm")
    parser.add_argument("--rerun", "-r", action="store_true", help="Rebuild or rerun existing work")
    parser.add_argument("--run-id", help="Operate on an explicit immutable SQS run")
    parser.add_argument("--dry-run", action="store_true", help="Print actions without running VASP or submitting Slurm")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="HT-VASP SQS mixing workflow",
        epilog=(
            "Examples: main_sqs.py generate; "
            "main_sqs.py endpoints --slurm; "
            "main_sqs.py points --slurm; "
            "main_sqs.py status"
        ),
    )
    _add_portal_arguments(parser)
    return parser


def main(argv: list[str] | None = None) -> int:
    try:
        args = build_parser().parse_args(argv)
        if _action_from_args(args) is not None:
            return int(run_portal(args))
        build_parser().print_help()
        return 0
    except (ValueError, FileNotFoundError, SqsOrchestrationError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

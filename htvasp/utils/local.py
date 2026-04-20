"""HT-VASP local job runner with resume support."""

from __future__ import annotations

import logging
import shutil
import traceback
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

from jobflow import SETTINGS, Flow, Job, JobStore, Response, initialize_logger
from jobflow.core.flow import get_flow
from jobflow.core.reference import OnMissing, OutputReference
from monty.os import cd

if TYPE_CHECKING:
    import jobflow

log = logging.getLogger(__name__)


def _get_job_dir_name(job: Job, prefix: int, counter: dict[str, int]) -> str:
    """Generate directory name for a job.

    Args:
        job: The job to generate directory name for.
        prefix: Global sequential prefix number.
        counter: Dictionary tracking job name occurrences.

    Returns:
        Directory name like "1-relax_1" or "2-relax_1-2".
    """
    name = job.name.replace(" ", "_")

    if name not in counter:
        counter[name] = 1
        return f"{prefix}-{name}"

    counter[name] += 1
    return f"{prefix}-{name}-{counter[name]}"


def _compute_job_dir_names(flow: Flow) -> dict[str, str]:
    """Compute directory names for all jobs in a flow.

    Args:
        flow: The flow to process.

    Returns:
        Mapping from job name to directory name.
    """
    dir_names: dict[str, str] = {}
    prefix = 1
    counter: dict[str, int] = {}

    for job, _ in flow.iterflow():
        dir_names[job.name] = _get_job_dir_name(job, prefix, counter)
        prefix += 1

    return dir_names


def run_locally_custom(
    flow: Flow | Job | list[Job],
    log_fmt: bool | str = "%(asctime)s[%(levelname)s]%(message)s",
    store: JobStore | None = None,
    root_dir: str | Path | None = None,
    ensure_success: bool = False,
    allow_external_references: bool = False,
    raise_immediately: bool = False,
    resume: bool = True,
) -> dict[str, dict[int, Response]]:
    """Run a Job or Flow locally with sequential directory naming and resume support.

    Directory naming: {prefix}-{name} or {prefix}-{name}-{suffix}
    - prefix: global sequential number (1, 2, 3...)
    - name: job name with spaces replaced by underscores
    - suffix: counter for duplicate names (omitted for first occurrence)

    Args:
        flow: A job or flow.
        log_fmt: Logging format. Defaults to standard format.
        store: A job store. If not specified, uses default JobStore.
        root_dir: The root directory to run the jobs in.
        ensure_success: Raise an error if the flow was not executed successfully.
        allow_external_references: If False, all references must be from jobs in the same flow.
        raise_immediately: If True, raise an exception immediately if a job fails.
        resume: Whether to resume from previously completed jobs by name matching.

    Returns:
        The responses of the jobs, keyed by job uuid and index.
    """
    if store is None:
        store = SETTINGS.JOB_STORE

    root_dir = Path(root_dir).resolve() if root_dir else Path.cwd()
    root_dir.mkdir(exist_ok=True, parents=True)
    store.connect()

    if log_fmt:
        initialize_logger(fmt=log_fmt if isinstance(log_fmt, str) else "")

    flow = get_flow(flow, allow_external_references=allow_external_references)

    # Resume preparation
    completed_uuids = _prepare_resume(flow, store, root_dir) if resume else set()

    # Execution state
    stopped_parents: set[str] = set()
    errored: set[str] = set()
    responses: dict[str, dict[int, jobflow.Response]] = defaultdict(dict)
    stop_jobflow = False
    prefix = 1
    counter: dict[str, int] = {}

    for job, parents in flow.iterflow():
        job_dir = root_dir / _get_job_dir_name(job, prefix, counter)
        prefix += 1

        # For incomplete jobs, clean the directory before running
        if job.uuid not in completed_uuids and job_dir.exists():
            shutil.rmtree(job_dir)
            log.debug(f"Cleaned directory for incomplete job: {job.name}")

        job_dir.mkdir(parents=True, exist_ok=True)

        with cd(job_dir):
            response, stop_jobflow = _run_job(
                job, parents, store, completed_uuids,
                stopped_parents, errored, raise_immediately, responses
            )

        if response is not None:
            response.job_dir = job_dir

        if stop_jobflow:
            break

    log.info("Finished executing jobs locally")

    if ensure_success and stop_jobflow:
        raise RuntimeError("Flow did not finish running successfully")

    return dict(responses)


def _run_job(
    job: Job,
    parents: list[str],
    store: JobStore,
    completed_uuids: set[str],
    stopped_parents: set[str],
    errored: set[str],
    raise_immediately: bool,
    responses: dict[str, dict[int, jobflow.Response]],
) -> tuple[Response | None, bool]:
    """Run a single job with all checks and error handling.

    Args:
        job: The job to run.
        parents: Parent job uuids.
        store: The job store.
        completed_uuids: Set of completed job uuids.
        stopped_parents: Set of stopped parent uuids.
        errored: Set of errored job uuids.
        raise_immediately: Whether to raise on error.
        responses: Response accumulator.

    Returns:
        Tuple of (response, should_stop_jobflow).
    """
    log.info(f"Starting job - {job.name} ({job.uuid})")

    # Check if parent stopped
    if set(parents).intersection(stopped_parents):
        log.info(f"{job.name} is a child of a job with stop_children=True, skipping...")
        stopped_parents.add(job.uuid)
        return None, False

    # Check if parent errored
    if set(parents).intersection(errored):
        if job.config.on_missing_references == OnMissing.ERROR:
            errored.add(job.uuid)
            return None, False

    # Resume: skip completed jobs
    if job.uuid in completed_uuids:
        if set(parents).issubset(completed_uuids):
            response = _load_completed_output(job, store)
            if response is not None:
                responses[job.uuid][job.index] = response
                log.info(f"Skipped completed job: {job.name}")
                return response, False

        log.warning(f"Cannot skip {job.name}: dependencies incomplete, will re-execute")
        completed_uuids.discard(job.uuid)

    # Execute job
    if raise_immediately:
        response = job.run(store=store)
    else:
        response = _execute_job_safe(job, store, errored)

    if response is None:
        return None, False

    responses[job.uuid][job.index] = response

    if response.stored_data is not None:
        log.warning("Response.stored_data is not supported with local manager.")

    if response.stop_children:
        stopped_parents.add(job.uuid)

    if response.stop_jobflow:
        return None, True

    # Handle diversions
    for diversion in (response.replace, response.detour, response.addition):
        if diversion is None:
            continue
        # Note: Diversion handling would require recursive call, simplified here
        log.warning(f"Diversion detected for {job.name}, but not fully supported in resume mode")

    return response, False


def _prepare_resume(flow: Flow, store: JobStore, root_dir: Path) -> set[str]:
    """Restore uuids for existing jobs and return set of completed uuids.

    Handles the case where Maker.make() generates new uuids on each run.
    By restoring the original uuids, output references can be resolved correctly.
    Also updates dir_name in store to point to current flow directory.

    Args:
        flow: The flow to prepare for resume.
        store: The job store with previous results.
        root_dir: Current flow directory.

    Returns:
        Set of completed job uuids.
    """
    # Compute directory names for all jobs
    job_dir_names = _compute_job_dir_names(flow)

    # Collect uuid mappings (new_uuid -> old_uuid)
    uuid_mapping: dict[str, str] = {}
    for job, _ in flow.iterflow():
        existing = store.query_one({"name": job.name, "index": job.index})
        if existing is not None:
            uuid_mapping[job.uuid] = existing["uuid"]

    if not uuid_mapping:
        return set()

    # Update dir_name in store for completed jobs
    _update_dir_names(store, uuid_mapping, root_dir, job_dir_names)

    # Update OutputReferences in job arguments BEFORE set_uuid
    _update_job_references(flow, uuid_mapping)

    # Restore uuids for jobs found in store
    for job, _ in flow.iterflow():
        if job.uuid in uuid_mapping:
            old_uuid = uuid_mapping[job.uuid]
            job.set_uuid(old_uuid)
            log.debug(f"Restored uuid for {job.name}: {job.uuid[:8]}")

    # Validate outputs and collect completed uuids
    completed_uuids: set[str] = set()
    for job, _ in flow.iterflow():
        existing = store.query_one({"name": job.name, "index": job.index})
        if existing is None:
            continue

        if not _validate_job_output(job, store):
            log.warning(f"Job {job.name} has incomplete output, will re-execute")
            continue

        completed_uuids.add(job.uuid)
        log.debug(f"Verified completed job: {job.name}")

    if completed_uuids:
        log.info(f"Found {len(completed_uuids)} previously completed jobs")

    return completed_uuids


def _update_dir_names(
    store: JobStore,
    uuid_mapping: dict[str, str],
    root_dir: Path,
    job_dir_names: dict[str, str],
) -> None:
    """Update dir_name in store to point to current flow directory.

    When resuming, the dir_name stored in the database points to an old
    temporary directory that no longer exists. This function updates the
    dir_name to point to the current flow directory.

    Args:
        store: The job store.
        uuid_mapping: Mapping from new_uuid to old_uuid.
        root_dir: Current flow directory.
        job_dir_names: Mapping from job name to directory name.
    """
    for old_uuid in uuid_mapping.values():
        doc = store.query_one({"uuid": old_uuid})
        if doc is None:
            continue

        job_name = doc.get("name")
        if job_name not in job_dir_names:
            continue

        output = doc.get("output")
        if not isinstance(output, dict):
            continue

        if "dir_name" not in output:
            continue

        new_dir_path = str(root_dir / job_dir_names[job_name])
        output["dir_name"] = new_dir_path
        store.update(doc, key="uuid")
        log.debug(f"Updated dir_name for {job_name}: {new_dir_path}")


def _update_job_references(flow: Flow, uuid_mapping: dict[str, str]) -> None:
    """Update OutputReferences in job arguments to use restored uuids.

    When Maker.make() is called, it creates new jobs with new uuids. When we
    restore old uuids via set_uuid(), the OutputReferences in other jobs'
    arguments still point to the new uuids. This function updates those
    references to point to the restored (old) uuids.

    Args:
        flow: The flow containing jobs to update.
        uuid_mapping: Mapping from new_uuid to old_uuid.
    """
    def update_reference(obj: object) -> None:
        """Recursively update OutputReferences in an object."""
        if isinstance(obj, OutputReference):
            if obj.uuid in uuid_mapping:
                new_uuid = obj.uuid
                old_uuid = uuid_mapping[new_uuid]
                obj.set_uuid(old_uuid, inplace=True)
                log.debug(f"Updated reference: {new_uuid[:8]} -> {old_uuid[:8]}")
            return

        if isinstance(obj, dict):
            for v in obj.values():
                update_reference(v)
            return

        if isinstance(obj, (list, tuple)):
            for item in obj:
                update_reference(item)

    for job, _ in flow.iterflow():
        for arg in job.function_args:
            update_reference(arg)
        for kwarg in job.function_kwargs.values():
            update_reference(kwarg)


def _validate_job_output(job: Job, store: JobStore) -> bool:
    """Validate that a job's output exists and can be loaded.

    Args:
        job: The job to validate.
        store: The job store.

    Returns:
        True if output exists and can be loaded, False otherwise.
    """
    try:
        output = store.get_output(uuid=job.uuid, which=job.index, load=True)
        return output is not None
    except Exception:
        return False


def _load_completed_output(job: Job, store: JobStore) -> Response | None:
    """Load output from completed job in store.

    Args:
        job: The job to load output for.
        store: The job store.

    Returns:
        Response with loaded output, or None if loading failed.
    """
    try:
        output = store.get_output(uuid=job.uuid, which=job.index, load=True)
        return Response(output=output)
    except Exception as e:
        log.warning(f"Failed to load output for {job.name}: {e}")
        return None


def _execute_job_safe(job: Job, store: JobStore, errored: set[str]) -> Response | None:
    """Execute job with error handling.

    Args:
        job: The job to execute.
        store: The job store.
        errored: Set of errored job uuids (will be updated on error).

    Returns:
        Response on success, None on error.
    """
    try:
        return job.run(store=store)
    except Exception:
        log.info(f"{job.name} failed with exception:\n{traceback.format_exc()}")
        errored.add(job.uuid)
        return None

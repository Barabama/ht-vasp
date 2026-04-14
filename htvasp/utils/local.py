"""HT-VASP local job runner with resume support."""

from __future__ import annotations

import logging
import traceback
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

from jobflow import SETTINGS, Flow, Job, JobStore, Response, initialize_logger
from jobflow.core.flow import get_flow
from jobflow.core.reference import OnMissing
from monty.os import cd

if TYPE_CHECKING:
    import jobflow

logger = logging.getLogger(__name__)


def run_locally_custom(
    flow: Flow | Job | list[Job],
    log: bool | str = True,
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
        log: Controls logging. Defaults to True.
        store: A job store. If not specified, uses default JobStore.
        root_dir: The root directory to run the jobs in.
        ensure_success: Raise an error if the flow was not executed successfully.
        allow_external_references: If False, all references must be from jobs in the same flow.
        raise_immediately: If True, raise an exception immediately if a job fails.
        resume: Whether to resume from previously completed jobs by name matching.

    Returns:
        The responses of the jobs.
    """
    if store is None:
        store = SETTINGS.JOB_STORE

    root_dir = Path.cwd() if root_dir is None else Path(root_dir).resolve()
    root_dir.mkdir(exist_ok=True, parents=True)

    store.connect()

    if log:
        initialize_logger(fmt=log if isinstance(log, str) else "")

    flow = get_flow(flow, allow_external_references=allow_external_references)

    # Directory naming state
    global_prefix = 1
    name_counter: dict[str, int] = {}
    completed_uuids: set[str] = set()

    if resume:
        completed_uuids = _prepare_resume(flow, store)

    stopped_parents: set[str] = set()
    errored: set[str] = set()
    responses: dict[str, dict[int, jobflow.Response]] = defaultdict(dict)
    stop_jobflow = False

    def _get_job_dir(job: Job) -> Path:
        """Generate directory for a job with sequential naming."""
        nonlocal global_prefix
        name = job.name.replace(" ", "_")

        if name not in name_counter:
            name_counter[name] = 1
            dir_name = f"{global_prefix}-{name}"
            global_prefix += 1
        else:
            name_counter[name] += 1
            dir_name = f"{global_prefix - 1}-{name}-{name_counter[name]}"

        job_dir = root_dir / dir_name
        job_dir.mkdir(parents=True, exist_ok=True)
        return job_dir

    def _run_job(job: Job, parents: list[str]) -> tuple[Response | None, bool]:
        """Run a single job."""
        nonlocal stop_jobflow

        if stop_jobflow:
            return None, True

        # Check if parent stopped
        if set(parents).intersection(stopped_parents):
            logger.info(f"{job.name} is a child of a job with stop_children=True, skipping...")
            stopped_parents.add(job.uuid)
            return None, False

        # Check if parent errored
        if set(parents).intersection(errored) and job.config.on_missing_references == OnMissing.ERROR:
            errored.add(job.uuid)
            return None, False

        # Resume: skip completed jobs
        if job.uuid in completed_uuids:
            if set(parents).issubset(completed_uuids):
                response = _load_completed_output(job, store)
                if response is not None:
                    responses[job.uuid][job.index] = response
                    logger.info(f"Skipped completed job: {job.name}")
                    return response, False
            logger.warning(f"Cannot skip {job.name}: dependencies incomplete, will re-execute")
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
            logger.warning("Response.stored_data is not supported with local manager.")

        if response.stop_children:
            stopped_parents.add(job.uuid)

        if response.stop_jobflow:
            stop_jobflow = True
            return None, True

        # Handle diversions
        diversion_responses = []
        if response.replace:
            diversion_responses.append(_run(response.replace))
        if response.detour:
            diversion_responses.append(_run(response.detour))
        if response.addition:
            diversion_responses.append(_run(response.addition))

        if not all(diversion_responses):
            return None, False
        return response, False

    def _run(root_flow) -> bool:
        """Run a flow."""
        encountered_bad_response = False
        for job, parents in root_flow.iterflow():
            job_dir = _get_job_dir(job)
            with cd(job_dir):
                response, jobflow_stopped = _run_job(job, parents)

            if response is not None:
                response.job_dir = job_dir

            encountered_bad_response = encountered_bad_response or response is None
            if jobflow_stopped:
                return False

        return not encountered_bad_response

    logger.info("Started executing jobs locally")
    finished_successfully = _run(flow)
    logger.info("Finished executing jobs locally")

    if ensure_success and not finished_successfully:
        raise RuntimeError("Flow did not finish running successfully")

    return dict(responses)


def _prepare_resume(flow: Flow, store: JobStore) -> set[str]:
    """Restore uuids for existing jobs and return set of completed uuids.

    Handles the case where Maker.make() generates new uuids on each run.
    By restoring the original uuids, output references can be resolved correctly.
    """
    completed_uuids: set[str] = set()

    # First pass: restore uuids
    for job, _ in flow.iterflow():
        existing = store.query_one({"name": job.name, "index": job.index})
        if existing is not None:
            job.set_uuid(existing["uuid"])
            logger.debug(f"Found existing record: {job.name} (uuid: {job.uuid[:8]})")

    # Second pass: validate outputs
    for job, _ in flow.iterflow():
        existing = store.query_one({"name": job.name, "index": job.index})
        if existing is None:
            continue

        if not _validate_job_output(job, store):
            logger.warning(f"Job {job.name} has incomplete output, will re-execute")
            continue

        completed_uuids.add(job.uuid)
        logger.debug(f"Verified completed job: {job.name}")

    if completed_uuids:
        logger.info(f"Found {len(completed_uuids)} previously completed jobs")

    return completed_uuids


def _validate_job_output(job: Job, store: JobStore) -> bool:
    """Validate that a job's output exists and can be loaded."""
    try:
        output = store.get_output(uuid=job.uuid, which=job.index, load=True)
        return output is not None
    except Exception:
        return False


def _load_completed_output(job: Job, store: JobStore) -> Response | None:
    """Load output from completed job in store."""
    try:
        output = store.get_output(uuid=job.uuid, which=job.index, load=True)
        return Response(output=output)
    except Exception as e:
        logger.warning(f"Failed to load output for {job.name}: {e}")
        return None


def _execute_job_safe(job: Job, store: JobStore, errored: set[str]) -> Response | None:
    """Execute job with error handling."""
    try:
        return job.run(store=store)
    except Exception:
        logger.info(f"{job.name} failed with exception:\n{traceback.format_exc()}")
        errored.add(job.uuid)
        return None

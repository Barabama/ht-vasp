"""
HT-VASP - Custom run_locally implementation

Custom run_locally with flexible directory naming and job control.
"""

from __future__ import annotations

import typing
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone
from random import randint

from monty.os import cd

from jobflow import (
    Flow,
    Job,
    JobStore,
    initialize_logger,
    OnMissing,
    Response,
    SETTINGS,
)
from jobflow.core.flow import get_flow
from jobflow.managers.local import logger

if typing.TYPE_CHECKING:
    import jobflow


def _find_existing_job_dir(job: Job, root_dir: Path, resume: bool) -> Path | None:
    """Find an existing directory for a job if resume is enabled."""
    if not resume:
        return None

    for item in root_dir.iterdir():
        if not item.is_dir():
            continue
        if job.name.replace(" ", "_") in item.name:
            logger.debug(f"Found existing directory for job {job.name}: {item}")
            return item
    return None


def _check_job_completed(job: Job, store: JobStore, resume: bool) -> dict | None:
    """Check if a job has already been completed in the store."""
    if not resume:
        return None

    try:
        # First try by uuid and index (exact match)
        existing_doc = store.query_one(criteria={"uuid": job.uuid, "index": job.index})

        # If not found by uuid, try by name and index (more flexible)
        if existing_doc is None:
            existing_doc = store.query_one(criteria={"name": job.name, "index": job.index})

        if existing_doc is not None:
            logger.debug(f"Found completed job: {job.name} (index {job.index})")
        return existing_doc
    except Exception as e:
        logger.debug(f"Error checking job completion: {e}")
        return None


def _load_existing_output(existing_doc: dict, store: JobStore) -> Response:
    """Load existing output from store and create a Response."""
    try:
        existing_output = store.get_output(
            uuid=existing_doc["uuid"], which=existing_doc["index"], load=True
        )
        return Response(output=existing_output)
    except Exception as load_error:
        logger.warning(f"Could not load existing output, using dummy response: {load_error}")
        return Response(output=None)


def _generate_job_dir(
    job: Job,
    root_dir: Path,
    dir_format: str,
    dir_prefix: str | None,
    dir_suffix: str,
    auto_prefix_enabled: bool,
    current_prefix: int | None,
) -> tuple[Path, int | None]:
    """Generate a directory name for a job."""
    # Build directory name from format
    dir_name = dir_format.format(
        name=job.name.replace(" ", "_"),
        uuid=job.uuid[:8],
        index=job.index,
    )

    # Add prefix and suffix
    if auto_prefix_enabled:
        current_dir_prefix = f"{current_prefix}-"
        new_prefix = current_prefix + 1
    else:
        current_dir_prefix = dir_prefix
        new_prefix = current_prefix

    dir_name = f"{current_dir_prefix}{dir_name}{dir_suffix}"

    # Handle duplicates by adding a sequential suffix
    idx = 1
    job_dir = root_dir.joinpath(dir_name)
    while job_dir.exists():
        tmp_dir_name = f"{dir_name}_{idx}"
        job_dir = root_dir.joinpath(tmp_dir_name)
        idx += 1

    job_dir.mkdir(parents=True, exist_ok=True)
    return job_dir, new_prefix


def run_locally_custom(
    flow: Flow | Job | list[Job],
    log: bool | str = True,
    store: JobStore | None = None,
    root_dir: Path | str = "",
    ensure_success: bool = False,
    allow_external_references: bool = False,
    raise_immediately: bool = False,
    dir_format: str = "{name}",
    dir_prefix: str | None = None,
    dir_suffix: str = "",
    resume: bool = True,
) -> dict[str, dict[int, Response]]:
    """
    Run a Job or Flow locally with custom directory naming.

    This is an enhanced version of jobflow's run_locally that supports:
    - Custom directory naming patterns
    - Prefix/suffix for job directories
    - Resume from previously completed jobs
    - Better control over job execution

    Args:
        flow: A job or flow
        log: Controls logging. Defaults to True
        store: A job store. If not specified, uses default JobStore
        root_dir: The root directory to run the jobs in. If None, uses current directory
        ensure_success: Raise an error if the flow was not executed successfully
        allow_external_references: If False all references to other outputs
            should be from other Jobs of the same Flow
        raise_immediately: If True, raise an exception immediately if a job fails
        dir_format: Format string for directory names. Available placeholders:
            - {name}: job name
            - {uuid}: job uuid (first 8 chars)
            - {index}: job index
            Default is "{name}"
        dir_prefix: Prefix to add to directory names. If None, uses auto-numbering starting from 1. Defaults to None
        dir_suffix: Suffix to add to directory names
        resume: Whether to resume from previously completed jobs. Defaults to True

    Returns:
        dict[str, dict[int, Response]]: The responses of the jobs.

    Examples
    >>> # Use job name as directory name
    >>> run_locally_custom(job, store=store, root_dir="./runs", dir_format="{name}")

    >>> # Use job name with prefix and suffix
    >>> run_locally_custom(job, store=store, root_dir="./runs", dir_prefix="calc_", dir_suffix="_v1")

    >>> # Use uuid as directory name
    >>> run_locally_custom(job, store=store, root_dir="./runs", dir_format="{uuid}")

    >>> # Combined format
    >>> run_locally_custom(job, store=store, root_dir="./runs", dir_format="{name}_{index}")
    """
    if store is None:
        store = SETTINGS.JOB_STORE

    root_dir = Path(root_dir).resolve()
    root_dir.mkdir(exist_ok=True, parents=True)

    # Initialize auto-numbering if dir_prefix is None
    auto_prefix_enabled = False
    current_prefix = None
    if dir_prefix is None:
        auto_prefix_enabled = True
        existing_prefixes = []
        for item in root_dir.iterdir():
            if not item.is_dir():
                continue
            try:
                parts = item.name.split("-", 1)
                if parts[0].isdigit():
                    existing_prefixes.append(int(parts[0]))
            except (ValueError, IndexError):
                pass

        if existing_prefixes:
            current_prefix = max(existing_prefixes) + 1
        else:
            current_prefix = 1
    else:
        if dir_prefix and not dir_prefix.endswith("-"):
            dir_prefix += "-"

    store.connect()

    if log:
        initialize_logger(fmt=log if isinstance(log, str) else "")

    flow = get_flow(flow, allow_external_references=allow_external_references)

    stopped_parents: set[str] = set()
    errored: set[str] = set()
    responses: dict[str, dict[int, Response]] = defaultdict(dict)
    stop_jobflow = False

    def _get_job_dir(job: Job) -> Path:
        """Generate or find a directory for a job."""
        nonlocal current_prefix

        existing_dir = _find_existing_job_dir(job, root_dir, resume)
        if existing_dir is not None:
            existing_dir.mkdir(parents=True, exist_ok=True)
            return existing_dir

        job_dir, new_prefix = _generate_job_dir(
            job,
            root_dir,
            dir_format,
            dir_prefix,
            dir_suffix,
            auto_prefix_enabled,
            current_prefix,
        )
        current_prefix = new_prefix
        return job_dir

    def _run_job(job: Job, parents) -> tuple[Response | None, bool]:
        """Run a single job."""
        nonlocal stop_jobflow

        if stop_jobflow:
            return None, True

        if len(set(parents).intersection(stopped_parents)) > 0:
            logger.info(f"{job.name} is a child of a job with stop_children=True, skipping...")
            stopped_parents.add(job.uuid)
            return None, False

        if (
            len(set(parents).intersection(errored)) > 0
            and job.config.on_missing_references == OnMissing.ERROR
        ):
            errored.add(job.uuid)
            return None, False

        # Check if job is already completed in store
        existing_doc = _check_job_completed(job, store, resume)
        if existing_doc is not None:
            logger.info(f"Job {job.name} (index {job.index}) already completed, skipping...")
            existing_response = _load_existing_output(existing_doc, store)
            responses[job.uuid][job.index] = existing_response
            return existing_response, False

        # Run the job
        if raise_immediately:
            response = job.run(store=store)
        else:
            try:
                response = job.run(store=store)
            except Exception:
                import traceback

                logger.info(f"{job.name} failed with exception:\n{traceback.format_exc()}")
                errored.add(job.uuid)
                return None, False

        responses[job.uuid][job.index] = response

        if response.stored_data is not None:
            logger.warning("Response.stored_data is not supported with local manager.")

        if response.stop_children:
            stopped_parents.add(job.uuid)

        if response.stop_jobflow:
            stop_jobflow = True
            return None, True

        diversion_responses = []
        if response.replace is not None:
            diversion_responses.append(_run(response.replace))

        if response.detour is not None:
            diversion_responses.append(_run(response.detour))

        if response.addition is not None:
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

"""HT-VASP local job runner - simplified run_locally implementation."""

from __future__ import annotations

import logging
import typing
from collections import defaultdict
from pathlib import Path

from monty.os import cd

from jobflow import (
    Flow,
    Job,
    JobStore,
    Response,
    SETTINGS,
    initialize_logger,
)
from jobflow.core.flow import get_flow
from jobflow.core.reference import OnMissing

if typing.TYPE_CHECKING:
    import jobflow

logger = logging.getLogger(__name__)


def run_locally(
    flow: Flow | Job | list[Job],
    log: bool | str = True,
    store: JobStore | None = None,
    root_dir: str | Path | None = None,
    ensure_success: bool = False,
    allow_external_references: bool = False,
    raise_immediately: bool = False,
) -> dict[str, dict[int, Response]]:
    """
    Run a Job or Flow locally with sequential directory naming.

    Directory naming format: {prefix}-{name} or {prefix}-{name}-{suffix}
    - prefix: global sequential number (1, 2, 3...)
    - name: job name with spaces replaced by underscores
    - suffix: counter for duplicate names (omitted for first occurrence)

    Args:
        flow: A job or flow.
        log: Controls logging. Defaults to True.
        store: A job store. If not specified, uses default JobStore.
        root_dir: The root directory to run the jobs in. Defaults to current directory.
        ensure_success: Raise an error if the flow was not executed successfully.
        allow_external_references: If False, all references must be from jobs in the same flow.
        raise_immediately: If True, raise an exception immediately if a job fails.

    Returns:
        dict[str, dict[int, Response]]: The responses of the jobs.

    Examples:
        >>> # Jobs create directories: 1-job_a, 2-job_b, 3-job_a-1
        >>> run_locally(flow, root_dir="./runs")
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

    stopped_parents: set[str] = set()
    errored: set[str] = set()
    responses: dict[str, dict[int, jobflow.Response]] = defaultdict(dict)
    stop_jobflow = False

    def _get_job_dir(job: Job) -> Path:
        """Generate directory for a job with sequential naming."""
        nonlocal global_prefix
        name = job.name.replace(" ", "_")

        if name not in name_counter:
            # New name: assign new prefix
            name_counter[name] = 1
            dir_name = f"{global_prefix}-{name}"
            global_prefix += 1
        else:
            # Existing name: increment suffix counter
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

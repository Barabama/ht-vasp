"""Custom run_locally with custom directory naming."""

from __future__ import annotations

import typing
from pathlib import Path

if typing.TYPE_CHECKING:
    import jobflow


def run_locally_custom(
    flow: jobflow.Flow | jobflow.Job | list[jobflow.Job],
    log: bool | str = True,
    store: jobflow.JobStore | None = None,
    root_dir: str | Path | None = None,
    ensure_success: bool = False,
    allow_external_references: bool = False,
    raise_immediately: bool = False,
    dir_format: str = "{name}",
    dir_prefix: str = "",
    dir_suffix: str = "",
) -> dict[str, dict[int, jobflow.Response]]:
    """
    Run a Job or Flow locally with custom directory naming.

    Parameters
    ----------
    flow : Flow | Job | list[Job]
        A job or flow.
    log : bool | str
        Controls logging. Defaults to True.
    store : JobStore
        A job store. If not specified, uses default JobStore.
    root_dir : str | Path | None
        The root directory to run the jobs in. If None, uses current directory.
    ensure_success : bool
        Raise an error if the flow was not executed successfully.
    allow_external_references : bool
        If False all references to other outputs should be from other Jobs
        of the same Flow.
    raise_immediately : bool
        If True, raise an exception immediately if a job fails.
    dir_format : str
        Format string for directory names. Available placeholders:
        - {name}: job name
        - {uuid}: job uuid (first 8 chars)
        - {index}: job index
        Default is "{name}".
    dir_prefix : str
        Prefix to add to directory names.
    dir_suffix : str
        Suffix to add to directory names.

    Returns
    -------
    dict[str, dict[int, Response]]
        The responses of the jobs.

    Examples
    --------
    >>> # Use job name as directory name
    >>> run_locally_custom(job, root_dir="./runs", dir_format="{name}")

    >>> # Use job name with prefix and suffix
    >>> run_locally_custom(job, root_dir="./runs", dir_prefix="calc_", dir_suffix="_v1")

    >>> # Use uuid as directory name
    >>> run_locally_custom(job, root_dir="./runs", dir_format="{uuid}")

    >>> # Combined format
    >>> run_locally_custom(job, root_dir="./runs", dir_format="{name}_{index}")
    """
    from collections import defaultdict
    from datetime import datetime, timezone
    from random import randint

    from monty.os import cd

    from jobflow import SETTINGS, initialize_logger
    from jobflow.core.flow import get_flow
    from jobflow.core.reference import OnMissing
    from jobflow.managers.local import logger

    if store is None:
        store = SETTINGS.JOB_STORE

    root_dir = Path.cwd() if root_dir is None else Path(root_dir).resolve()
    root_dir.mkdir(exist_ok=True, parents=True)

    store.connect()

    if log:
        initialize_logger(fmt=log if isinstance(log, str) else "")

    flow = get_flow(flow, allow_external_references=allow_external_references)

    stopped_parents: set[str] = set()
    errored: set[str] = set()
    responses: dict[str, dict[int, jobflow.Response]] = defaultdict(dict)
    stop_jobflow = False

    def _get_job_dir(job):
        """Generate custom directory name for a job."""
        # Build directory name from format
        dir_name = dir_format.format(
            name=job.name,
            uuid=job.uuid[:8],
            index=job.index,
        )
        # Add prefix and suffix
        dir_name = f"{dir_prefix}{dir_name}{dir_suffix}"
        # Handle duplicates by adding a random suffix
        job_dir = root_dir / dir_name
        if job_dir.exists():
            dir_name = f"{dir_name}_{randint(1000, 9999)}"
            job_dir = root_dir / dir_name
        job_dir.mkdir(parents=True, exist_ok=True)
        return job_dir

    def _run_job(job: jobflow.Job, parents):
        nonlocal stop_jobflow

        if stop_jobflow:
            return None, True

        if len(set(parents).intersection(stopped_parents)) > 0:
            logger.info(
                f"{job.name} is a child of a job with stop_children=True, skipping..."
            )
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

                logger.info(
                    f"{job.name} failed with exception:\n{traceback.format_exc()}"
                )
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

        return response, False

    def _run(root_flow):
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

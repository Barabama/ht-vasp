"""HT-VASP CLI Runner - Concurrent job submission and monitoring."""

import time
import logging
from pathlib import Path
from threading import Thread, Semaphore
from typing import Any, Callable

from htvasp.slurm import SlurmJobManager

log = logging.getLogger(__name__)

# Default settings
DEFAULT_MAX_CONCURRENT = 2
DEFAULT_POLL_INTERVAL = 60  # seconds


def submit_and_monitor(
    name: str,
    manager: SlurmJobManager,
    semaphore: Semaphore,
    force: bool,
    results: dict[str, Any],
    command_template: Callable[[str, bool], str],
    config_kwargs: dict,
    workdir: str | Path = ".",
    poll_interval: int = DEFAULT_POLL_INTERVAL,
) -> None:
    """Submit job and monitor until completion.

    Args:
        name: Structure name
        manager: SlurmJobManager instance
        semaphore: Semaphore for concurrency control
        force: Force run flag
        results: Shared dict to store results
        command_template: Function(name, force) -> command string
        config_kwargs: Slurm config kwargs (job_name, output_log, ntasks, etc.)
        workdir: Working directory
        poll_interval: Polling interval in seconds
    """
    job_id = None
    try:
        # Build config
        config = manager.get_cpu_config(**config_kwargs)

        # Submit job
        command = command_template(name, force)
        job_id = manager.submit_command(
            command=command,
            config=config,
            workdir=workdir,
        )

        if not job_id:
            log.error(f"Failed to submit job for {name}")
            results[name] = {"status": "failed", "error": "submit failed"}
            return

        log.info(f"Submitted job for {name} with ID: {job_id}")

        # Monitor until completion
        while not manager.is_job_completed(job_id):
            time.sleep(poll_interval)

        log.info(f"Job {job_id} for {name} completed")
        results[name] = {"status": "completed", "job_id": job_id}

    except Exception as e:
        log.error(f"Error monitoring job for {name}: {e}")
        results[name] = {"status": "failed", "error": str(e), "job_id": job_id}
    finally:
        semaphore.release()


def submit_jobs(
    struct_names: list[str],
    command_template: Callable[[str, bool], str],
    manager: SlurmJobManager | None = None,
    force: bool = False,
    max_concurrent: int = DEFAULT_MAX_CONCURRENT,
    poll_interval: int = DEFAULT_POLL_INTERVAL,
    workdir: str | Path = ".",
    **config_kwargs,
) -> dict[str, Any]:
    """Submit jobs with concurrency limit.

    Args:
        struct_names: List of structure names
        command_template: Function(name, force) -> command string
        manager: SlurmJobManager instance (created if None)
        force: Force run flag
        max_concurrent: Maximum concurrent jobs
        poll_interval: Polling interval in seconds
        workdir: Working directory
        **config_kwargs: Slurm config kwargs (job_name, output_log, ntasks, etc.)
                         Can be values or callables that take name as argument.

    Returns:
        Dict mapping name -> {status, job_id, error?}
    """
    manager = manager or SlurmJobManager()
    semaphore = Semaphore(max_concurrent)
    threads = []
    results: dict[str, Any] = {}

    log.info(f"Starting job submission with max {max_concurrent} concurrent jobs")

    for name in struct_names:
        semaphore.acquire()  # Block if max jobs running

        # Build config_kwargs for this name (support callables)
        name_config = {}
        for key, value in config_kwargs.items():
            if callable(value):
                name_config[key] = value(name)
            else:
                name_config[key] = value

        t = Thread(
            target=submit_and_monitor,
            args=(name, manager, semaphore, force, results, command_template,
                  name_config, Path(workdir), poll_interval),
        )
        t.start()
        threads.append(t)

    # Wait for all threads to complete
    for t in threads:
        t.join()

    # Summary
    completed = sum(1 for r in results.values() if r.get("status") == "completed")
    failed = sum(1 for r in results.values() if r.get("status") == "failed")
    log.info(f"All jobs finished: {completed} completed, {failed} failed")

    return results


def run_batch(
    struct_names: list[str],
    run_func: Callable[[str, bool], None],
    force: bool = False,
) -> None:
    """Run all structures locally.

    Args:
        struct_names: List of structure names
        run_func: Function(name, force) to run a single structure
        force: Force run flag
    """
    for name in struct_names:
        run_func(name, force)

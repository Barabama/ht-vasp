import json
import shutil
import logging
import argparse
import time
from pathlib import Path
from threading import Thread, Semaphore
from typing import Any

from htvasp.model import Endmember
from htvasp.workflows import OJWorker
from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

# Concurrency settings
MAX_CONCURRENT_JOBS = 2
POLL_INTERVAL = 60  # seconds


CLUSTER_VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
}

# oj_incar to override INCAR settings for OJ calculations
global_incar = {
    "ENCUT": 400,
    "ISTART": 0,
    "ICHARG": 2,
    # Electronic
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ALGO": "Fast",
    "NELM": 100,
    "NELMIN": 6,
    "NELMDL": -6,
    # Ionic
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 50,
    "POTIM": 0.2,
    "EDIFF": 1e-5,
    "EDIFFG": -0.05,
    # Magnetic
    "ISPIN": 2,
    # Precision
    "ISYM": 0,
    "LREAL": "Auto",
    "PREC": "Normal",
    "SYMPREC": 1e-5,
    # Output
    "LWAVE": False,
    "LCHARG": False,
    "LORBIT": 11,
    "GGA": "PE",
    "KPAR": 2,
    "NCORE": 4,
    "GGA": "PE",
    "AMIX": 0.2,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.8,
    "BMIX_MAG": 1e-4,
}
struct_names = [
    "SER-Co",
    "SER-Fe",
    "SER-Mn",
    "SER-Ni",
    # "BCC-Co-Co",
    # "BCC-Co-Fe",
    # "BCC-Co-Mn",
    # "BCC-Co-Ni",
    # "BCC-Fe-Fe",
    # "BCC-Fe-Mn",
    # "BCC-Fe-Ni",
    # "BCC-Mn-Mn",
    # "BCC-Mn-Ni",
    # "BCC-Ni-Ni",
    "FCC-Co-Co",
    # "FCC-Co-Fe",
    # "FCC-Co-Mn",
    # "FCC-Co-Ni",
    # "FCC-Fe-Co",
    # "FCC-Fe-Fe",
    # "FCC-Fe-Mn",
    # "FCC-Fe-Ni",
    # "FCC-Mn-Co",
    # "FCC-Mn-Fe",
    # "FCC-Mn-Mn",
    # "FCC-Mn-Ni",
    # "FCC-Ni-Co",
    # "FCC-Ni-Fe",
    # "FCC-Ni-Mn",
    # "FCC-Ni-Ni",
]


def run_tick(name: str, force=False):
    """Run a single structure locally."""
    endmember = Endmember()
    flow_name = f"{name}-oj"
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("data/endmembers") / name / "ojflow"
    poscar_dir = Path("data/poscars")
    json_path = store_dir / f"{flow_name}.json"

    # Skip if already done (unless force is True)
    if not force and json_path.exists():
        with open(json_path, "r", encoding="utf-8") as jf:
            result = json.load(jf)
        if result.get("state", "failed") == "successful":
            log.info(f"Structure {name} already done")
            return

    # Run OJWorker
    log.info(f"Structure {name} start")

    if force and store_dir.exists():
        shutil.rmtree(store_dir)

    struct = endmember.get_poscar(name, poscar_dir)

    try:
        worker = OJWorker(
            vasp_args=CLUSTER_VASP_ARGS,
            global_incar=global_incar,
            j_count=4,
            extend_poscar=(2, 2, 2),
        )

        worker.run_flow(
            name,
            struct,
            flow_dir,
            store_dir,
            resume=not force,
        )
        output = worker.get_result()
        worker.write_result(data=output, json_path=json_path)

    except Exception as e:
        log.error(f"Structure {name} failed: {e}")
        return

    log.info(f"Structure {name} done")



def run_batch(force=False):
    """Run all structures locally."""
    for name in struct_names:
        run_tick(name, force=force)


def submit_and_monitor(
    name: str,
    manager: SlurmJobManager,
    semaphore: Semaphore,
    force: bool,
    results: dict[str, Any],
):
    """Submit job and monitor until completion.

    Args:
        name: Structure name
        manager: SlurmJobManager instance
        semaphore: Semaphore for concurrency control
        force: Force run flag
        results: Shared dict to store results
    """
    job_id = None
    try:
        config = manager.get_cpu_config(
            job_name=f"{name}-oj",
            output_log=f"logs/{name}-oj.log",
            ntasks=48,
        )
        job_id = manager.submit_command(
            command=f"python {__file__} --tick {name} {'--force' if force else ''}",
            config=config,
            conda_env="htvasp",
            workdir=".",
        )

        if not job_id:
            log.error(f"Failed to submit job for {name}")
            results[name] = {"status": "failed", "error": "submit failed"}
            return

        log.info(f"Submitted job for {name} with ID: {job_id}")

        # Monitor until completion
        while not manager.is_job_completed(job_id):
            time.sleep(POLL_INTERVAL)

        log.info(f"Job {job_id} for {name} completed")
        results[name] = {"status": "completed", "job_id": job_id}

    except Exception as e:
        log.error(f"Error monitoring job for {name}: {e}")
        results[name] = {"status": "failed", "error": str(e), "job_id": job_id}
    finally:
        semaphore.release()


def submit_jobs(force=False):
    """Submit jobs with concurrency limit.

    Limits concurrent Slurm jobs to MAX_CONCURRENT_JOBS (default: 2).
    Blocks when queue is full, resumes when a job completes.
    """
    manager = SlurmJobManager()
    semaphore = Semaphore(MAX_CONCURRENT_JOBS)
    threads = []
    results: dict[str, Any] = {}

    log.info(f"Starting job submission with max {MAX_CONCURRENT_JOBS} concurrent jobs")

    for name in struct_names:
        semaphore.acquire()  # Block if max jobs running

        t = Thread(
            target=submit_and_monitor,
            args=(name, manager, semaphore, force, results),
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run VASP task")
    parser.add_argument("--tick", type=str, help="run a single structure")
    parser.add_argument("--batch", action="store_true", help="run batch")
    parser.add_argument("--slurm", action="store_true", help="submit to slurm")
    parser.add_argument("--force", action="store_true", help="force run")
    args = parser.parse_args()

    if args.tick:
        run_tick(args.tick, force=args.force)
    elif args.batch:
        run_batch(args.force)
    elif args.slurm:
        submit_jobs(args.force)
    else:
        print("specify run mode: --batch or --slurm")

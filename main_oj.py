"""HT-VASP OJ Workflow - Magnetic exchange calculation."""

import json
import shutil
import logging
import argparse
from pathlib import Path

from htvasp.cli import submit_jobs, run_batch
from htvasp.model import Endmember
from htvasp.workflows import OJWorker

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

# VASP configuration
VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
}

GLOBAL_INCAR = {
    "ENCUT": 400,
    "ISTART": 0,
    "ICHARG": 2,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ALGO": "Fast",
    "NELM": 100,
    "NELMIN": 6,
    "NELMDL": -6,
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 50,
    "POTIM": 0.2,
    "EDIFF": 1e-5,
    "EDIFFG": -0.05,
    "ISPIN": 2,
    "ISYM": 0,
    "LREAL": "Auto",
    "PREC": "Normal",
    "SYMPREC": 1e-5,
    "LWAVE": False,
    "LCHARG": False,
    "LORBIT": 11,
    "GGA": "PE",
    "KPAR": 2,
    "NCORE": 4,
    "AMIX": 0.2,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.8,
    "BMIX_MAG": 1e-4,
}

# Structure names to process
STRUCT_NAMES = [
    "SER-Co",
    "SER-Fe",
    "SER-Mn",
    "SER-Ni",
    "FCC-Co-Co",
]


def run_tick(name: str, force: bool = False) -> None:
    """Run a single structure locally."""
    endmember = Endmember()
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("data/endmembers") / name / "ojflow"
    json_path = store_dir / f"{name}-oj.json"

    # Skip if already done
    if not force and json_path.exists():
        with open(json_path) as f:
            result = json.load(f)
        if result.get("state") == "successful":
            log.info(f"Structure {name} already done")
            return

    log.info(f"Structure {name} start")

    if force and store_dir.exists():
        shutil.rmtree(store_dir)

    struct = endmember.get_poscar(name, Path("data/poscars"))

    try:
        worker = OJWorker(
            vasp_args=VASP_ARGS,
            global_incar=GLOBAL_INCAR,
            j_count=4,
            extend_poscar=(2, 2, 2),
        )
        worker.run_flow(name, struct, flow_dir, store_dir, resume=not force)
        output = worker.get_result()
        worker.write_result(data=output, json_path=json_path)
    except Exception as e:
        log.error(f"Structure {name} failed: {e}")

    log.info(f"Structure {name} done")


def command_template(name: str, force: bool) -> str:
    """Generate command for Slurm submission."""
    return f"python {__file__} --tick {name} {'--force' if force else ''}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OJ workflow")
    parser.add_argument("--tick", type=str, help="Run single structure")
    parser.add_argument("--batch", action="store_true", help="Run all locally")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--force", action="store_true", help="Force re-run")
    args = parser.parse_args()

    if args.tick:
        run_tick(args.tick, force=args.force)
    elif args.batch:
        run_batch(STRUCT_NAMES, run_tick, force=args.force)
    elif args.slurm:
        submit_jobs(
            struct_names=STRUCT_NAMES,
            command_template=command_template,
            force=args.force,
            job_name=lambda n: f"{n}-oj",
            output_log=lambda n: f"logs/{n}-oj.log",
            ntasks=48,
            conda_env="htvasp",
        )
    else:
        parser.print_help()

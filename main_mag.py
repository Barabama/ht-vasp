"""HT-VASP mag Workflow - Magnetic exchange calculation."""

import json
import shutil
import logging
import argparse
from pathlib import Path

from htvasp.model import Endmember
from htvasp.workflows import StaticWorker
from htvasp.slurm import SlurmJobManager

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
    "EDIFF": 1e-6,
    "EDIFFG": -0.02,
    "ISPIN": 2,
    "ISYM": 0,
    "LREAL": "Auto",
    "PREC": "Normal",
    "SYMPREC": 1e-5,
    "LWAVE": False,
    "LCHARG": False,
    "LORBIT": 10,
    "GGA": "PE",
    "KPAR": 4,
    "NCORE": 2,
    "AMIX": 0.2,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.8,
    "BMIX_MAG": 1e-4,
}

# Structure names to process
STRUCT_NAMES = [
    "SER-Co",
    "SER-Fe",
    # "SER-Mn",
    # "SER-Ni",
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
    # "FCC-Co-Co",
    "FCC-Co-Fe",
    # "FCC-Co-Mn",
    # "FCC-Co-Ni",
    # "FCC-Fe-Co",
    # "FCC-Fe-Fe",
    # "FCC-Fe-Mn",
    "FCC-Fe-Ni",
    # "FCC-Mn-Co",
    # "FCC-Mn-Fe",
    # "FCC-Mn-Mn",
    # "FCC-Mn-Ni",
    # "FCC-Ni-Co",
    # "FCC-Ni-Fe",
    # "FCC-Ni-Mn",
    # "FCC-Ni-Ni",
]


def run_tick(name: str, force: bool = False):
    """Run a single structure locally."""
    endmember = Endmember()
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("data/endmembers") / name / "magflow"
    json_path = store_dir / f"{name}-mag.json"

    # Skip if already done
    if not force and json_path.exists():
        log.info(f"Structure {name} already done")
        return

    log.info(f"Structure {name} start")

    if force and store_dir.exists():
        shutil.rmtree(store_dir)

    struct = endmember.get_poscar(name, Path("data/poscars"))

    try:
        worker = StaticWorker(
            vasp_args=VASP_ARGS,
            global_incar=GLOBAL_INCAR,
            j_count=4,
            extend_poscar=(2, 2, 2),
            base_spin=2.0,
        )
        worker.run_flow(
            name=name, 
            structure=struct, 
            flow_dir=flow_dir, 
            store_dir=store_dir, 
            resume=not force,
        )
        output = worker.get_result()
        worker.write_result(data=output, json_path=json_path)
    except Exception as e:
        log.error(f"Structure {name} failed: {e}")

    log.info(f"Structure {name} done")


def run_batch(force: bool = False):
    """Run all structures locally."""
    for name in STRUCT_NAMES:
        run_tick(name, force=force)


def submit_jobs(force: bool = False) -> None:
    manager = SlurmJobManager()
    for name in STRUCT_NAMES:
        config = manager.get_cpu_config(
            job_name=f"{name}-mag",
            output_log=f"logs/{name}-mag.log",
            ntasks=32,
            memory="20G",
        )
        job_id = manager.submit_command(
            command=f"python {__file__} --tick {name} {'--force' if force else ''}",
            config=config,
            conda_env="htvasp",
            workdir=".",
        )
        if not job_id:
            log.error(f"Failed to submit job for {name}")
        else:
            log.info(f"Submitted job for {name} with ID {job_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Magnetic workflow")
    parser.add_argument("--tick", type=str, help="Run single structure")
    parser.add_argument("--batch", action="store_true", help="Run all locally")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--force", action="store_true", help="Force re-run")
    args = parser.parse_args()

    if args.tick:
        run_tick(args.tick, force=args.force)
    elif args.batch:
        run_batch(force=args.force)
    elif args.slurm:
        submit_jobs(force=args.force)
    else:
        parser.print_help()

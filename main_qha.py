"""HT-VASP QHA Workflow - Quasi-harmonic approximation."""

import json
import shutil
import logging
import argparse
from datetime import datetime
from pathlib import Path

from htvasp.cli import submit_jobs, run_batch
from htvasp.model import Endmember
from htvasp.workflows import QhaWorker

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    """JSON encoder for datetime objects."""

    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


# VASP configuration
VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
}

GLOBAL_INCAR = {
    "ENCUT": 450,
    "ISTART": 0,
    "ICHARG": 2,
    "ISMEAR": 1,
    "SIGMA": 0.1,
    "ALGO": "Normal",
    "NELM": 200,
    "NELMIN": 6,
    "NELMDL": -6,
    "IBRION": 2,
    "ISIF": 2,
    "NSW": 100,
    "POTIM": 0.2,
    "EDIFF": 1e-6,
    "EDIFFG": -0.01,
    "ISPIN": 2,
    "ISYM": 0,
    "LREAL": "Auto",
    "PREC": "Accurate",
    "SYMPREC": 1e-5,
    "LWAVE": False,
    "LCHARG": False,
    "LORBIT": None,
    "LOPTICS": False,
    "LVTOT": False,
    "KPAR": 4,
    "NCORE": 2,
    "GGA": "PE",
    "AMIX": 0.2,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.8,
    "BMIX_MAG": 1e-4,
}

RELAX_INCAR = {"KPAR": 4, "NCORE": 2}
PHONON_INCAR = {"KPAR": 2, "NCORE": 4}

# Structure names to process
STRUCT_NAMES = [
    "SER-Co",
    "SER-Fe",
]


def run_tick(name: str, force: bool = False) -> None:
    """Run a single structure locally."""
    endmember = Endmember()
    workdir = Path("data/endmembers") / name
    flowdir = workdir / "qhaflow"
    json_path = workdir / f"{name}-qha.json"

    workdir.mkdir(parents=True, exist_ok=True)

    # Skip if already done
    if not force and json_path.exists():
        with open(json_path) as f:
            result = json.load(f)
        if result.get("state") == "successful":
            log.info(f"Structure {name} already done")
            return

    log.info(f"Structure {name} start")

    if force and flowdir.exists():
        shutil.rmtree(flowdir)
    flowdir.mkdir(parents=True, exist_ok=True)

    struct = endmember.get_poscar(name, Path("data/poscars"))

    try:
        worker = QhaWorker(
            worker_name=f"{name}-qha",
            vasp_args=VASP_ARGS,
            global_incar=GLOBAL_INCAR,
            relax_incar=RELAX_INCAR,
            eos_incar=RELAX_INCAR,
            phonon_incar=PHONON_INCAR,
            temperature_range=(0, 3100, 50),
            supercell_matrix=((2, 0, 0), (0, 2, 0), (0, 0, 2)),
        )
        worker.run_flow(name, struct, flowdir, resume=not force)
        qha_data = worker.get_result()

        if qha_data:
            result = {"name": name, "state": "successful", **qha_data}
        else:
            result = {"name": name, "state": "failed", "struct": struct.as_dict()}

    except Exception as e:
        log.error(f"Structure {name} failed: {e}")
        result = {"name": name, "state": "failed", "struct": struct.as_dict()}

    with open(json_path, "w") as f:
        json.dump(result, f, ensure_ascii=False, indent=2, cls=DateTimeEncoder)

    log.info(f"Structure {name} done")


def command_template(name: str, force: bool) -> str:
    """Generate command for Slurm submission."""
    return f"python {__file__} --tick {name} {'--force' if force else ''}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="QHA workflow")
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
            job_name=lambda n: f"{n}-qha",
            output_log=lambda n: f"logs/{n}-qha.log",
            ntasks=48,
            conda_env="htvasp",
        )
    else:
        parser.print_help()

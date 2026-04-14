"""HT-VASP Static Workflow - Static calculation for magnetism."""

import json
import shutil
import logging
import argparse
from datetime import datetime
from pathlib import Path

from htvasp.cli import submit_jobs, run_batch
from htvasp.model import Endmember
from htvasp.workflows import StaticWorker

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
    "KPAR": 2,
    "NCORE": 2,
    "GGA": "PE",
}

# Structure names to process
STRUCT_NAMES = [
    "SER-Co",
    "SER-Fe",
    "SER-Mn",
    "SER-Ni",
    "BCC-Co-Co",
    "BCC-Co-Fe",
    "BCC-Co-Mn",
    "BCC-Co-Ni",
    "BCC-Fe-Fe",
    "BCC-Fe-Mn",
    "BCC-Fe-Ni",
    "BCC-Mn-Mn",
    "BCC-Mn-Ni",
    "BCC-Ni-Ni",
    "FCC-Co-Co",
    "FCC-Co-Fe",
    "FCC-Co-Mn",
    "FCC-Co-Ni",
    "FCC-Fe-Co",
    "FCC-Fe-Fe",
    "FCC-Fe-Mn",
    "FCC-Fe-Ni",
    "FCC-Mn-Co",
    "FCC-Mn-Fe",
    "FCC-Mn-Mn",
    "FCC-Mn-Ni",
    "FCC-Ni-Co",
    "FCC-Ni-Fe",
    "FCC-Ni-Mn",
    "FCC-Ni-Ni",
]


def run_tick(name: str, force: bool = False) -> None:
    """Run a single structure locally."""
    endmember = Endmember()
    workdir = Path("data/endmembers") / name
    flowdir = workdir / "staticflow"
    json_path = workdir / f"{name}-static.json"

    workdir.mkdir(parents=True, exist_ok=True)

    # Skip if already done
    if not force and json_path.exists():
        with open(json_path) as f:
            result = json.load(f)
        if result.get("state") == "successful":
            log.info(f"Structure {name} already done")
            return

    log.info(f"Structure {name} start")

    if flowdir.exists():
        shutil.rmtree(flowdir)
    flowdir.mkdir(parents=True, exist_ok=True)

    struct = endmember.get_poscar(name, Path("data/poscars"))

    try:
        worker = StaticWorker(
            worker_name=f"{name}-static",
            vasp_args=VASP_ARGS,
            global_incar=GLOBAL_INCAR,
        )
        worker.run_flow(name, struct, flowdir, resume=not force)
        static_data = worker.get_result()

        if static_data:
            result = {"name": name, "state": "successful", **static_data}
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
    parser = argparse.ArgumentParser(description="Static workflow")
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
            job_name=lambda n: f"{n}-static",
            output_log=lambda n: f"logs/{n}-static.log",
            ntasks=16,
            conda_env="htvasp",
        )
    else:
        parser.print_help()

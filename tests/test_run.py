import os
import sys
import json
import logging
import argparse
from datetime import datetime
from pathlib import Path

from pymatgen.core import Structure

from htvasp.workflows import RelaxWorker
from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


def run_locally():
    vasp_wrap = (
        "/bin/bash -c '"
        ". /etc/profile.d/modules.sh && "
        "module load vasp-cpu && "
        "srun vasp_std'"
    )
    vasp_gam_wrap = (
        "/bin/bash -c '"
        ". /etc/profile.d/modules.sh && "
        "module load vasp-cpu && "
        "srun vasp_gam'"
    )

    vasp_args = {
        "handlers": [],
        "vasp_cmd": vasp_wrap,
        "vasp_gamma_cmd": vasp_gam_wrap,
    }

    si_structure = Structure(
        lattice=[[0, 2.73, 2.73], [2.73, 0, 2.73], [2.73, 2.73, 0]],
        species=["Si", "Si"],
        coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
    )
    flow_dir = Path("temp", "relax-si")
    json_path = flow_dir.joinpath("si_relax.json")

    worker = RelaxWorker(
        worker_name="relax-si",
        vasp_args=vasp_args,
        incar_settings={"GGA": "PE"},
    )

    output = worker.run_flow(
        struct_name="si",
        structure=si_structure,
        flow_dir=flow_dir,
    )
    if output:
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump(output, jf, indent=2, cls=DateTimeEncoder)


def submit_to_slurm():
    manager = SlurmJobManager()
    config = manager.get_cpu_config(
        tasks_per_node=4,
        memory="4G",
    )
    job_id = manager.submit_command(
        command=f"python {__file__} --local",
        config=config,
        conda_env="ht-vasp",
        workdir=".",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run VASP task")
    parser.add_argument("--local", action="store_true", help="run locally")
    parser.add_argument("--slurm", action="store_true", help="submit to slurm")
    args = parser.parse_args()

    if args.local:
        run_locally()
    elif args.slurm:
        submit_to_slurm()
    else:
        print("specify run mode: --local or --slurm")

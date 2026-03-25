import json
import shutil
import logging
import argparse
from datetime import datetime
from pathlib import Path

from pymatgen.core import Structure

from htvasp.workflows import StaticWorker
from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


vasp_args = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
}

incar_settings = {
    "KPAR": 2,
    "NCORE": 2,
    "GGA": "PE",
}


def run_locally():

    struct = Structure(
        lattice=[[2.73, 0, 0], [0, 2.73, 0], [0, 0, 2.73]],
        species=["Al", "Al"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    flow_dir = Path("temp", "static-Al")
    json_path = flow_dir.joinpath("static_Al.json")
    if flow_dir.exists():
        shutil.rmtree(flow_dir)

    worker = StaticWorker(
        worker_name="static-Al",
        vasp_args=vasp_args,
        incar_settings={
            "GGA": "PE",
        },
    )

    output = worker.run_flow(
        name="Al",
        structure=struct,
        flow_dir=flow_dir,
    )
    if output:
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump(output, jf, indent=2, cls=DateTimeEncoder)


def submit_job():
    manager = SlurmJobManager()
    config = manager.get_cpu_config(
        ntasks=32,
        memory="4G",
    )
    job_id = manager.submit_command(
        command=f"python {__file__} --local",
        config=config,
        conda_env="htvasp",
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
        submit_job()
    else:
        print("specify run mode: --local or --slurm")
